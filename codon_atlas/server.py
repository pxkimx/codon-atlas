"""Local web server: the page, a JSON API, and one background job at a time. Standard library only."""
from __future__ import annotations

import io
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np

from . import __version__, pfamdb, report, sites
from .gencode import CODONS, TABLES, THREE
from .genome import GENOME_SUFFIXES, read_genome
from .lifecycle import Lifecycle

WEB = Path(__file__).resolve().parent / "web"
PKG_PARENT = Path(__file__).resolve().parent.parent
DEFAULT_PORT = 8784
MAX_UPLOAD = 2 * 1024 ** 3
RUNS = pfamdb.HOME / "runs"
JOBS = pfamdb.HOME / "jobs"
UPLOADS = pfamdb.HOME / "uploads"
TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript", ".css": "text/css",
         ".svg": "image/svg+xml", ".png": "image/png", ".ico": "image/x-icon"}
LIFE: Lifecycle | None = None


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    return o


# ---------------------------------------------------------------- jobs
class Jobs:
    """At most one background process (analysis, batch, database preparation or download)."""

    def __init__(self):
        self.lock = threading.Lock()
        self.proc: subprocess.Popen | None = None
        self.folder: Path | None = None
        self.kind = ""
        self.run_id = ""

    def busy(self) -> bool:
        return self.proc is not None and self.proc.poll() is None

    def start(self, kind: str, job: dict, folder: Path, run_id: str = "") -> dict:
        with self.lock:
            if self.busy():
                raise ValueError("Another job is running. Wait for it to finish or stop it first.")
            folder.mkdir(parents=True, exist_ok=True)
            (folder / "job.json").write_text(json.dumps({"kind": kind, **job}))
            for f in ("progress.json",):
                (folder / f).unlink(missing_ok=True)
            log = open(folder / "worker.log", "ab")
            self.proc = subprocess.Popen([sys.executable, "-m", "codon_atlas.worker", str(folder)], cwd=str(PKG_PARENT),
                                         stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            self.folder, self.kind, self.run_id = folder, kind, run_id
        return self.status()

    def status(self) -> dict:
        if not self.folder:
            return {"active": False}
        p = self.folder / "progress.json"
        prog = {}
        try:
            prog = json.loads(p.read_text()) if p.exists() else {}
        except Exception:
            pass
        running = self.busy()
        if not running and prog.get("status") == "running":
            prog["status"] = "error"
            prog["error"] = prog.get("error") or "The analysis stopped unexpectedly. See worker.log in the run folder."
        return {"active": running, "kind": self.kind, "run_id": self.run_id, **prog}

    def stop(self) -> dict:
        with self.lock:
            if self.busy():
                try:
                    os.killpg(self.proc.pid, signal.SIGTERM)
                except OSError:
                    pass
                try:
                    self.proc.wait(5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(self.proc.pid, signal.SIGKILL)
                    except OSError:
                        pass
                if self.folder:
                    p = self.folder / "progress.json"
                    try:
                        st = json.loads(p.read_text())
                        st.update(status="stopped", message="Stopped")
                        p.write_text(json.dumps(st))
                    except Exception:
                        pass
        return self.status()

    def kill_all(self):
        if self.proc and self.proc.poll() is None:
            try:
                os.killpg(self.proc.pid, signal.SIGKILL)
            except OSError:
                pass


JOBS_MGR = Jobs()


# ---------------------------------------------------------------- runs
def _new_id(prefix: str) -> str:
    return time.strftime("%Y%m%d-%H%M%S") + "-" + prefix + "-" + uuid.uuid4().hex[:4]


def run_dir(run_id: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", run_id or ""):
        raise ValueError("Unknown run.")
    d = RUNS / run_id
    if not d.is_dir():
        raise ValueError("That run no longer exists.")
    return d


def _meta(d: Path) -> dict:
    try:
        return json.loads((d / "meta.json").read_text())
    except Exception:
        return {}


def list_runs() -> dict:
    out = []
    if RUNS.exists():
        for d in sorted(RUNS.iterdir(), reverse=True):
            if not d.is_dir():
                continue
            m = _meta(d)
            if not m:
                continue
            prog = {}
            try:
                prog = json.loads((d / "progress.json").read_text())
            except Exception:
                pass
            status = prog.get("status", "unknown")
            if status == "running" and not (JOBS_MGR.busy() and JOBS_MGR.run_id == d.name):
                status = "stopped"
            row = {**m, "status": status, "error": prog.get("error")}
            if m.get("kind") == "single" and (d / "result.json").exists():
                try:
                    r = json.loads((d / "result.json").read_text())
                    row["summary"] = {"organism": r["genome"]["organism"], "length": r["genome"]["length"],
                                      "code": r["inference"]["code"],
                                      "reassignments": [f"{h['codon']}→{h['aa']}" for h in r["hypotheses"]],
                                      "leads": [f"{l['codon']}→{l['aa']}" for l in r["leads"][:3]]}
                except Exception:
                    pass
            if m.get("kind") == "batch" and (d / "batch.json").exists():
                try:
                    b = json.loads((d / "batch.json").read_text())
                    row["summary"] = {"done": len(b["rows"]), "total": b["total"],
                                      "hits": sum(1 for x in b["rows"] if x.get("reassignments"))}
                except Exception:
                    pass
            out.append(row)
    return {"runs": out, "root": str(RUNS)}


def start_single(path: str, settings: dict, name: str = "") -> dict:
    p = Path(path).expanduser()
    if not p.is_file():
        raise ValueError(f"File not found: {path}")
    _check_db(settings)
    read_genome(p)  # fail fast with a clear message
    rid = _new_id("g")
    d = RUNS / rid
    d.mkdir(parents=True)
    meta = {"id": rid, "kind": "single", "name": name or re.sub(r"\.gz$", "", p.name), "file": str(p),
            "created": time.time(), "settings": settings, "note": ""}
    (d / "meta.json").write_text(json.dumps(meta))
    return JOBS_MGR.start("single", {"genome": str(p), "settings": settings}, d, rid)


def start_batch(paths: list[str], settings: dict, name: str = "") -> dict:
    files = []
    for x in paths:
        p = Path(x).expanduser()
        if p.is_dir():
            files += sorted(str(f) for f in p.iterdir() if f.is_file() and _is_genome(f.name))
        elif p.is_file():
            files.append(str(p))
    if not files:
        raise ValueError("No genome files found (FASTA or GenBank, optionally .gz).")
    _check_db(settings)
    rid = _new_id("b")
    d = RUNS / rid
    d.mkdir(parents=True)
    meta = {"id": rid, "kind": "batch", "name": name or f"Batch of {len(files)} genomes", "files": files,
            "created": time.time(), "settings": settings, "note": ""}
    (d / "meta.json").write_text(json.dumps(meta))
    return JOBS_MGR.start("batch", {"genomes": files, "settings": settings}, d, rid)


def _is_genome(name: str) -> bool:
    n = name.lower()
    n = n[:-3] if n.endswith(".gz") else n
    return n.endswith(GENOME_SUFFIXES)


def _check_db(settings: dict):
    which = settings.get("profiles") or pfamdb.status()["default"]
    if not which:
        raise ValueError("No profile database is ready yet. Open the Database page to set one up.")
    settings["profiles"] = which
    if not pfamdb.status()["sets"][which]["ready"]:
        raise ValueError(f"The '{which}' profile set is not ready. Open the Database page.")


def result_path(run_id: str, sub: str = "") -> Path:
    d = run_dir(run_id)
    if sub:
        if not re.fullmatch(r"\d{4}", sub):
            raise ValueError("Unknown genome in batch.")
        d = d / "genomes" / sub
    return d


def load_result(run_id: str, sub: str = "") -> dict:
    d = result_path(run_id, sub)
    f = d / "result.json"
    if not f.exists():
        raise ValueError("This run has no result yet.")
    r = json.loads(f.read_text())
    r["meta"] = _meta(run_dir(run_id))
    r["run_id"], r["sub"] = run_id, sub
    if sub:
        try:
            files = r["meta"].get("files", [])
            r["meta"] = {**r["meta"], "file": files[int(sub) - 1], "name": Path(files[int(sub) - 1]).name}
        except Exception:
            pass
    return r


def load_batch(run_id: str) -> dict:
    d = run_dir(run_id)
    b = json.loads((d / "batch.json").read_text()) if (d / "batch.json").exists() else {"rows": [], "total": 0}
    try:
        prog = json.loads((d / "progress.json").read_text())
    except Exception:
        prog = {}
    b["meta"] = _meta(d)
    b["progress"] = prog
    b["rows"].sort(key=lambda r: -r.get("score", 0))
    return b


_CACHE: dict = {}


def _genome_for(run_id: str, sub: str):
    r = load_result(run_id, sub)
    path = r["meta"].get("file")
    key = (run_id, sub)
    if _CACHE.get("key") != key:
        g = read_genome(path)
        from .genome import ensure_genes
        ensure_genes(g, r["genome"].get("gene_table", 11))
        obs = np.load(result_path(run_id, sub) / "columns.npy")
        prof = pfamdb.ProfileSet(r["profiles"]["set"])
        _CACHE.clear()
        _CACHE.update(key=key, genome=g, obs=obs, prof=prof, result=r)
    return _CACHE["genome"], _CACHE["obs"], _CACHE["prof"], _CACHE["result"]


def site_payload(run_id: str, sub: str, codon: str, aa: str) -> dict:
    g, obs, prof, r = _genome_for(run_id, sub)
    table_id = r["settings"].get("reference_table", 11)
    rows = sites.site_list(g, codon, aa, obs, prof, table_id)
    for x in rows:
        x.pop("gene_i", None)
    n_dom = sum(1 for x in rows if x["pfam"])
    return {"codon": codon, "aa": aa, "sites": rows, "proteins": len({x["protein_id"] or x["locus_tag"] for x in rows}),
            "in_domains": n_dom, "gene_source": g.gene_source}


def package(run_id: str, sub: str, hyp: str) -> tuple[str, bytes]:
    pairs = []
    for part in (hyp or "").split(","):
        m = re.fullmatch(r"\s*([ACGTU]{3})\s*[:>→2]\s*([A-Z])\s*", part.upper())
        if m:
            pairs.append((m.group(1).replace("U", "T"), m.group(2)))
    if not pairs:
        raise ValueError("Choose at least one codon and amino acid, e.g. CGG → W.")
    g, obs, prof, r = _genome_for(run_id, sub)
    note = f"- Codon Atlas run: {r['meta'].get('name', run_id)} (profiles: {r['profiles']['set']}, Pfam {r['profiles']['release']})\n"
    return sites.database_package(g, pairs, r["settings"].get("reference_table", 11), obs=obs, prof=prof, note=note)


def tables_zip(run_id: str, sub: str) -> tuple[str, bytes]:
    r = load_result(run_id, sub)
    buf = io.BytesIO()
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", r["meta"].get("name", "genome"))
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in report.tables(r).items():
            z.writestr(f"{stem}_codon_atlas/{name}", text)
        z.writestr(f"{stem}_codon_atlas/report.html", report.html(r))
        z.writestr(f"{stem}_codon_atlas/result.json", json.dumps(clean({k: v for k, v in r.items() if k != 'meta'}), indent=1))
    return f"{stem}_codon_atlas.zip", buf.getvalue()


def delete_run(run_id: str) -> dict:
    d = run_dir(run_id)
    if JOBS_MGR.busy() and JOBS_MGR.run_id == run_id:
        raise ValueError("This run is still going. Stop it first.")
    shutil.rmtree(d)
    _CACHE.clear()
    return list_runs()


def rename_run(run_id: str, name: str, note: str | None = None) -> dict:
    d = run_dir(run_id)
    m = _meta(d)
    if name:
        m["name"] = name[:200]
    if note is not None:
        m["note"] = note[:4000]
    (d / "meta.json").write_text(json.dumps(m))
    return list_runs()


def save_upload(name: str, body: bytes) -> dict:
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", Path(name).name) or "genome.fasta"
    d = UPLOADS / uuid.uuid4().hex[:10]
    d.mkdir(parents=True, exist_ok=True)
    p = d / name
    p.write_bytes(body)
    try:
        g = read_genome(p)
    except ValueError as e:
        shutil.rmtree(d, ignore_errors=True)
        raise
    return {"path": str(p), "name": name, "info": _genome_info(g)}


def inspect_path(path: str) -> dict:
    p = Path(path).expanduser()
    if p.is_dir():
        files = sorted(f for f in p.iterdir() if f.is_file() and _is_genome(f.name))
        return {"path": str(p), "folder": True, "files": [{"path": str(f), "name": f.name, "mb": round(f.stat().st_size / 1e6, 1)} for f in files]}
    if not p.is_file():
        raise ValueError(f"Not found: {path}")
    g = read_genome(p)
    return {"path": str(p), "name": p.name, "info": _genome_info(g)}


def _genome_info(g) -> dict:
    return {"organism": g.organism, "length": g.length, "contigs": len(g.contigs), "genes": len(g.genes),
            "format": g.source_format, "gc": g.gc(), "table": g.transl_table}


def pick(kind: str) -> dict:
    if sys.platform != "darwin":
        raise ValueError("The file picker is only available on macOS. Drag the file onto the page instead.")
    if kind == "folder":
        script = 'POSIX path of (choose folder with prompt "Choose a folder of genome files")'
    elif kind == "files":
        script = ('set fs to choose file with prompt "Choose genome files (FASTA or GenBank)" with multiple selections allowed\n'
                  'set out to ""\nrepeat with f in fs\nset out to out & POSIX path of f & linefeed\nend repeat\nreturn out')
    else:
        script = 'POSIX path of (choose file with prompt "Choose a genome (FASTA or GenBank)")'
    r = subprocess.run(["osascript", "-e", script], capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise ValueError("Nothing was chosen.")
    paths = [x for x in r.stdout.strip().splitlines() if x.strip()]
    return {"paths": paths}


def state() -> dict:
    db = pfamdb.status()
    return {"version": __version__, "db": db, "job": JOBS_MGR.status(), "tables": {k: v[0] for k, v in TABLES.items()},
            "classes": pfamdb.CLASSES, "cpus": os.cpu_count()}


# ---------------------------------------------------------------- HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = "CodonAtlas"
    protocol_version = "HTTP/1.1"

    def log_message(self, *a):
        pass

    def _send(self, code, body: bytes, ctype="application/json", extra=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(clean(obj), allow_nan=False).encode())

    def _body(self) -> bytes:
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_UPLOAD:
            raise ValueError("File is larger than 2 GB.")
        return self.rfile.read(n) if n else b""

    def _guard(self, fn):
        try:
            self._json(fn())
        except (ValueError, FileNotFoundError) as e:
            self._json({"error": str(e)}, 400)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._json({"error": f"{type(e).__name__}: {e}"}, 500)

    def _file(self, fn):
        try:
            name, data = fn()
        except (ValueError, FileNotFoundError) as e:
            return self._json({"error": str(e)}, 400)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return self._json({"error": f"{type(e).__name__}: {e}"}, 500)
        ctype = "application/zip" if name.endswith(".zip") else ("text/html; charset=utf-8" if name.endswith(".html") else "text/tab-separated-values")
        self._send(200, data, ctype, {"Content-Disposition": f'attachment; filename="{name}"'})

    def do_GET(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        if u.path == "/api/state":
            return self._guard(state)
        if u.path == "/api/window":
            return self._window()
        if u.path == "/api/job":
            return self._guard(JOBS_MGR.status)
        if u.path == "/api/runs":
            return self._guard(list_runs)
        if u.path == "/api/result":
            return self._guard(lambda: load_result(q.get("id", ""), q.get("sub", "")))
        if u.path == "/api/batch":
            return self._guard(lambda: load_batch(q.get("id", "")))
        if u.path == "/api/sites":
            return self._guard(lambda: site_payload(q.get("id", ""), q.get("sub", ""), q.get("codon", ""), q.get("aa", "")))
        if u.path == "/api/export/package":
            return self._file(lambda: package(q.get("id", ""), q.get("sub", ""), q.get("hyp", "")))
        if u.path == "/api/export/tables":
            return self._file(lambda: tables_zip(q.get("id", ""), q.get("sub", "")))
        if u.path == "/api/export/report":
            return self._file(lambda: (re.sub(r"[^A-Za-z0-9_.-]+", "_", load_result(q.get("id", ""), q.get("sub", ""))["meta"].get("name", "genome")) + "_codon_atlas_report.html",
                                       report.html(load_result(q.get("id", ""), q.get("sub", ""))).encode()))
        if u.path == "/api/export/batch":
            return self._file(lambda: report.batch_tsv(load_batch(q.get("id", ""))))
        rel = "index.html" if u.path in ("/", "") else u.path.lstrip("/")
        f = (WEB / rel).resolve()
        if WEB.resolve() not in f.parents or not f.is_file():
            return self._send(404, b"not found", "text/plain")
        self._send(200, f.read_bytes(), TYPES.get(f.suffix, "application/octet-stream"))

    def do_POST(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            body = self._body()
        except ValueError as e:
            return self._json({"error": str(e)}, 413)
        js = lambda: json.loads(body or b"{}")  # noqa: E731
        if u.path == "/api/upload":
            return self._guard(lambda: save_upload(q.get("name", "genome.fasta"), body))
        if u.path == "/api/inspect":
            return self._guard(lambda: inspect_path(js().get("path", "")))
        if u.path == "/api/pick":
            return self._guard(lambda: pick(js().get("kind", "file")))
        if u.path == "/api/analyze":
            d = js()
            return self._guard(lambda: start_single(d.get("path", ""), d.get("settings", {}), d.get("name", "")))
        if u.path == "/api/batch/start":
            d = js()
            return self._guard(lambda: start_batch(d.get("paths", []), d.get("settings", {}), d.get("name", "")))
        if u.path == "/api/job/stop":
            return self._guard(JOBS_MGR.stop)
        if u.path == "/api/runs/delete":
            return self._guard(lambda: delete_run(js().get("id", "")))
        if u.path == "/api/runs/rename":
            d = js()
            return self._guard(lambda: rename_run(d.get("id", ""), d.get("name", ""), d.get("note")))
        if u.path == "/api/db/download":
            d = js()
            which = d.get("which", "")
            if which not in ("hmm", "seed"):
                return self._json({"error": "unknown file"}, 400)
            return self._guard(lambda: JOBS_MGR.start("download", {"which": which}, JOBS / _new_id("dl")))
        if u.path == "/api/db/prepare":
            d = js()
            which = d.get("which", "")
            if which not in pfamdb.SETS:
                return self._json({"error": "unknown profile set"}, 400)
            return self._guard(lambda: JOBS_MGR.start("prepare", {"which": which}, JOBS / _new_id("db")))
        if u.path == "/api/reveal":
            d = js()
            return self._guard(lambda: _reveal(d.get("id", "")))
        if u.path == "/api/quit":
            self._json({"ok": True})
            if LIFE:
                LIFE.quit_now()
            return
        self._send(404, b"not found", "text/plain")

    def _window(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Connection", "close")
        self.end_headers()
        if LIFE:
            LIFE.window_opened()
        try:
            while True:
                self.wfile.write(b": hold\n\n")
                self.wfile.flush()
                threading.Event().wait(2.0)
        except (BrokenPipeError, ConnectionResetError, OSError):
            pass
        finally:
            if LIFE:
                LIFE.window_closed()
        self.close_connection = True


def _reveal(run_id: str) -> dict:
    d = run_dir(run_id) if run_id else RUNS
    subprocess.run(["open", str(d)], check=False)
    return {"ok": True}


def serve(port: int = DEFAULT_PORT, open_browser: bool = True) -> None:
    global LIFE
    for d in (RUNS, JOBS, UPLOADS):
        d.mkdir(parents=True, exist_ok=True)
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    httpd.daemon_threads = True

    def quit_():
        sys.stdout.flush()
        JOBS_MGR.kill_all()  # never leave a search running after the program quits
        os._exit(0)

    LIFE = Lifecycle(on_quit=quit_)
    LIFE.watch()
    url = f"http://localhost:{port}"
    print(f"Codon Atlas {__version__} at {url}", flush=True)
    if open_browser:
        import webbrowser
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        quit_()


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(prog="codon-atlas", description="Codon Atlas local app")
    ap.add_argument("--port", type=int, default=int(os.environ.get("CA_PORT", DEFAULT_PORT)))
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    serve(a.port, not a.no_browser)


if __name__ == "__main__":
    main()
