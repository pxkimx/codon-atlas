"""Background jobs, run in their own process group so closing the app can stop them cleanly.

    python -m codon_atlas.worker <job-folder>

The folder holds job.json ({"kind": "single"|"batch"|"prepare"|"download", ...}); the worker
writes progress.json as it goes and status "done" or "error" at the end.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time
import traceback
import urllib.request
from pathlib import Path

from . import pfamdb, pipeline

PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"


class Progress:
    def __init__(self, folder: Path):
        self.path = folder / "progress.json"
        self.state = {"status": "running", "message": "Starting", "fraction": 0.0, "started": time.time(),
                      "pid": os.getpid(), "items": []}
        self._last = 0.0
        self.write()

    def __call__(self, message: str, fraction: float, eta=None, force=False):
        self.state.update(message=message, fraction=float(max(0, min(1, fraction))), eta=eta, updated=time.time())
        now = time.time()
        if force or now - self._last > 0.4:
            self._last = now
            self.write()

    def write(self):
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.state))
        os.replace(tmp, self.path)


def run(folder: Path) -> int:
    job = json.loads((folder / "job.json").read_text())
    prog = Progress(folder)
    try:
        kind = job["kind"]
        if kind == "single":
            pipeline.analyze(job["genome"], job.get("settings", {}), folder, prog)
        elif kind == "batch":
            _batch(job, folder, prog)
        elif kind == "prepare":
            pfamdb.prepare(job["which"], prog)
        elif kind == "download":
            _download(job["which"], prog)
        else:
            raise ValueError(f"unknown job kind {kind}")
        prog.state["status"] = "done"
        prog("Done", 1.0, force=True)
        return 0
    except Exception as e:
        prog.state.update(status="error", error=f"{type(e).__name__}: {e}", trace=traceback.format_exc())
        prog.write()
        return 1


def _batch(job: dict, folder: Path, prog: Progress):
    genomes = job["genomes"]
    rows = []
    prof = pfamdb.ProfileSet(job.get("settings", {}).get("profiles", "enone"))
    n = len(genomes)
    prog.state["items"] = [{"name": Path(g).name, "status": "waiting"} for g in genomes]
    for i, g in enumerate(genomes):
        item = prog.state["items"][i]
        item["status"] = "running"
        sub = folder / "genomes" / f"{i + 1:04d}"
        t = time.time()

        def say(m, f, eta=None, **k):
            prog(f"Genome {i + 1} of {n} · {Path(g).name}: {m}", (i + f) / n, eta=None)

        try:
            res = pipeline.analyze(g, job.get("settings", {}), sub, say, prof)
            row = pipeline.summary_row(res)
            row.update(index=i + 1, folder=sub.name, status="done")
            item.update(status="done", seconds=round(time.time() - t, 1))
        except Exception as e:
            row = {"index": i + 1, "name": Path(g).name, "folder": sub.name, "status": "error",
                   "error": f"{type(e).__name__}: {e}", "score": -1}
            item.update(status="error", error=row["error"])
        rows.append(row)
        (folder / "batch.json").write_text(json.dumps({"rows": rows, "total": n}, default=str))
        prog.write()


def _download(which: str, prog: Progress):
    """Fetch a Pfam file from EBI into the data folder (resumable)."""
    names = {"hmm": "Pfam-A.hmm.gz", "seed": "Pfam-A.seed.gz", "version": "Pfam.version.gz"}
    pfamdb.PFAM_DIR.mkdir(parents=True, exist_ok=True)
    for nm in ([names[which], names["version"]] if which in ("hmm", "seed") else [names[which]]):
        dest = pfamdb.PFAM_DIR / nm
        part = dest.with_suffix(dest.suffix + ".part")
        have = part.stat().st_size if part.exists() else 0
        req = urllib.request.Request(PFAM_URL + nm, headers={"Range": f"bytes={have}-"} if have else {})
        with urllib.request.urlopen(req, timeout=60) as r:
            total = int(r.headers.get("Content-Length") or 0) + have
            mode = "ab" if have and r.status == 206 else "wb"
            done = have if mode == "ab" else 0
            t0, d0 = time.time(), done
            with open(part, mode) as fo:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fo.write(chunk)
                    done += len(chunk)
                    rate = (done - d0) / max(time.time() - t0, 1e-3)
                    eta = (total - done) / rate if rate > 0 and total else None
                    prog(f"Downloading {nm} · {done / 1e6:,.0f} of {total / 1e6:,.0f} MB", done / total if total else 0, eta=eta)
        part.rename(dest)
        if nm == names["version"]:
            import gzip
            (pfamdb.PFAM_DIR / "Pfam.version").write_bytes(gzip.decompress(dest.read_bytes()))
            dest.unlink()


if __name__ == "__main__":
    sys.exit(run(Path(sys.argv[1])))
