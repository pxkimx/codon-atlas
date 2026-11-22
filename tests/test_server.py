"""The HTTP server, run for real on a free port with a throw-away data folder (CA_HOME)."""
import json
import os
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from codon_atlas import __version__

ROOT = Path(__file__).resolve().parents[1]
RUN = "20260101-000000-run-abcd"


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture
def server(tmp_path):
    home = tmp_path / "home"
    (home / "runs" / RUN).mkdir(parents=True)
    (home / "runs" / RUN / "meta.json").write_text(json.dumps({"name": "fake run"}))
    (home / "pfam").mkdir()
    (home / "pfam" / "keep.txt").write_text("the database must survive")
    env = {k: v for k, v in os.environ.items() if k not in ("CA_PORT", "CA_PARENT_PID")}
    env["CA_HOME"] = str(home)
    port = free_port()
    p = subprocess.Popen([sys.executable, "-m", "codon_atlas", "--port", str(port), "--no-browser"], cwd=ROOT, env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        for _ in range(150):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/api/state", timeout=1)
                break
            except Exception:  # noqa: BLE001
                time.sleep(0.1)
        else:
            raise RuntimeError("the server did not start")
        yield port, home
    finally:
        p.kill()
        p.wait(5)


def call(port, path, body=None, headers=None):
    data = None if body is None else (body if isinstance(body, bytes) else json.dumps(body).encode())
    req = urllib.request.Request(f"http://127.0.0.1:{port}{path}", data=data,
                                 headers={"Content-Type": "application/json", **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=5) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def runs_of(body: bytes) -> list:
    data = json.loads(body)
    return data["runs"] if isinstance(data, dict) else data


def test_state_from_this_mac(server):
    port, _ = server
    code, body = call(port, "/api/state")
    assert code == 200
    assert json.loads(body)["version"] == __version__
    code, body = call(port, "/api/runs")
    assert code == 200 and [r["name"] for r in runs_of(body)] == ["fake run"]


def test_another_site_cannot_drive_the_app(server):
    port, home = server
    # a page on another site, posting as text/plain so the browser sends it without asking
    code, _ = call(port, "/api/runs/delete", b'{"id": ".."}', {"Origin": "https://evil.example", "Content-Type": "text/plain"})
    assert code == 403
    # a request that reached us through another host name (DNS rebinding)
    code, _ = call(port, "/api/runs", headers={"Host": "evil.example"})
    assert code == 403
    code, _ = call(port, "/api/state", headers={"Origin": "http://localhost:9999"})
    assert code == 200, "the page served from another local port, e.g. inside Toolbox Launcher, is fine"
    assert (home / "runs" / RUN).is_dir() and (home / "pfam" / "keep.txt").exists()


def test_a_run_id_is_one_folder_name(server):
    port, home = server
    for bad in ("..", ".", "...", "../pfam", "runs/..", "", "no-such-run"):
        code, _ = call(port, "/api/runs/delete", {"id": bad})
        assert code == 400, bad
        code, _ = call(port, "/api/reveal", {"id": bad}) if bad else (400, b"")
        assert code == 400, bad
    assert (home / "runs" / RUN).is_dir() and (home / "pfam" / "keep.txt").exists()
    code, body = call(port, "/api/runs/delete", {"id": RUN})
    assert code == 200 and runs_of(body) == []
    assert not (home / "runs" / RUN).exists() and (home / "pfam" / "keep.txt").exists()
