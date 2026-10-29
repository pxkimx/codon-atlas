"""Quit completely when the user closes the tool.

The page holds a server-sent-events connection open (GET /api/window). A dropped connection is
what counts as "window closed": a reload drops and re-opens it within a second, so the server
waits a few seconds before deciding. The server also quits if the launcher that started it is
gone, and if no window ever connects.
"""
from __future__ import annotations

import os
import threading
import time

GRACE_S = 6.0           # after the last window closes
NEVER_OPENED_S = 180.0  # nobody ever connected
POLL_S = 0.5


class Lifecycle:
    def __init__(self, on_quit=None, grace: float = GRACE_S, never_opened: float = NEVER_OPENED_S):
        self.windows = 0
        self.ever_opened = False
        self.last_zero = time.time()
        self.started = time.time()
        self.grace, self.never_opened = grace, never_opened
        self.lock = threading.Lock()
        self.stopping = False
        self.on_quit = on_quit or (lambda: os._exit(0))
        self.parent = int(os.environ.get("CA_PARENT_PID", "0") or 0)

    def window_opened(self) -> None:
        with self.lock:
            self.windows += 1
            self.ever_opened = True

    def window_closed(self) -> None:
        with self.lock:
            self.windows = max(0, self.windows - 1)
            if self.windows == 0:
                self.last_zero = time.time()

    def _parent_gone(self) -> bool:
        if not self.parent:
            return False
        try:
            os.kill(self.parent, 0)
            return False
        except OSError:
            return True

    def should_quit(self, now: float | None = None) -> bool:
        now = now or time.time()
        with self.lock:
            if self.windows > 0:
                return False
            if self.ever_opened:
                return now - self.last_zero >= self.grace
            return now - self.started >= self.never_opened

    def watch(self) -> threading.Thread:
        def loop():
            while not self.stopping:
                if self._parent_gone() or self.should_quit():
                    self.stopping = True
                    self.on_quit()
                    return
                time.sleep(POLL_S)
        t = threading.Thread(target=loop, daemon=True, name="lifecycle")
        t.start()
        return t

    def quit_now(self) -> None:
        self.stopping = True
        threading.Thread(target=lambda: (time.sleep(0.2), self.on_quit()), daemon=True).start()
