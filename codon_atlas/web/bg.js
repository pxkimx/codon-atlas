/* Background: mRNA strands drifting across the window, ribosomes translating them, nascent
   peptide chains trailing behind, and tRNAs arriving at the A site. Kept faint so text reads. */
(function () {
  const cv = document.getElementById("bg");
  if (!cv) return;
  const ctx = cv.getContext("2d");
  let W = 0, H = 0, dpr = 1, strands = [], t0 = performance.now(), last = 0, running = true;
  const reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const BASECOL = { A: [52, 199, 89], U: [255, 69, 58], G: [255, 159, 10], C: [10, 132, 255] };
  const AAC = [[142, 107, 224], [31, 158, 137], [217, 67, 47], [47, 111, 217]];

  function dark() { return document.documentElement.getAttribute("data-theme") === "dark"; }
  function rnd(a, b) { return a + Math.random() * (b - a); }
  function seq(n) { let s = ""; for (let i = 0; i < n; i++) s += "AUGC"[(Math.random() * 4) | 0]; return s; }

  function resize() {
    dpr = Math.min(2, window.devicePixelRatio || 1);
    W = window.innerWidth; H = window.innerHeight;
    cv.width = W * dpr; cv.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    build();
  }

  function build() {
    strands = [];
    const n = Math.max(3, Math.round(H / 230));
    for (let i = 0; i < n; i++) {
      const depth = rnd(0.45, 1);
      const s = {
        y: (i + 0.5) * H / n + rnd(-40, 40), amp: rnd(14, 34), k: rnd(0.004, 0.009), ph: rnd(0, 6.28),
        speed: rnd(8, 16) * depth, depth, seq: seq(400), offset: rnd(0, 300), ribos: []
      };
      const nr = Math.random() < 0.7 ? 2 : 1;
      for (let r = 0; r < nr; r++) s.ribos.push({ x: rnd(0.15, 0.85) * W, pep: [], next: 0, trna: null });
      strands.push(s);
    }
  }

  function yAt(s, x, t) { return s.y + s.amp * Math.sin(s.k * x + s.ph + t * 0.25); }

  function frame(now) {
    if (!running) return;
    requestAnimationFrame(frame);
    if (now - last < 33) return;  // ~30 fps
    const dt = Math.min(0.1, (now - last) / 1000); last = now;
    const t = (now - t0) / 1000;
    const dk = dark();
    ctx.clearRect(0, 0, W, H);
    const baseA = dk ? 0.22 : 0.16;
    for (const s of strands) {
      const a = baseA * s.depth;
      const step = 9;  // px per nucleotide
      s.offset += (reduce ? 0 : s.speed * dt);
      // backbone
      ctx.lineWidth = 1.4 * s.depth;
      ctx.strokeStyle = dk ? `rgba(200,210,230,${a})` : `rgba(60,70,90,${a})`;
      ctx.beginPath();
      for (let x = -10; x <= W + 10; x += 8) { const y = yAt(s, x, t); x === -10 ? ctx.moveTo(x, y) : ctx.lineTo(x, y); }
      ctx.stroke();
      // bases as short ticks, grouped in codons
      const first = Math.floor(s.offset / step);
      for (let i = 0; ; i++) {
        const x = (first + i) * step - s.offset;
        if (x > W + 10) break;
        if (x < -10) continue;
        const b = s.seq[(first + i) % s.seq.length];
        const c = BASECOL[b];
        const y = yAt(s, x, t);
        const len = 7 * s.depth * ((first + i) % 3 === 0 ? 1.25 : 1);
        ctx.strokeStyle = `rgba(${c[0]},${c[1]},${c[2]},${a * 1.5})`;
        ctx.lineWidth = 2 * s.depth;
        ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x, y + len); ctx.stroke();
      }
      // ribosomes walk 5'->3' relative to the strand (strand drifts left, ribosomes hold position and translate)
      for (const r of s.ribos) {
        r.x += (reduce ? 0 : rnd(-0.3, 0.3));
        const y = yAt(s, r.x, t);
        const sc = 1.25 * s.depth;
        // small subunit (below) and large subunit (above)
        ctx.fillStyle = dk ? `rgba(170,190,230,${a * 0.9})` : `rgba(90,110,150,${a * 0.75})`;
        ctx.beginPath(); ctx.ellipse(r.x, y + 9 * sc, 20 * sc, 9 * sc, 0, 0, 6.283); ctx.fill();
        ctx.fillStyle = dk ? `rgba(150,175,225,${a * 1.1})` : `rgba(70,95,140,${a * 0.9})`;
        ctx.beginPath(); ctx.ellipse(r.x, y - 12 * sc, 26 * sc, 16 * sc, 0, 0, 6.283); ctx.fill();
        // tRNA delivery
        r.next -= dt;
        if (!r.trna && r.next <= 0 && !reduce) r.trna = { p: 0, col: AAC[(Math.random() * 4) | 0] };
        if (r.trna) {
          r.trna.p += dt / 1.6;
          const p = Math.min(1, r.trna.p), ease = 1 - Math.pow(1 - p, 3);
          const tx = r.x + 8 * sc + (1 - ease) * 60, ty = y - 14 * sc - (1 - ease) * 70;
          ctx.strokeStyle = dk ? `rgba(220,230,255,${a * 1.4})` : `rgba(40,60,100,${a * 1.2})`;
          ctx.lineWidth = 2 * s.depth;
          ctx.beginPath(); ctx.moveTo(tx - 6 * sc, ty - 10 * sc); ctx.lineTo(tx, ty); ctx.lineTo(tx + 7 * sc, ty - 4 * sc); ctx.stroke();
          const c = r.trna.col;
          ctx.fillStyle = `rgba(${c[0]},${c[1]},${c[2]},${a * 2.2})`;
          ctx.beginPath(); ctx.arc(tx - 6 * sc, ty - 12 * sc, 3.2 * sc, 0, 6.283); ctx.fill();
          if (r.trna.p >= 1.25) { r.pep.unshift({ c, w: rnd(0, 6.28) }); if (r.pep.length > 16) r.pep.pop(); r.trna = null; r.next = rnd(0.6, 2.2); }
        }
        // nascent chain leaving the exit tunnel
        let px = r.x - 10 * sc, py = y - 26 * sc;
        for (let j = 0; j < r.pep.length; j++) {
          const q = r.pep[j];
          px += (-5 + 3 * Math.sin(t * 1.3 + q.w)) * sc; py -= 5.5 * sc;
          ctx.fillStyle = `rgba(${q.c[0]},${q.c[1]},${q.c[2]},${a * 2 * (1 - j / 18)})`;
          ctx.beginPath(); ctx.arc(px, py, 3 * sc, 0, 6.283); ctx.fill();
        }
      }
    }
  }

  document.addEventListener("visibilitychange", () => {
    running = !document.hidden;
    if (running) { last = performance.now(); requestAnimationFrame(frame); }
  });
  window.addEventListener("resize", resize);
  resize();
  requestAnimationFrame(frame);
})();
