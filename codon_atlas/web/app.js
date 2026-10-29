/* Codon Atlas page logic. */
"use strict";
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const BASES = "TCAG", AA20 = "ACDEFGHIKLMNPQRSTVWY", MODELS = [...AA20, "?"];
const CODONS = []; for (const a of BASES) for (const b of BASES) for (const c of BASES) CODONS.push(a + b + c);
const THREE = { A: "Ala", C: "Cys", D: "Asp", E: "Glu", F: "Phe", G: "Gly", H: "His", I: "Ile", K: "Lys", L: "Leu", M: "Met", N: "Asn", P: "Pro", Q: "Gln", R: "Arg", S: "Ser", T: "Thr", V: "Val", W: "Trp", Y: "Tyr", "*": "Stop", "?": "?" };
const CLS = { A: "np", V: "np", L: "np", I: "np", M: "np", F: "np", W: "np", P: "np", G: "np", S: "po", T: "po", C: "po", Y: "po", N: "po", Q: "po", D: "ac", E: "ac", K: "ba", R: "ba", H: "ba" };
const U = (c) => c.replace(/T/g, "U");
const fmtN = (n) => (n == null ? "—" : Number(n).toLocaleString());
const fmtP = (p) => (p == null ? "—" : p === 0 ? "< 1e-300" : p >= 0.9999 ? p.toFixed(6) : p >= 0.01 ? p.toFixed(4) : p.toExponential(1));

const S = { state: null, file: null, result: null, runId: "", sub: "", codon: null, view: "wheel", pkg: [], batchFiles: [], batchRun: "", poll: null };

async function api(path, opts = {}) {
  const r = await fetch(path, opts);
  const ct = r.headers.get("content-type") || "";
  const body = ct.includes("json") ? await r.json() : await r.text();
  if (!r.ok) throw new Error((body && body.error) || r.statusText);
  return body;
}
const post = (p, d) => api(p, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(d || {}) });
function toast(msg, ms = 2600) { const t = $("#toast"); t.textContent = msg; t.classList.add("on"); clearTimeout(t._h); t._h = setTimeout(() => t.classList.remove("on"), ms); }
function fail(e) { console.error(e); toast(e.message || String(e), 5000); }

/* ---------------- lifecycle: the server quits when the last window closes */
const QS = new URLSearchParams(location.search);  // ?run=<id>&page=code&codon=TGA&shot=1 (shot: no keep-alive, no transitions, for screenshots)
if (QS.get("shot")) document.documentElement.classList.add("shot");
(function hold() { if (QS.get("shot")) return; try { const es = new EventSource("/api/window"); es.onerror = () => {}; } catch (e) {} })();
window.addEventListener("beforeunload", (e) => { if (S.state && S.state.job && S.state.job.active) { e.preventDefault(); e.returnValue = ""; } });

/* ---------------- theme */
function setTheme(t) {
  document.documentElement.setAttribute("data-theme", t);
  try { localStorage.setItem("ca-theme", t); } catch (e) {}
  $$("#themePill button").forEach((b) => b.classList.toggle("on", b.dataset.t === t));
  if (S.result) drawWheel();
}
$$("#themePill button").forEach((b) => b.onclick = () => setTheme(b.dataset.t));
setTheme(document.documentElement.getAttribute("data-theme") || "light");

/* ---------------- navigation */
function go(p) {
  $$("#tabs button").forEach((b) => b.classList.toggle("on", b.dataset.p === p));
  $$("section.page").forEach((s) => s.classList.toggle("on", s.id === "p-" + p));
  if (p === "runs") loadRuns();
  if (p === "db") renderDb();
  try { localStorage.setItem("ca-page", p); } catch (e) {}
  window.scrollTo({ top: 0, behavior: "smooth" });
}
$$("#tabs button").forEach((b) => b.onclick = () => !b.disabled && go(b.dataset.p));
function needResult(on) { $$("#tabs button[data-need]").forEach((b) => b.disabled = !on); }
needResult(false);
$("#jobchip").onclick = () => go(S.state?.job?.kind === "batch" ? "batch" : (["prepare", "download"].includes(S.state?.job?.kind) ? "db" : "analyze"));

/* ---------------- state + job polling */
async function refresh() {
  S.state = await api("/api/state");
  $("#ver").textContent = "v" + S.state.version; $("#footVer").textContent = "version " + S.state.version;
  renderSettings();
  renderJob(S.state.job);
}
function renderSettings() {
  const db = S.state.db, sel = $("#sProfiles"), cur = sel.value;
  sel.innerHTML = Object.entries(db.sets).map(([k, v]) => `<option value="${k}" ${v.ready ? "" : "disabled"}>${esc(v.label)}${v.ready ? ` · ${fmtN(v.families)} families` : " (not prepared)"}</option>`).join("");
  sel.value = cur && db.sets[cur]?.ready ? cur : (db.default || "enone");
  const ts = $("#sTable");
  if (!ts.options.length) { ts.innerHTML = Object.entries(S.state.tables).map(([k, v]) => `<option value="${k}">${k} · ${esc(v)}</option>`).join(""); ts.value = "11"; }
  const ex = $("#sExclude");
  if (!ex.children.length) ex.innerHTML = Object.entries(S.state.classes).map(([k, v]) => `<label class="chk"><input type="checkbox" value="${k}" checked>${esc(v)}</label>`).join("");
  $("#dbNote").innerHTML = db.default ? "" : `No profile database yet · <a href="#" onclick="go('db');return false">set it up</a>`;
  updateRunBtn();
}
function settings() {
  return {
    profiles: $("#sProfiles").value, reference_table: +$("#sTable").value, evalue: +$("#sEvalue").value || 1e-10,
    threshold: +$("#sThr").value || 0.9999, exclude: $$("#sExclude input:checked").map((i) => i.value),
    bootstrap: +$("#sBoot").value, column_cap: +$("#sCap").value || 0.01, min_pp: +$("#sPP").value || 0.95, cpus: +$("#sCpus").value || 0,
  };
}
function updateRunBtn() { $("#runBtn").disabled = !(S.file && S.state?.db?.default && !S.state?.job?.active); $("#bStart").disabled = !(S.batchFiles.length && S.state?.db?.default && !S.state?.job?.active); }

const RIBBON = (() => { let s = ""; for (let i = 0; i < 160; i++) s += CODONS[(Math.random() * 64) | 0] + " "; return s + s; })();
$("#ribbon span").textContent = RIBBON;
let lastJobKey = "";
function renderJob(j) {
  const chip = $("#jobchip");
  chip.classList.toggle("on", !!j.active);
  if (j.active) $("#jobchipText").textContent = `${Math.round(100 * (j.fraction || 0))}% · ${(j.message || "").slice(0, 46)}`;
  const box = $("#progress");
  if (j.kind === "single" && (j.active || j.status)) {
    box.classList.add("on");
    $("#pMsg").textContent = j.status === "error" ? "The analysis failed" : j.status === "stopped" ? "Stopped" : (j.message || "Working");
    $("#pBar").style.width = `${Math.round(100 * (j.fraction || 0))}%`;
    $("#pEta").textContent = j.active && j.eta ? `about ${Math.ceil(j.eta / 60)} min left` : "";
    $("#stopBtn").style.display = j.active ? "" : "none";
    $("#pErr").classList.toggle("on", j.status === "error"); $("#pErr").textContent = j.error || "";
  }
  if (j.kind === "batch" && S.batchRun === j.run_id) loadBatch(j.run_id, true);
  if (["prepare", "download"].includes(j.kind) && $("#p-db").classList.contains("on")) renderDb();
  const key = `${j.run_id}|${j.status}`;
  if (j.status === "done" && key !== lastJobKey && lastJobKey) {
    if (j.kind === "single") { openRun(j.run_id); toast("Analysis finished"); }
    if (j.kind === "batch") toast("Batch finished");
    if (["prepare", "download"].includes(j.kind)) { toast("Database step finished"); refresh(); }
  }
  lastJobKey = key;
  updateRunBtn();
  clearTimeout(S.poll);
  S.poll = setTimeout(async () => { try { const jj = await api("/api/job"); S.state.job = jj; renderJob(jj); } catch (e) {} }, j.active ? 900 : 4000);
}
$("#stopBtn").onclick = () => post("/api/job/stop").then(renderJob).catch(fail);

/* ---------------- choosing a genome */
const drop = $("#drop");
["dragenter", "dragover"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.add("over"); }));
["dragleave", "drop"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
drop.addEventListener("drop", async (e) => {
  const f = e.dataTransfer.files[0]; if (!f) return;
  $(".big", drop).textContent = `Reading ${f.name}…`;
  try { const r = await api(`/api/upload?name=${encodeURIComponent(f.name)}`, { method: "POST", body: f }); setFile(r); }
  catch (err) { fail(err); } finally { $(".big", drop).textContent = "Drop a genome here"; }
});
drop.addEventListener("click", (e) => { if (e.target.id !== "pickBtn") $("#pickBtn").click(); });
$("#pickBtn").onclick = async (e) => {
  e.stopPropagation();
  try { const r = await post("/api/pick", { kind: "file" }); const info = await post("/api/inspect", { path: r.paths[0] }); setFile(info); } catch (err) { if (!/Nothing/.test(err.message)) fail(err); }
};
function setFile(r) {
  S.file = r;
  const i = r.info;
  drop.style.display = "none"; $("#fileinfo").classList.add("on");
  $("#fiName").textContent = r.name;
  $("#fiOrg").innerHTML = `${i.organism ? `<i>${esc(i.organism)}</i> · ` : ""}${i.format === "genbank" ? "GenBank" : "FASTA"}${i.table ? ` · annotated table ${i.table}` : ""}`;
  $("#fiStats").innerHTML = [[fmtN(i.length), "base pairs"], [fmtN(i.contigs), "records"], [(100 * i.gc).toFixed(1) + "%", "GC"],
    [i.genes ? fmtN(i.genes) : "predict", "genes"]].map(([b, s]) => `<div class="stat"><b>${b}</b><span>${s}</span></div>`).join("");
  updateRunBtn();
}
$("#fiClear").onclick = () => { S.file = null; drop.style.display = ""; $("#fileinfo").classList.remove("on"); updateRunBtn(); };
$("#runBtn").onclick = async () => {
  try {
    const j = await post("/api/analyze", { path: S.file.path, settings: settings(), name: S.file.name });
    lastJobKey = `${j.run_id}|running`; renderJob(j); toast("Analysis started");
  } catch (e) { fail(e); }
};

/* ---------------- results */
async function openRun(id, sub = "") {
  try {
    const r = await api(`/api/result?id=${encodeURIComponent(id)}&sub=${sub}`);
    S.result = r; S.runId = id; S.sub = sub; S.codon = null;
    S.pkg = r.hypotheses.map((h) => [h.codon, h.aa]);
    needResult(true); renderResult(); go("code");
  } catch (e) { fail(e); }
}
function byCodon() { const m = {}; for (const c of S.result.inference.codons) m[c.codon] = c; return m; }

function renderResult() {
  const r = S.result, g = r.genome, inf = r.inference;
  $("#rRun").textContent = `${r.meta?.name || g.name} · ${new Date(r.created * 1000).toLocaleString()} · ${r.profiles.set === "enone" ? "Pfam seed (enone)" : "Pfam"} ${r.profiles.release || ""}`;
  $("#rTitle").innerHTML = g.organism ? `<span class="org">${esc(g.organism)}</span>` : esc(g.name);
  $("#rSub").textContent = `${fmtN(g.length)} bp · ${g.contigs} record(s) · GC ${(100 * g.gc).toFixed(1)}% · ${fmtN(g.genes)} genes (${g.gene_source === "annotation" ? "annotated" : "predicted"})`;
  $("#rStats").innerHTML = [[fmtN(inf.n_columns), "aligned columns"], [fmtN(inf.n_families), "Pfam families"], [r.hypotheses.length, "reassignments"], [Math.round(r.seconds) + " s", "run time"]]
    .map(([b, s]) => `<div class="stat"><b>${b}</b><span>${s}</span></div>`).join("");
  const bc = byCodon();
  $("#rCode").innerHTML = CODONS.map((c) => { const x = bc[c]; const ch = esc(x.inferred); return x.status === "reassigned" || x.status === "stop-to-sense" ? `<b title="${c}">${ch}</b>` : x.inferred === "?" ? `<i title="${c}">?</i>` : `<span title="${c}">${ch}</span>`; }).join("");
  $("#rHyp").innerHTML = r.hypotheses.map((h) => `<span class="chip re" data-c="${h.codon}">${U(h.codon)} ${THREE[h.standard]} → ${THREE[h.aa]}</span>`).join("") || `<span class="chip ok">None: matches table ${r.settings.reference_table}</span>`;
  $("#rLeads").innerHTML = r.leads.map((l) => `<span class="chip lead" data-c="${l.codon}">${U(l.codon)} → ${THREE[l.aa]} · 10<sup>${l.log10.toFixed(1)}</sup>, n=${l.n}</span>`).join("") || `<span class="muted small">None</span>`;
  $("#rTflags").innerHTML = r.trna_flags.map((t) => `<span class="chip lead" data-c="${t.reads[0] || ""}">tRNA-${esc(t.amino_acid)}(${t.anticodon}) · ${esc(t.identity[0])}</span>`).join("") || `<span class="muted small">None</span>`;
  $$("#p-code .chip[data-c]").forEach((ch) => ch.onclick = () => ch.dataset.c && selectCodon(ch.dataset.c));
  drawWheel(); drawTable();
  selectCodon(r.hypotheses[0]?.codon || r.leads[0]?.codon || "TGG");
  renderStops(); renderTrnas(); renderSitesForm();
}

function colorFor(x) {
  const css = getComputedStyle(document.documentElement);
  if (x.status === "reassigned" || x.status === "stop-to-sense") return css.getPropertyValue("--warn").trim();
  const aa = x.inferred;
  if (aa === "?" || aa === "*") return css.getPropertyValue("--stop").trim();
  return css.getPropertyValue("--" + CLS[aa]).trim();
}

function drawWheel() {
  if (!S.result) return;
  const bc = byCodon(), R = [70, 128, 186, 232], cx = 260, cy = 260;
  const css = getComputedStyle(document.documentElement), ink = css.getPropertyValue("--ink").trim(), line = css.getPropertyValue("--line2").trim(), panel = css.getPropertyValue("--panel-solid").trim();
  const arc = (r0, r1, a0, a1) => { const p = (r, a) => [cx + r * Math.sin(a), cy - r * Math.cos(a)]; const [x0, y0] = p(r1, a0), [x1, y1] = p(r1, a1), [x2, y2] = p(r0, a1), [x3, y3] = p(r0, a0); const lg = a1 - a0 > Math.PI ? 1 : 0; return `M${x0},${y0}A${r1},${r1} 0 ${lg} 1 ${x1},${y1}L${x2},${y2}A${r0},${r0} 0 ${lg} 0 ${x3},${y3}Z`; };
  const lab = (r, a, t, sz, col, w = 600) => { const x = cx + r * Math.sin(a), y = cy - r * Math.cos(a); return `<text x="${x}" y="${y}" text-anchor="middle" dominant-baseline="central" font-size="${sz}" font-weight="${w}" fill="${col}">${t}</text>`; };
  let svg = `<svg class="wheel" viewBox="0 0 520 520">`;
  const TAU = Math.PI * 2;
  BASES.split("").forEach((a, i) => {
    const a0 = i * TAU / 4, a1 = (i + 1) * TAU / 4;
    svg += `<path d="${arc(26, R[0], a0, a1)}" fill="${panel}" stroke="${line}"/>` + lab((26 + R[0]) / 2, (a0 + a1) / 2, U(a), 20, ink, 700);
    BASES.split("").forEach((b, j) => {
      const b0 = a0 + j * TAU / 16, b1 = b0 + TAU / 16;
      svg += `<path d="${arc(R[0], R[1], b0, b1)}" fill="${panel}" stroke="${line}"/>` + lab((R[0] + R[1]) / 2, (b0 + b1) / 2, U(b), 15, ink, 600);
      BASES.split("").forEach((c, k) => {
        const cod = a + b + c, x = bc[cod], c0 = b0 + k * TAU / 64, c1 = c0 + TAU / 64, mid = (c0 + c1) / 2;
        const col = colorFor(x), sel = S.codon === cod;
        const op = x.inferred === "?" ? 0.28 : 0.2 + 0.7 * Math.min(1, Math.log10(1 + x.n_used) / 3);
        svg += `<g class="seg" data-c="${cod}"><path d="${arc(R[1], R[2], c0, c1)}" fill="${panel}" stroke="${line}"/>` + lab((R[1] + R[2]) / 2 + 4, mid, U(c), 11, ink, 500);
        svg += `<path d="${arc(R[2] + 2, R[3], c0, c1)}" fill="${col}" fill-opacity="${op}" stroke="${sel ? ink : panel}" stroke-width="${sel ? 2.5 : 1}"/>`;
        svg += lab((R[2] + R[3]) / 2 + 1, mid, esc(x.inferred), 13, ink, 700);
        if (x.status === "reassigned" || x.status === "stop-to-sense") svg += lab(R[3] + 14, mid, x.standard === "*" ? "×" : x.standard, 10, col, 700);
        svg += `</g>`;
      });
    });
  });
  svg += `<text x="${cx}" y="${cy}" text-anchor="middle" dominant-baseline="central" font-size="11" fill="${css.getPropertyValue("--muted").trim()}">5′→3′</text></svg>`;
  $("#wheelHost").innerHTML = svg;
  $$("#wheelHost .seg").forEach((g) => g.onclick = () => selectCodon(g.dataset.c));
}
function drawTable() {
  const bc = byCodon(); let h = "";
  for (const a of BASES) for (const c of BASES) for (const b of BASES) {
    const cod = a + b + c, x = bc[cod];
    const cls = x.status === "reassigned" || x.status === "stop-to-sense" ? "re" : x.inferred === "?" ? "un" : "";
    h += `<div class="ct ${cls} ${S.codon === cod ? "sel" : ""}" data-c="${cod}"><b>${U(cod)}</b><span class="aa">${esc(x.inferred)}</span><div class="faint">${x.standard} · n=${fmtN(x.n_used)}</div></div>`;
  }
  $("#tableHost").innerHTML = `<div class="codontable">${h}</div>`;
  $$("#tableHost .ct").forEach((d) => d.onclick = () => selectCodon(d.dataset.c));
}
$$("#viewTabs button").forEach((b) => b.onclick = () => {
  $$("#viewTabs button").forEach((x) => x.classList.toggle("on", x === b));
  $("#wheelHost").style.display = b.dataset.v === "wheel" ? "" : "none"; $("#tableHost").style.display = b.dataset.v === "table" ? "" : "none";
});

function selectCodon(cod) {
  S.codon = cod; drawWheel(); drawTable();
  const x = byCodon()[cod], r = S.result;
  const verdict = {
    standard: `Read as <b>${THREE[x.inferred]}</b>, as in the reference code.`,
    reassigned: `Read as <b style="color:var(--warn)">${THREE[x.inferred]}</b>, not ${THREE[x.standard]} as in table ${r.settings.reference_table}.`,
    "stop-to-sense": `A stop codon in the reference, here read as <b style="color:var(--warn)">${THREE[x.inferred]}</b>.`,
    stop: `Stop codon in the reference; no amino acid passes the threshold (expected).`,
    uninferred: `Not inferred: ${x.n_used ? "evidence is too weak or conflicting" : "no aligned columns"}.`,
  }[x.status];
  const post = x.posterior, order = post.map((p, i) => [p, i]).sort((a, b) => b[0] - a[0]).slice(0, 8);
  const ll = x.loglik, llmax = Math.max(...ll);
  const bars = order.map(([p, i]) => {
    const m = MODELS[i], w = Math.max(0.5, 100 * p);
    const cls = m === x.standard ? "std" : (i === order[0][1] && m !== x.standard && m !== "?" ? "win" : "");
    return `<span class="mono">${m}</span><div class="b"><i class="${cls}" style="width:${w}%"></i></div><span class="mono small">${fmtP(p)}</span>`;
  }).join("");
  const cons = Object.entries(x.consensus || {}), ctot = cons.reduce((s, [, v]) => s + v, 0) || 1;
  const consBars = cons.map(([a, v]) => `<span class="mono">${a}</span><div class="b"><i class="${a === x.standard ? "std" : ""}" style="width:${100 * v / ctot}%"></i></div><span class="mono small">${Math.round(v)}</span>`).join("");
  const fams = (x.families || []).map((f) => `<tr><td class="mono">${esc(f.name)}</td><td>${esc(f.desc)}</td><td>${f.columns}</td></tr>`).join("");
  const trn = (x.trnas || []).map((t) => `<span class="tag ${t.how === "Watson-Crick" ? "a" : ""}">tRNA-${esc(t.amino_acid)}(${t.anticodon}) · ${esc(t.how)}</span>${(t.identity || []).map((s) => `<span class="tag ${/Arg-type|unclass/.test(s) ? "" : "w"}">${esc(s)}</span>`).join("")}`).join("<br>") || `<span class="muted small">No tRNA with a pairing anticodon.</span>`;
  $("#detail").innerHTML = `
    <div class="row"><div class="big-codon">${U(cod)}</div><div><div class="muted small">DNA ${cod} · reference ${THREE[x.standard]}</div><div class="verdict">${verdict}</div></div></div>
    <dl class="kv" style="margin:14px 0">
      <dt>Columns used</dt><dd>${fmtN(x.n_used)} <span class="faint">(${fmtN(x.n_raw)} aligned, ${fmtN(x.n_unique)} unique columns)</span></dd>
      <dt>Best model</dt><dd>${THREE[x.best]} · ${fmtP(x.best_prob)}</dd>
      ${x.log10_alt_vs_standard != null ? `<dt>Best alternative vs ${x.standard}</dt><dd>${THREE[x.alternative]} · 10<sup>${x.log10_alt_vs_standard.toFixed(1)}</sup></dd>` : ""}
      <dt>Bootstrap support</dt><dd>${x.support == null ? "—" : Math.round(100 * x.support) + "% of resamples keep the same best model"}</dd>
      <dt>Use in genes</dt><dd>${fmtN(x.usage)} codons</dd>
    </dl>
    <h3>Decoding probability</h3><div class="bars">${bars}</div>
    ${cons.length ? `<h3 style="margin-top:16px">Most expected residue at its aligned columns</h3><div class="bars">${consBars}</div>` : ""}
    ${fams ? `<h3 style="margin-top:16px">Families contributing most</h3><div class="scroll" style="max-height:200px"><table class="tbl"><tr><th>Family</th><th>Description</th><th>Columns</th></tr>${fams}</table></div>` : ""}
    <h3 style="margin-top:16px">tRNAs that can read it</h3><div>${trn}</div>
    ${x.status !== "standard" && x.status !== "stop" ? `<div class="row" style="margin-top:16px"><button class="btn sm" id="toSites">List ${U(cod)} sites →</button></div>` : ""}`;
  const ts = $("#toSites"); if (ts) ts.onclick = () => { $("#hCodon").value = cod; $("#hAa").value = x.inferred !== "?" ? x.inferred : (x.alternative || x.best); go("sites"); listSites(); };
}

/* ---------------- stops */
function renderStops() {
  const s = S.result.stops, bc = byCodon();
  $("#stopCards").innerHTML = ["TAA", "TAG", "TGA"].map((c) => {
    const x = s.summary[c], inf = bc[c];
    const exp = Object.entries(x.expected).map(([a, n]) => `<span class="tag">${a} ${n}</span>`).join("");
    const kinds = Object.entries(x.kinds).map(([k, n]) => `<div class="small">${esc(k)} <b>${n}</b></div>`).join("");
    return `<div class="card"><div class="row"><div class="big-codon" style="font-size:30px">${U(c)}</div><span class="spacer"></span>
      ${inf.inferred !== "?" ? `<span class="chip re">read as ${THREE[inf.inferred]}</span>` : `<span class="chip">stop</span>`}</div>
      <div class="stats" style="margin:10px 0"><div class="stat"><b>${fmtN(x.events)}</b><span>in-domain</span></div><div class="stat"><b>${fmtN(x.genes)}</b><span>genes</span></div><div class="stat"><b>${fmtN(inf.usage)}</b><span>in genes</span></div></div>
      ${exp ? `<div class="small muted">Columns expect</div><div>${exp}</div>` : ""}<div style="margin-top:8px">${kinds || '<span class="muted small">None found.</span>'}</div></div>`;
  }).join("");
  drawStopTable();
}
function drawStopTable() {
  const f = $("#stopFilter").value, q = $("#stopSearch").value.toLowerCase();
  const ev = S.result.stops.events.filter((e) => (!f || e.codon === f) && (!q || `${e.family} ${e.desc} ${e.gene?.id || ""} ${e.gene?.product || ""}`.toLowerCase().includes(q)));
  $("#stopTbl").innerHTML = `<tr><th>Codon</th><th>Position</th><th>Gene</th><th>Family</th><th>Column expects</th><th>Kind</th><th>E-value</th></tr>` +
    ev.slice(0, 1500).map((e) => `<tr class="${/Seleno|Pyrro/.test(e.kind) ? "flag" : ""}"><td class="mono">${U(e.codon)}</td><td class="mono">${esc(e.contig)}:${fmtN(e.pos)} ${e.strand}</td>
    <td>${e.gene ? `${esc(e.gene.locus_tag || e.gene.id)} <span class="faint">aa ${e.gene.aa_pos}</span><div class="faint">${esc(e.gene.product)}</div>` : '<span class="faint">intergenic</span>'}</td>
    <td><span class="mono">${esc(e.family)}</span><div class="faint">${esc(e.desc)}</div></td><td class="mono">${e.expected} ${e.expected_p.toFixed(2)} · ${e.second} ${e.second_p.toFixed(2)}</td><td>${esc(e.kind)}</td><td class="mono">${e.evalue.toExponential(0)}</td></tr>`).join("") +
    (ev.length ? "" : `<tr><td colspan="7" class="muted">No stop codons on aligned columns.</td></tr>`);
}
$("#stopFilter").onchange = drawStopTable; $("#stopSearch").oninput = drawStopTable;

/* ---------------- tRNAs */
function renderTrnas() {
  const r = S.result;
  const rows = [...r.trnas].sort((a, b) => (b.identity.length - a.identity.length) || a.amino_acid.localeCompare(b.amino_acid));
  $("#trnaTbl").innerHTML = `<tr><th>tRNA</th><th>Anticodon</th><th>Reads</th><th>Location</th><th>N20</th><th>N73</th><th>1:72</th><th>Identity (Arg anticodons)</th><th>Sequence</th></tr>` +
    rows.map((t) => { const flag = t.identity.some((s) => !/Arg-type|unclass/.test(s));
      return `<tr class="${flag ? "flag" : ""}"><td><b>${esc(t.amino_acid)}</b></td><td class="mono">${U(t.anticodon)}</td><td class="mono">${t.reads.map((x) => `<span title="${esc(x.how)}">${U(x.codon)}</span>`).join(" ")}</td>
      <td class="mono small">${esc(t.contig)}:${fmtN(t.begin)}–${fmtN(t.end)} ${t.strand}</td><td class="mono">${t.n20}</td><td class="mono">${t.n73}</td><td class="mono">${t.pair_1_72}</td>
      <td>${t.identity.map((s) => `<span class="tag ${/Arg-type|unclass/.test(s) ? "" : "w"}">${esc(s)}</span>`).join("")}</td><td class="seq">${esc(t.sequence)}</td></tr>`; }).join("") +
    (rows.length ? "" : `<tr><td colspan="9" class="muted">${esc(r.trna_error || "No tRNAs found.")}</td></tr>`);
}

/* ---------------- sites + export */
function renderSitesForm() {
  const bc = byCodon();
  $("#hCodon").innerHTML = CODONS.map((c) => `<option value="${c}">${U(c)} (${THREE[bc[c].standard]})</option>`).join("");
  $("#hAa").innerHTML = [...AA20].map((a) => `<option value="${a}">${THREE[a]} (${a})</option>`).join("");
  const h = S.result.hypotheses[0] || S.result.leads[0];
  if (h) { $("#hCodon").value = h.codon; $("#hAa").value = h.aa; }
  renderPkg();
  $("#siteTbl").innerHTML = ""; $("#siteTitle").textContent = "Sites"; $("#siteNote").textContent = "";
}
function renderPkg() {
  $("#hList").innerHTML = S.pkg.map(([c, a], i) => `<span class="chip re" data-i="${i}" title="Remove">${U(c)} → ${THREE[a]} ✕</span>`).join("") || `<span class="muted small">Add a hypothesis above.</span>`;
  $$("#hList .chip").forEach((ch) => ch.onclick = () => { S.pkg.splice(+ch.dataset.i, 1); renderPkg(); });
  $("#dlPkg").disabled = !S.pkg.length;
}
$("#hAdd").onclick = () => { const c = $("#hCodon").value, a = $("#hAa").value; if (!S.pkg.some(([x]) => x === c)) S.pkg.push([c, a]); else S.pkg = S.pkg.map(([x, y]) => [x, x === c ? a : y]); renderPkg(); };
$("#hShow").onclick = listSites;
async function listSites() {
  const c = $("#hCodon").value, a = $("#hAa").value;
  $("#siteTitle").textContent = `${U(c)} → ${THREE[a]}: loading…`;
  try {
    const r = await api(`/api/sites?id=${S.runId}&sub=${S.sub}&codon=${c}&aa=${a}`);
    $("#siteTitle").textContent = `${U(c)} → ${THREE[a]}: ${fmtN(r.sites.length)} sites in ${fmtN(r.proteins)} proteins`;
    $("#siteNote").textContent = `${r.in_domains} in Pfam-aligned positions · genes ${r.gene_source === "annotation" ? "from the annotation" : "predicted"}`;
    $("#siteTbl").innerHTML = `<tr><th>Protein</th><th>Product</th><th>Position</th><th>Context</th><th>Pfam column</th><th>P(${c === "" ? "" : r.sites[0]?.standard || "std"})</th><th>P(${a})</th><th>log₁₀ ratio</th></tr>` +
      r.sites.map((s) => `<tr class="${s.log10_variant_vs_standard > 0 ? "flag" : ""}"><td class="mono">${esc(s.locus_tag || s.protein_id)}<div class="faint">${esc(s.protein_id)}</div></td><td>${esc(s.product)}</td>
      <td class="mono">aa ${s.aa_pos}<div class="faint">${esc(s.record)}:${fmtN(s.genome_pos)} ${s.strand}</div></td><td class="mono">${esc(s.context).replace(/\[(.)\]/, '<b style="color:var(--warn)">$1</b>')}</td>
      <td>${s.pfam ? `<span class="mono">${esc(s.pfam)}</span> <span class="faint">#${s.column} · expects ${s.column_consensus}</span>` : '<span class="faint">—</span>'}</td>
      <td class="mono">${s.p_standard == null ? "" : s.p_standard.toFixed(3)}</td><td class="mono">${s.p_variant == null ? "" : s.p_variant.toFixed(3)}</td><td class="mono">${s.log10_variant_vs_standard == null ? "" : s.log10_variant_vs_standard.toFixed(2)}</td></tr>`).join("");
  } catch (e) { fail(e); $("#siteTitle").textContent = "Sites"; }
}
function download(url) { const a = document.createElement("a"); a.href = url; a.download = ""; document.body.appendChild(a); a.click(); a.remove(); }
$("#dlPkg").onclick = () => { download(`/api/export/package?id=${S.runId}&sub=${S.sub}&hyp=${encodeURIComponent(S.pkg.map(([c, a]) => c + ":" + a).join(","))}`); toast("Building the database package…"); };
$("#dlTables").onclick = () => download(`/api/export/tables?id=${S.runId}&sub=${S.sub}`);
$("#dlReport").onclick = () => download(`/api/export/report?id=${S.runId}&sub=${S.sub}`);

/* ---------------- batch */
const bdrop = $("#bdrop");
["dragenter", "dragover"].forEach((ev) => bdrop.addEventListener(ev, (e) => { e.preventDefault(); bdrop.classList.add("over"); }));
["dragleave", "drop"].forEach((ev) => bdrop.addEventListener(ev, (e) => { e.preventDefault(); bdrop.classList.remove("over"); }));
bdrop.addEventListener("drop", async (e) => {
  const files = [...e.dataTransfer.files];
  for (const f of files) {
    try { const r = await api(`/api/upload?name=${encodeURIComponent(f.name)}`, { method: "POST", body: f }); S.batchFiles.push({ path: r.path, name: r.name }); }
    catch (err) { toast(`${f.name}: ${err.message}`, 4000); }
    renderQueue();
  }
});
$("#bFiles").onclick = async () => { try { const r = await post("/api/pick", { kind: "files" }); r.paths.forEach((p) => S.batchFiles.push({ path: p, name: p.split("/").pop() })); renderQueue(); } catch (e) { if (!/Nothing/.test(e.message)) fail(e); } };
$("#bFolder").onclick = async () => { try { const r = await post("/api/pick", { kind: "folder" }); const info = await post("/api/inspect", { path: r.paths[0] }); info.files.forEach((f) => S.batchFiles.push(f)); renderQueue(); if (!info.files.length) toast("No genome files in that folder"); } catch (e) { if (!/Nothing/.test(e.message)) fail(e); } };
function renderQueue() {
  $("#bQueue").innerHTML = S.batchFiles.length ? `<b>${S.batchFiles.length} genome(s)</b> <a href="#" id="bClear">clear</a><div class="scroll" style="max-height:160px;margin-top:8px;padding:6px 10px">${S.batchFiles.map((f) => `<div class="mono">${esc(f.name)}</div>`).join("")}</div>` : "No genomes chosen.";
  const c = $("#bClear"); if (c) c.onclick = (e) => { e.preventDefault(); S.batchFiles = []; renderQueue(); };
  updateRunBtn();
}
$("#bStart").onclick = async () => {
  try { const j = await post("/api/batch/start", { paths: S.batchFiles.map((f) => f.path), settings: settings(), name: $("#bName").value }); S.batchRun = j.run_id; S.batchFiles = []; renderQueue(); renderJob(j); toast("Batch started"); }
  catch (e) { fail(e); }
};
async function loadBatch(id, quiet) {
  try {
    const b = await api(`/api/batch?id=${id}`);
    S.batchRun = id;
    $("#bTitle").textContent = b.meta.name || "Results";
    const p = b.progress || {};
    $("#bProg").textContent = p.status === "running" ? `${b.rows.length} of ${b.total} done · ${(p.message || "").slice(0, 70)}` : `${b.rows.length} of ${b.total} genomes · ${p.status || ""}`;
    $("#bExport").disabled = !b.rows.length;
    $("#bTbl").innerHTML = `<tr><th>#</th><th>Genome</th><th>Called</th><th>Leads</th><th>tRNA flags</th><th>In-domain stops</th><th>Inferred code</th><th>Columns</th></tr>` +
      b.rows.map((r) => r.status === "error" ? `<tr><td>${r.index}</td><td>${esc(r.name)}</td><td colspan="6" style="color:var(--warn)">${esc(r.error)}</td></tr>` :
        `<tr class="${r.reassignments.length ? "flag" : ""}" data-sub="${r.folder}" style="cursor:pointer"><td>${r.index}</td><td><b>${esc(r.name)}</b><div class="faint"><i>${esc(r.organism)}</i> · ${fmtN(r.length)} bp · GC ${(100 * r.gc).toFixed(0)}%</div></td>
        <td>${r.reassignments.map((x) => `<span class="tag w">${U(x)}</span>`).join("") || '<span class="faint">—</span>'}</td><td class="small">${r.leads.map((x) => `<div>${U(esc(x))}</div>`).join("") || '<span class="faint">—</span>'}</td>
        <td class="small">${r.trna_flags.map((x) => `<div>${esc(x)}</div>`).join("") || '<span class="faint">—</span>'}</td><td class="mono small">${Object.entries(r.in_domain_stops).map(([k, v]) => `${U(k)} ${v}`).join(" · ")}</td>
        <td class="mono small" style="max-width:260px;word-break:break-all">${esc(r.code)}</td><td class="mono">${fmtN(r.columns)}</td></tr>`).join("") +
      (b.rows.length ? "" : `<tr><td colspan="8" class="muted">Waiting for the first genome…</td></tr>`);
    $$("#bTbl tr[data-sub]").forEach((tr) => tr.onclick = () => openRun(id, tr.dataset.sub));
    if (!quiet) go("batch");
  } catch (e) { if (!quiet) fail(e); }
}
$("#bExport").onclick = () => download(`/api/export/batch?id=${S.batchRun}`);

/* ---------------- runs */
async function loadRuns() {
  try {
    const r = await api("/api/runs");
    $("#runsTbl").innerHTML = `<tr><th>Name</th><th>When</th><th>Kind</th><th>Status</th><th>Summary</th><th></th></tr>` +
      r.runs.map((x) => `<tr><td><b>${esc(x.name)}</b>${x.summary?.organism ? `<div class="faint"><i>${esc(x.summary.organism)}</i></div>` : ""}</td><td class="small">${new Date(x.created * 1000).toLocaleString()}</td>
      <td>${x.kind === "batch" ? "Batch" : "Genome"}</td><td><span class="tag ${x.status === "done" ? "g" : x.status === "error" ? "w" : "a"}">${esc(x.status)}</span>${x.error ? `<div class="faint small">${esc(x.error).slice(0, 120)}</div>` : ""}</td>
      <td class="small">${x.kind === "single" && x.summary ? (x.summary.reassignments.length ? x.summary.reassignments.map((s) => `<span class="tag w">${U(s)}</span>`).join("") : "standard code") + (x.summary.leads.length ? ` · leads ${x.summary.leads.map(U).join(", ")}` : "") : x.summary ? `${x.summary.done}/${x.summary.total} genomes · ${x.summary.hits} with reassignments` : ""}</td>
      <td style="white-space:nowrap"><button class="btn sm" data-open="${x.id}" data-kind="${x.kind}" ${x.status === "done" || x.kind === "batch" ? "" : "disabled"}>Open</button> <button class="btn sm" data-ren="${x.id}">Rename</button> <button class="btn sm danger" data-del="${x.id}">Delete</button></td></tr>`).join("") +
      (r.runs.length ? "" : `<tr><td colspan="6"><div class="empty"><div class="big">No runs yet</div>Analyze a genome to start.</div></td></tr>`);
    $$("[data-open]").forEach((b) => b.onclick = () => b.dataset.kind === "batch" ? loadBatch(b.dataset.open) : openRun(b.dataset.open));
    $$("[data-ren]").forEach((b) => b.onclick = async () => { const n = prompt("New name"); if (n) { await post("/api/runs/rename", { id: b.dataset.ren, name: n }).catch(fail); loadRuns(); } });
    $$("[data-del]").forEach((b) => b.onclick = async () => { if (confirm("Delete this run and its files?")) { await post("/api/runs/delete", { id: b.dataset.del }).catch(fail); loadRuns(); } });
  } catch (e) { fail(e); }
}
$("#revealRuns").onclick = () => post("/api/reveal", {}).catch(fail);

/* ---------------- database */
async function renderDb() {
  try { S.state = await api("/api/state"); } catch (e) { return; }
  const db = S.state.db, j = S.state.job, busy = j.active && ["prepare", "download"].includes(j.kind);
  const src = db.sources;
  const setRow = (k, v, need, needTxt) => `<div class="dbset"><span class="dot ${v.ready ? "ok" : "no"}"></span><div style="flex:1"><b>${esc(v.label)}</b>${k === "enone" ? ' <span class="tag a">recommended</span>' : ""}
    <div class="small muted">${v.ready ? `${fmtN(v.families)} families · Pfam ${esc(v.release || "")} · prepared ${esc(v.built)} · ${fmtN(v.size_mb)} MB` : needTxt}</div></div>
    <button class="btn ${v.ready ? "" : "primary"} sm" data-prep="${k}" ${need || j.active ? "disabled" : ""}>${v.ready ? "Rebuild" : "Prepare"}</button></div>`;
  const fileRow = (k, nm, s, mb) => `<div class="dbset"><span class="dot ${s.path ? "ok" : "no"}"></span><div style="flex:1"><b>${nm}</b><div class="small muted">${s.path ? `${fmtN(s.mb)} MB · ${esc(s.path)}` : `Not downloaded (about ${mb} MB from EBI)`}</div></div>
    <button class="btn sm" data-dl="${k}" ${j.active ? "disabled" : ""}>${s.path ? "Download again" : "Download"}</button></div>`;
  $("#dbCard").innerHTML = `<h2>Source files</h2>${fileRow("seed", "Pfam-A.seed (seed alignments)", src.seed, 195)}${fileRow("hmm", "Pfam-A.hmm (profile HMMs)", src.hmm, 420)}
    <h2 style="margin-top:18px">Profile sets</h2>${setRow("enone", db.sets.enone, !src.seed.path, "Needs Pfam-A.seed. Building takes about 5–10 minutes.")}${setRow("pfam", db.sets.pfam, !src.hmm.path, "Needs Pfam-A.hmm. Preparing takes a few minutes.")}
    ${busy ? `<div style="margin-top:14px"><b>${esc(j.message || "")}</b><div class="track" style="margin-top:8px"><i style="width:${Math.round(100 * (j.fraction || 0))}%"></i></div><div class="row" style="margin-top:8px"><span class="muted small">${j.eta ? `about ${Math.ceil(j.eta / 60)} min left` : ""}</span><span class="spacer"></span><button class="btn sm danger" id="dbStop">Stop</button></div></div>` : ""}
    ${j.status === "error" && ["prepare", "download"].includes(j.kind) ? `<div class="err on">${esc(j.error)}</div>` : ""}
    <p class="small muted" style="margin-top:14px">Data folder: ${esc(db.home)}</p>`;
  $$("[data-prep]").forEach((b) => b.onclick = () => post("/api/db/prepare", { which: b.dataset.prep }).then((x) => { renderJob(x); renderDb(); }).catch(fail));
  $$("[data-dl]").forEach((b) => b.onclick = () => post("/api/db/download", { which: b.dataset.dl }).then((x) => { renderJob(x); renderDb(); }).catch(fail));
  const st = $("#dbStop"); if (st) st.onclick = () => post("/api/job/stop").then(() => renderDb()).catch(fail);
}

/* ---------------- start */
(async function init() {
  try { await refresh(); } catch (e) { fail(e); return; }
  const j = S.state.job;
  if (j.kind === "batch" && j.run_id) S.batchRun = j.run_id;
  if (!S.state.db.default) { go("db"); return; }
  if (QS.get("run")) {
    await openRun(QS.get("run"), QS.get("sub") || "");
    if (QS.get("codon")) selectCodon(QS.get("codon"));
    if (QS.get("page")) go(QS.get("page"));
    return;
  }
  try {
    const r = await api("/api/runs");
    const last = r.runs.find((x) => x.kind === "single" && x.status === "done");
    if (last) { const p = localStorage.getItem("ca-page"); await openRun(last.id); if (p && !["code"].includes(p)) go(p); }
  } catch (e) {}
})();
