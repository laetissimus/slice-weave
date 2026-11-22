/* SliceWeave - local viewer. No third-party code; talks only to this server. */
(() => {
  "use strict";
  const $ = (s, root = document) => root.querySelector(s);
  const $$ = (s, root = document) => Array.from(root.querySelectorAll(s));
  const el = (tag, attrs = {}, ...kids) => {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (k === "class") e.className = v;
      else if (k === "text") e.textContent = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else if (v !== null && v !== undefined && v !== false) e.setAttribute(k, v === true ? "" : v);
    }
    for (const k of kids) if (k !== null && k !== undefined && k !== false) e.append(k instanceof Node ? k : document.createTextNode(String(k)));
    return e;
  };
  const fmt = (x, d = 2) => (x === null || x === undefined || Number.isNaN(x)) ? "-" : (typeof x === "number" ? x.toFixed(d) : String(x));
  const store = { get(k) { try { return localStorage.getItem("sw." + k); } catch (e) { return null; } },
                  set(k, v) { try { localStorage.setItem("sw." + k, v); } catch (e) { /* storage unavailable */ } } };
  const S = { csrf: null, session: null, cid: null, case: null, poll: null, view: "three", which: "overview",
              grid: null, layers: [], base: null, cmp: null, cursor: [0, 0, 0], tool: "nav", cross: true, showRoi: true,
              report: null, jobKind: null };

  // ---------------------------------------------------------------- API
  async function api(method, path, body) {
    const opt = { method, headers: {}, credentials: "same-origin" };
    if (method !== "GET") opt.headers["X-SW-CSRF"] = S.csrf || "";
    if (body instanceof FormData) opt.body = body;
    else if (body !== undefined) { opt.body = JSON.stringify(body); opt.headers["Content-Type"] = "application/json"; }
    const r = await fetch(path, opt);
    if (r.status === 401) { document.body.textContent = "Session expired: open the login link printed in the SliceWeave window again."; throw new Error("not logged in"); }
    const ct = r.headers.get("content-type") || "";
    const data = ct.includes("json") ? await r.json() : await r.arrayBuffer();
    if (!r.ok) { const e = new Error((data && data.error) || r.statusText); e.data = data; throw e; }
    return { data, headers: r.headers };
  }
  const msg = (id, text, kind) => { const m = $("#" + id); m.textContent = text || ""; m.className = "msg " + (kind || ""); };
  let toastT = null;
  function toast(text, err) {
    const t = $("#toast"); t.textContent = text; t.className = "toast" + (err ? " err" : ""); t.hidden = false;
    clearTimeout(toastT); toastT = setTimeout(() => { t.hidden = true; }, 4000);
  }

  // ---------------------------------------------------------------- views
  function setView(v) {
    S.view = v; store.set("view", v);
    $$("#views button").forEach(b => b.classList.toggle("on", b.dataset.view === v));
    const rep = v === "report";
    $("#reportPage").hidden = !rep; $("#viewer").hidden = rep; $("#toolbar").hidden = rep;
    if (rep) { renderReport(); return; }
    const vw = $("#viewer");
    vw.className = "viewer " + (v === "three" ? "three" : "single");
    $$(".view").forEach(x => x.classList.toggle("only", x.dataset.view === v));
    draw();
  }
  $$("#views button").forEach(b => b.addEventListener("click", () => setView(b.dataset.view)));

  async function init() {
    const s = (await api("GET", "/api/session")).data;
    S.csrf = s.csrf; S.session = s;
    $("#banner").title = s.banner;
    const enc = s.encryption || {};
    const sys = $("#sysInfo"); sys.textContent = "";
    sys.append(el("span", { class: "k", text: "Version" }), el("span", { class: "v", text: `${s.name} ${s.version}` }),
      el("span", { class: "k", text: "Data folder" }), el("span", { class: "v", text: s.data_dir }),
      el("span", { class: "k", text: "Disk encryption" }), el("span", { class: "v " + (enc.status === "encrypted" ? "good" : "warn"), text: `${enc.status} (${enc.detail || ""})` }),
      el("span", { class: "k", text: "Network" }), el("span", { class: "v good", text: s.network_guard && s.network_guard.active ? "offline (guard on)" : "guard off" }),
      el("span", { class: "k", text: "JPEG DICOM" }), el("span", { class: "v", text: s.hooks && s.hooks.pydicom ? `pydicom ${s.hooks.pydicom}` : "not available" }));
    if (s.over_retention && s.over_retention.length) {
      const r = $("#retention"); r.hidden = false;
      r.textContent = `${s.over_retention.length} case(s) older than ${s.retention_days} days: ${s.over_retention.join(", ")}. `;
      r.append(el("button", { class: "danger", text: "Delete them", onclick: async () => { if (!confirm("Delete these cases?")) return; await api("POST", "/api/retention/purge"); r.hidden = true; loadCases(); } }));
    }
    setView(store.get("view") || "three");
    await loadCases();
    const last = store.get("case");
    if (last && $$("#caseSel option").some(o => o.value === last)) await selectCase(last);
  }

  // ---------------------------------------------------------------- cases
  async function loadCases() {
    const list = (await api("GET", "/api/cases")).data.cases;
    const sel = $("#caseSel"); const cur = S.cid;
    sel.textContent = ""; sel.append(el("option", { value: "", text: "No case selected" }));
    for (const c of list) sel.append(el("option", { value: c.case_id, text: `${c.case_id}${c.simulated ? " (simulated)" : ""} - ${c.status}` }));
    if (cur) sel.value = cur;
  }
  $("#caseSel").addEventListener("change", e => selectCase(e.target.value));
  $("#newCase").addEventListener("submit", async ev => {
    ev.preventDefault();
    const cid = ev.target.case_id.value.trim();
    try { await api("POST", "/api/cases", { case_id: cid }); msg("newCaseMsg", "Created.", "ok"); ev.target.reset(); $("#newCaseBox").open = false; await loadCases(); await selectCase(cid); }
    catch (e) { msg("newCaseMsg", e.message, "err"); }
  });
  async function demo(kind) {
    try { const r = (await api("POST", "/api/demo", { kind })).data; await loadCases(); await selectCase(r.case_id); toast("Simulating the demo stacks (about a minute)..."); }
    catch (e) { toast(e.message, true); }
  }
  $("#demoAdult").addEventListener("click", () => demo("adult"));
  $("#demoFetal").addEventListener("click", () => demo("fetal"));
  $("#demoFlair").addEventListener("click", () => demo("flair"));

  async function selectCase(cid) {
    S.cid = cid || null; store.set("case", cid || ""); $("#caseSel").value = cid || "";
    S.grid = null; S.layers = []; S.report = null; S.which = "overview";
    if (!cid) { S.case = null; renderAll(); draw(); return; }
    await refreshCase();
    const job = S.case.job || {};
    if (job.status === "running" || job.status === "queued") startPoll();
    await loadGrid(S.case.result ? "result" : "overview");
    if (S.case.result) loadReport();
  }
  async function refreshCase() {
    if (!S.cid) return;
    S.case = (await api("GET", `/api/cases/${S.cid}`)).data;
    renderAll();
  }
  $("#caseLabel").textContent = "";

  function renderAll() {
    const c = S.case;
    $("#caseLabel").textContent = c ? c.case_id + (c.simulated ? "  (SIMULATED)" : "") : "";
    $("#caseState").textContent = c ? c.status : "";
    $("#caseMeta").textContent = c ? `created ${c.created.slice(0, 16).replace("T", " ")}` : "";
    renderStaged(); renderStacks(); renderPlan(); renderRun(); renderResult();
  }

  // ---------------------------------------------------------------- import
  function renderStaged() {
    const st = (S.case && S.case.staged) || { files: 0, bytes: 0 };
    const ul = $("#staged"); ul.textContent = "";
    if (st.files) ul.append(el("li", {}, el("span", { text: `${st.files} file(s) staged` }), el("span", { text: `${(st.bytes / 1e6).toFixed(1)} MB` })));
    $("#btnImport").disabled = !S.cid || !st.files;
  }
  async function uploadFiles(files) {
    if (!S.cid) { toast("Create or open a case first.", true); return; }
    const arr = Array.from(files); if (!arr.length) return;
    const batch = 200;
    for (let i = 0; i < arr.length; i += batch) {
      const fd = new FormData();
      for (const f of arr.slice(i, i + batch)) fd.append("files", f, f.name);
      msg("importMsg", `Uploading ${Math.min(i + batch, arr.length)} / ${arr.length} files...`);
      try { await api("POST", `/api/cases/${S.cid}/upload`, fd); } catch (e) { msg("importMsg", e.message, "err"); return; }
    }
    msg("importMsg", `${arr.length} file(s) staged.`, "ok"); await refreshCase();
  }
  function dropZone(zone, onFiles) {
    zone.addEventListener("dragover", e => { e.preventDefault(); zone.classList.add("over"); });
    zone.addEventListener("dragleave", () => zone.classList.remove("over"));
    zone.addEventListener("drop", async e => { e.preventDefault(); zone.classList.remove("over"); onFiles(await filesFromDrop(e.dataTransfer)); });
  }
  async function filesFromDrop(dt) {
    const out = [];
    const walk = (entry) => new Promise(res => {
      if (entry.isFile) entry.file(f => { out.push(f); res(); }, () => res());
      else if (entry.isDirectory) {
        const rd = entry.createReader(); const all = [];
        const read = () => rd.readEntries(async ents => { if (!ents.length) { for (const x of all) await walk(x); res(); } else { all.push(...ents); read(); } }, () => res());
        read();
      } else res();
    });
    const items = Array.from(dt.items || []).map(i => i.webkitGetAsEntry && i.webkitGetAsEntry()).filter(Boolean);
    if (items.length) { for (const it of items) await walk(it); return out; }
    return Array.from(dt.files || []);
  }
  dropZone($("#drop"), uploadFiles);
  $("#pickFolder").addEventListener("click", () => $("#folderInput").click());
  $("#pickFiles").addEventListener("click", () => $("#filesInput").click());
  $("#folderInput").addEventListener("change", e => { uploadFiles(e.target.files); e.target.value = ""; });
  $("#filesInput").addEventListener("change", e => { uploadFiles(e.target.files); e.target.value = ""; });
  $("#btnClearStage").addEventListener("click", async () => { if (!S.cid) return; await api("DELETE", `/api/cases/${S.cid}/staging`); msg("importMsg", ""); refreshCase(); });
  $("#btnImport").addEventListener("click", async () => {
    try { await api("POST", `/api/cases/${S.cid}/import`); msg("importMsg", "Importing..."); startPoll(); }
    catch (e) { msg("importMsg", e.message, "err"); }
  });

  // ---------------------------------------------------------------- stacks
  function renderStacks() {
    const box = $("#stackTable"); box.textContent = "";
    const c = S.case;
    const sts = (c && c.stacks) || [];
    $("#stackState").textContent = sts.length ? `${sts.filter(s => s.include).length} of ${sts.length} used` : "";
    if (!sts.length) { box.append(el("p", { class: "muted small", text: "No stacks yet." })); return; }
    const plan = c.plan || {};
    const t = el("table", { class: "stacks" });
    t.append(el("tr", {}, ...["Use", "Ref", "Stack", "Pixel", "Thick.", "Sl.", "Motion"].map(h => el("th", { text: h, title: h === "Motion" ? "Correlation of neighbouring slices: higher = less motion" : "" }))));
    for (const s of sts) {
      const use = el("input", { type: "checkbox" }); use.checked = !!s.include;
      use.addEventListener("change", () => changeStacks([{ id: s.id, include: use.checked }]));
      const ref = el("input", { type: "radio", name: "refStack", disabled: !s.include }); ref.checked = plan.reference === s.id;
      ref.addEventListener("change", () => changePlan({ reference: s.id }));
      const th = el("input", { class: "th", type: "number", step: "0.1", min: "0.3", value: s.thickness_override || s.thickness_mm, title: `Slice thickness (mm). Source: ${s.thickness_source}` });
      th.addEventListener("change", () => changeStacks([{ id: s.id, thickness_override: Number(th.value) === s.thickness_mm ? null : Number(th.value) }]));
      const tags = [];
      if (s.haste_like) tags.push(el("span", { class: "tag on", text: "SSh" , title: "single-shot (HASTE-like) series"}));
      if (s.flair_like) tags.push(el("span", { class: "tag on", text: "FLAIR", title: "2-D FLAIR (multi-shot inversion recovery)" + (s.ti_ms ? `, TI ${Math.round(s.ti_ms)} ms` : "") }));
      if (s.has_mask) tags.push(el("span", { class: "tag ok", text: "mask" }));
      if (s.warnings && s.warnings.length) tags.push(el("span", { class: "tag", text: "!", title: s.warnings.join("; ") }));
      t.append(el("tr", { class: s.include ? "" : "off" },
        el("td", {}, use), el("td", {}, ref),
        el("td", {}, el("span", { class: "lab", text: s.label }), ...tags, el("span", { class: "desc", text: s.description, title: s.description })),
        el("td", { class: "num", text: `${fmt(s.pixel_mm[0], 2)}`,
          title: s.resolution_mm ? `pixel ${fmt(s.pixel_mm[0], 2)} x ${fmt(s.pixel_mm[1], 2)} mm; acquired resolution ${fmt(s.resolution_mm[0], 2)} x ${fmt(s.resolution_mm[1], 2)} mm (${s.resolution_source})` : "" },
          ...(s.resolution_mm && Math.max(...s.resolution_mm) > 1.05 * Math.max(...s.pixel_mm) ? [el("small", { class: "acq", text: ` acq ${fmt(Math.max(...s.resolution_mm), 1)}` })] : [])), el("td", {}, th),
        el("td", { class: "num", text: s.shape[2] }), el("td", { class: "num", text: s.motion_score === undefined ? "-" : fmt(s.motion_score, 2) })));
    }
    box.append(t);
  }
  async function changeStacks(changes) {
    try { await api("PUT", `/api/cases/${S.cid}/stacks`, { changes }); await refreshCase(); await changePlan({}); }
    catch (e) { msg("stackMsg", e.message, "err"); await refreshCase(); }
  }

  // ---------------------------------------------------------------- plan
  function renderPlan() {
    const c = S.case, plan = c && c.plan;
    $$("#modePick button").forEach(b => { b.classList.toggle("on", !!plan && b.dataset.mode === plan.mode); b.disabled = !c || !(c.stacks || []).length; });
    $("#planState").textContent = plan ? plan.roi_source || "" : "";
    $("#modeHint").textContent = !plan ? "" : plan.mode === "fetal"
      ? "Fetal: slice-to-volume motion correction. Place the region (ellipsoid) on the fetal brain in the full view: Draw region, or Centre on cursor. A brain mask imported with the stacks is used if present."
      : plan.mode === "flair"
      ? "2-D FLAIR: the head is found automatically; motion between series is corrected. FLAIR is multi-shot, so slices cannot be corrected one by one: exclude a series degraded by motion. Use the three planes of one FLAIR protocol (same TR, TE, TI)."
      : "Adult / paediatric single-shot T2: the head is found automatically; motion between series is corrected. Choose Slice-to-volume for a moving child.";
    const warn = $("#planWarn"); warn.textContent = "";
    const ws = (plan && plan.warnings) || [];
    ws.forEach(w => warn.append(el("li", { text: w })));
    warn.hidden = !ws.length;
    const sliceOpt = $$("#setMotion option").find(o => o.value === "slice");
    if (sliceOpt) { sliceOpt.disabled = !!plan && plan.mode === "flair"; sliceOpt.title = sliceOpt.disabled ? "Not for multi-shot FLAIR" : ""; }
    const info = $("#roiInfo"); info.textContent = "";
    if (plan) {
      const r = plan.roi;
      info.append(el("span", { class: "k", text: "Region" }), el("span", { class: "v", text: `${r.shape}, ${r.half.map(h => (2 * h).toFixed(0)).join(" x ")} mm` }),
        el("span", { class: "k", text: "Centre (RAS)" }), el("span", { class: "v", text: r.center.map(v => v.toFixed(1)).join(", ") }));
      $("#roiShape").value = r.shape;
      const p = plan.params;
      $("#setRes").value = String(p.resolution.toFixed(1)).replace(/\.0$/, ".0");
      if (!$$("#setRes option").some(o => o.value === $("#setRes").value)) $("#setRes").append(el("option", { text: String(p.resolution) }));
      $("#setRes").value = [...$$("#setRes option")].find(o => Math.abs(Number(o.value) - p.resolution) < 1e-6).value;
      $("#setMotion").value = p.motion; $("#setCycles").value = p.cycles;
      $("#cyclesHint").textContent = cyclesHint(p.motion, p.cycles);
      $("#setAlphaMode").value = p.alpha === null || p.alpha === undefined ? "auto" : "manual";
      $("#alphaField").hidden = $("#setAlphaMode").value === "auto"; if (p.alpha) $("#setAlpha").value = p.alpha;
      $("#setReg").value = p.regularizer; $("#setOutlier").checked = !!p.outlier_rejection; $("#setBias").checked = !!p.bias_correction;
    }
    $$("#secSettings select, #secSettings input, #secPlan button, #roiShape").forEach(x => { if (!x.closest("#modePick")) x.disabled = !plan; });
    $("#setCycles").disabled = !plan || $("#setMotion").value === "none";
    if (!plan) $("#cyclesHint").textContent = "";
  }
  async function changePlan(body) {
    if (!S.cid || !S.case.plan) return;
    try {
      const r = (await api("PUT", `/api/cases/${S.cid}/plan`, body)).data;
      S.case.plan = r.plan; renderGridInfo(r);
      await refreshCase();
      if (body.reference) await loadGrid("overview"); else draw();
    } catch (e) { msg("planMsg", e.message, "err"); }
  }
  function renderGridInfo(r) {
    const g = $("#gridInfo"); g.textContent = "";
    if (!r) return;
    g.append(el("span", { class: "k", text: "Output grid" }), el("span", { class: "v", text: `${r.grid.shape.join(" x ")} at ${r.grid.spacing.toFixed(2)} mm (${r.voxels_m} M voxels)` }));
  }
  $$("#modePick button").forEach(b => b.addEventListener("click", async () => {
    try { const r = (await api("POST", `/api/cases/${S.cid}/plan`, { mode: b.dataset.mode })).data; renderGridInfo(r); await refreshCase(); await loadGrid("overview"); }
    catch (e) { msg("planMsg", e.message, "err"); }
  }));
  $("#btnRoiReset").addEventListener("click", async () => { const m = S.case.plan.mode; const r = (await api("POST", `/api/cases/${S.cid}/plan`, { mode: m })).data; renderGridInfo(r); await refreshCase(); draw(); });
  $("#roiShape").addEventListener("change", e => changePlan({ roi: { shape: e.target.value } }));
  $("#btnRoiCenter").addEventListener("click", () => { if (!S.grid) return; changePlan({ roi: { center: voxToWorld(S.cursor) } }); });
  $("#btnRoiTool").addEventListener("click", () => {
    if (S.which !== "overview") { setWhich("overview"); }
    S.tool = S.tool === "roi" ? "nav" : "roi";
    $("#btnRoiTool").classList.toggle("primary", S.tool === "roi");
    $$(".view canvas").forEach(c => c.classList.toggle("roi", S.tool === "roi"));
    if (S.tool === "roi") toast("Drag a rectangle on a view to set the region in that plane.");
  });
  const params = () => ({ resolution: Number($("#setRes").value), motion: $("#setMotion").value, cycles: Number($("#setCycles").value),
    alpha: $("#setAlphaMode").value === "auto" ? null : Number($("#setAlpha").value || 1), regularizer: $("#setReg").value,
    outlier_rejection: $("#setOutlier").checked, bias_correction: $("#setBias").checked });
  $$("#secSettings select, #secSettings input").forEach(x => x.addEventListener("change", () => {
    $("#alphaField").hidden = $("#setAlphaMode").value === "auto";
    changePlan({ params: params() });
  }));

  function cyclesHint(motion, n) {
    if (motion === "none") return "Cycles: not used without motion correction (intensity correction and the outlier check run once).";
    if (motion === "stack") return "Cycles: rounds of re-registering each series to the latest volume, each followed by intensity correction, " +
      "the outlier check and a new reconstruction. 1 is enough in most cases; 2 if the head moved a lot between series." +
      (n === 0 ? " 0: only the first registration, against the reference stack." : "");
    return "Cycles: rounds of slice-to-volume correction, coarse to fine (1: whole stacks and packets of slices; 2: smaller packets, " +
      "then single slices; 3 and more: single slices, finer each round). Default 6; 8 for strong motion." +
      (n < 4 ? " Below 4, slice motion stays partly uncorrected." : "");
  }

  // ---------------------------------------------------------------- run
  function renderRun() {
    const c = S.case, job = (c && c.job) || { status: "idle" };
    const running = job.status === "running" || job.status === "queued";
    $("#btnRun").disabled = !c || !c.plan || running;
    $("#btnCancel").hidden = !(running && job.kind === "run");
    $("#runState").textContent = running ? `${job.kind} ${job.status} ${job.elapsed ? Math.round(job.elapsed) + " s" : ""}` : (c && c.result ? "done" : "");
    const ul = $("#steps"); ul.textContent = "";
    const steps = (S.session && S.session.steps) || [];
    if (job.kind === "run" && job.steps && Object.keys(job.steps).length) {
      for (const [k, name] of steps) {
        const st = job.steps[k] || {};
        const cls = st.status === "running" ? "active" : st.status === "done" ? "done" : st.status === "error" ? "err" : "";
        ul.append(el("li", { class: cls }, el("span", { text: name }), el("span", { class: "t", text: st.msg || "" })));
      }
    }
    if (job.status === "error") msg("runMsg", job.error, "err");
    else if (job.status === "cancelled") msg("runMsg", "Cancelled.");
    $("#log").textContent = (job.log || []).join("\n");
    $("#log").scrollTop = 1e9;
  }
  $("#btnRun").addEventListener("click", async () => {
    try { await api("POST", `/api/cases/${S.cid}/run`); msg("runMsg", ""); startPoll(); } catch (e) { msg("runMsg", e.message, "err"); }
  });
  $("#btnCancel").addEventListener("click", () => api("POST", `/api/cases/${S.cid}/cancel`));
  function startPoll() {
    clearInterval(S.poll);
    S.poll = setInterval(async () => {
      if (!S.cid) { clearInterval(S.poll); return; }
      let job;
      try { job = (await api("GET", `/api/cases/${S.cid}/job`)).data; } catch (e) { clearInterval(S.poll); return; }
      S.case.job = job; renderRun();
      if (!["running", "queued"].includes(job.status)) {
        clearInterval(S.poll);
        await refreshCase();
        if (job.status === "done" && job.kind === "import") { msg("importMsg", "Imported. Check the stacks, the mode and the region.", "ok"); await loadCases(); await loadGrid("overview"); }
        if (job.status === "done" && job.kind === "run") { toast("Reconstruction finished."); await loadGrid("result"); await loadReport(); }
        if (job.status === "error" && job.kind === "import") msg("importMsg", job.error, "err");
      }
    }, 1000);
  }

  // ---------------------------------------------------------------- result
  function renderResult() {
    const c = S.case, r = c && c.result;
    const box = $("#resInfo"); box.textContent = "";
    const has = !!r;
    ["expNifti", "expAvg", "expDicom", "expDicomStudy", "expBundle", "expCsv", "btnReorient", "btnReorientReset"].forEach(id => { $("#" + id).disabled = !has; });
    $("#whichPick").querySelector('[data-which="result"]').disabled = !has;
    $("#resState").textContent = has ? r.finished.slice(0, 16).replace("T", " ") : "";
    if (!has) { box.append(el("span", { class: "k", text: "-" }), el("span", { class: "v", text: "no reconstruction yet" })); return; }
    const s = r.summary;
    const kv = [["Voxel", `${fmt(s.voxel_mm, 2)} mm isotropic`], ["Grid", s.shape.join(" x ")], ["Stacks", s.stacks],
                ["Slices used", `${s.slices_used} of ${s.slices}`], ["Regularisation", `alpha ${fmt(s.alpha, 2)}`], ["Time", `${Math.round(s.seconds)} s`],
                ["Reoriented", r.reorient ? `${r.reorient.rx}, ${r.reorient.ry}, ${r.reorient.rz} deg` : "no"]];
    for (const [k, v] of kv) box.append(el("span", { class: "k", text: k }), el("span", { class: "v", text: v }));
    if (c.simulated) box.append(el("span", { class: "k bad", text: "Data" }), el("span", { class: "v bad", text: "SIMULATED" }));
    const ro = r.reorient || { rx: 0, ry: 0, rz: 0 };
    for (const a of ["rx", "ry", "rz"]) { $("#" + a).value = ro[a]; $("#" + a + "v").textContent = ro[a]; }
  }
  for (const a of ["rx", "ry", "rz"]) $("#" + a).addEventListener("input", e => { $("#" + a + "v").textContent = e.target.value; });
  async function doReorient(zero) {
    const body = zero ? { rx: 0, ry: 0, rz: 0 } : { rx: Number($("#rx").value), ry: Number($("#ry").value), rz: Number($("#rz").value) };
    msg("reoMsg", "Resampling...");
    try { await api("POST", `/api/cases/${S.cid}/reorient`, body); msg("reoMsg", ""); await refreshCase(); await loadGrid("result"); }
    catch (e) { msg("reoMsg", e.message, "err"); }
  }
  $("#btnReorient").addEventListener("click", () => doReorient(false));
  $("#btnReorientReset").addEventListener("click", () => doReorient(true));
  const dl = kind => { window.location.href = `/api/cases/${S.cid}/download/${kind}`; };
  $("#expNifti").addEventListener("click", () => dl("nifti"));
  $("#expAvg").addEventListener("click", () => dl("average"));
  $("#expDicom").addEventListener("click", () => dl("dicom"));
  $("#expBundle").addEventListener("click", () => dl("bundle"));
  $("#expCsv").addEventListener("click", () => dl("csv"));
  $("#expDicomStudy").addEventListener("click", () => $("#refDicom").click());
  $("#refDicom").addEventListener("change", async e => {
    const f = e.target.files[0]; e.target.value = ""; if (!f) return;
    const fd = new FormData(); fd.append("reference", f, "reference.dcm"); fd.append("which", "srr");
    msg("expMsg", "Writing the DICOM series...");
    try {
      const r = await api("POST", `/api/cases/${S.cid}/dicom-attached`, fd);
      const url = URL.createObjectURL(new Blob([r.data], { type: "application/zip" }));
      const a = el("a", { href: url, download: `${S.cid}_sliceweave_dicom_study.zip` }); document.body.append(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 5000);
      msg("expMsg", "DICOM series attached to the original study (contains identifiers).", "ok");
    } catch (err) { msg("expMsg", err.message, "err"); }
  });

  async function loadReport() {
    try { S.report = (await api("GET", `/api/cases/${S.cid}/report`)).data; } catch (e) { S.report = null; }
    renderQC(); if (S.view === "report") renderReport();
  }
  function renderQC() {
    const box = $("#qc"); box.textContent = "";
    const r = S.report;
    if (!r) { box.append(el("p", { class: "muted small", text: "After reconstruction: per stack, how well each slice matches the volume (NCC; red = excluded) and how far it moved." })); return; }
    const prm = r.params || {};
    // many excluded slices: the patient probably moved during a series; a moving reference stack
    // also distorts the volume the other stacks are registered to
    const nEx = s => s.slices - s.slices_used, frac = s => nEx(s) / Math.max(s.slices, 1);
    const bad = r.stacks.filter(s => nEx(s) >= 3 && frac(s) > 0.1);
    const ref = r.stacks.find(s => s.reference), best = r.stacks.slice().sort((a, b) => frac(a) - frac(b))[0];
    const refTip = ref && best && !best.reference && bad.includes(ref) && frac(ref) - frac(best) > 0.1
      ? ` ${ref.name} is the reference stack (the frame the others are registered to): choose ${best.name} (fewest excluded slices) as reference.` : "";
    if (bad.length && prm.motion !== "slice") {
      const which = bad.map(s => `${s.name}: ${nEx(s)} of ${s.slices}`).join("; ");
      box.append(el("p", { class: "warn small", text: prm.mode === "flair"
        ? `Slices that did not match the volume (${which}). Motion within a FLAIR series cannot be corrected: compare that series with the result and exclude it if it is degraded.${refTip}`
        : `Slices that did not match the volume (${which}): the patient probably moved during the series. Choose Motion correction: Slice-to-volume (6 cycles) and reconstruct again.${refTip}` }));
    } else if (refTip) {
      box.append(el("p", { class: "warn small", text: `${ref.name}: ${nEx(ref)} of ${ref.slices} slices excluded.${refTip} Reconstruct again.` }));
    }
    for (const s of r.stacks) {
      box.append(el("div", { class: "qclabel" }, el("span", { text: s.name + (s.reference ? " (ref)" : "") }), el("span", { text: `${s.slices_used}/${s.slices} used` })));
      const cv = el("canvas", { class: "qcchart" }); box.append(cv);
      requestAnimationFrame(() => drawQC(cv, s.slices_detail));
    }
    box.append(el("p", { class: "legendline", text: "Bars: NCC of each slice with the volume (0-1). Line: displacement from the stack's mean pose (top = 10 mm)." }));
  }
  function drawQC(cv, rows) {
    const W = cv.clientWidth, H = cv.clientHeight; if (!W) return;
    cv.width = W * devicePixelRatio; cv.height = H * devicePixelRatio;
    const ctx = cv.getContext("2d"); ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);
    const n = rows.length, bw = W / Math.max(n, 1);
    rows.forEach((r, i) => {
      const v = Math.max(0, r.ncc || 0);
      ctx.fillStyle = r.included ? "#59c27a" : "#ff6b5e";
      ctx.fillRect(i * bw + 1, H - v * (H - 4), Math.max(1, bw - 2), v * (H - 4));
    });
    ctx.strokeStyle = "#f0bc5a"; ctx.lineWidth = 1.2; ctx.beginPath();
    rows.forEach((r, i) => { const y = H - Math.min(1, (r.displacement_mm || 0) / 10) * (H - 4) - 2; const x = i * bw + bw / 2; i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
    ctx.stroke();
  }
  function renderReport() {
    const d = $("#reportDoc"); d.textContent = "";
    const r = S.report;
    if (!r) { d.append(el("h1", { text: "Report" }), el("p", { class: "muted", text: "No reconstruction yet." })); return; }
    d.append(el("h1", { text: `Reconstruction report - ${r.case_id}` }), el("div", { class: "sub", text: `SliceWeave ${r.version}, ${r.finished}` }));
    if (r.simulated) d.append(el("div", { class: "banner-sim", text: "SIMULATED DATA: every number below is a software test, not a patient result." }));
    const p = r.params;
    d.append(card("Settings", kvTable([["Mode", p.mode], ["Output voxel", `${r.grid.spacing} mm isotropic`], ["Grid", r.grid.shape.join(" x ")],
      ["Motion correction", { none: "none", stack: "between series (rigid, per stack)", slice: "slice-to-volume (stacks, packets, then single slices)" }[p.motion]],
      ["Cycles", p.cycles], ["Regularisation", `${p.regularizer}, alpha ${fmt(r.alpha_used, 2)}${p.alpha === null ? " (automatic from the noise" + (r.alpha_motion && Math.abs(r.alpha_motion - r.alpha_used) > 1e-6 ? `; ${fmt(r.alpha_motion, 2)} for the motion-correction volumes` : "") + ")" : ""}`],
      ["Outlier rejection", p.outlier_rejection ? `on (NCC threshold ${p.ncc_threshold})` : "off"], ["Intensity correction", p.bias_correction ? "on" : "off"],
      ["Region", `${r.roi.shape}, ${r.roi.half.map(h => (2 * h).toFixed(0)).join(" x ")} mm`]]),
      el("p", { class: "how" }, el("b", { text: "How to read this. " }), "The volume is estimated from the acquired slices; detail finer than the acquired in-plane resolution cannot be created, and through-plane detail is recovered from the other orientations. Check the result against the original 2-D series, especially where slices were excluded.")));
    const t = el("table", { class: "t" });
    t.append(el("tr", {}, ...["Stack", "Slices used", "Pixel (mm)", "Thickness (mm)", "Noise (rel.)", "Stack rotation (deg)", "Stack shift (mm)", "Median slice motion (mm)"].map(h => el("th", { text: h }))));
    for (const s of r.stacks) {
      const mot = s.slices_detail.map(x => x.displacement_mm).sort((a, b) => a - b);
      t.append(el("tr", {}, el("td", { text: s.name + (s.reference ? " (reference)" : "") }), el("td", { class: "num", text: `${s.slices_used} / ${s.slices}` }),
        el("td", { class: "num", text: s.pixel_mm.join(" x ") + (s.binned > 1 ? ` (binned ${s.binned})` : "") }), el("td", { class: "num", text: s.thickness_mm }),
        el("td", { class: "num", text: fmt(s.noise_rel, 3) }), el("td", { class: "num", text: fmt(s.stack_rotation_deg, 2) }),
        el("td", { class: "num", text: fmt(s.stack_translation_mm, 2) }), el("td", { class: "num", text: fmt(mot[Math.floor(mot.length / 2)] || 0, 2) })));
    }
    d.append(card("Stacks", t, el("p", { class: "how" }, el("b", { text: "Stack rotation and shift: " }), "rigid correction of each stack relative to the reference (patient motion between series). ",
      el("b", { text: "Slice motion: " }), "how far each slice was moved relative to its stack's mean pose (fetal mode). ",
      el("b", { text: "Slices used: " }), "slices whose match with the volume (NCC) was too low are excluded.")));
    d.append(card("Timing (s)", kvTable(Object.entries(r.timing_s))));
  }
  function kvTable(pairs) { const t = el("table", { class: "t compact" }); for (const [k, v] of pairs) t.append(el("tr", {}, el("td", { class: "muted", text: k }), el("td", { text: String(v) }))); return t; }
  function card(title, ...kids) { return el("div", { class: "card" }, el("h2", { text: title }), ...kids); }

  // ---------------------------------------------------------------- viewer
  function setWhich(w) {
    if (w === "result" && !(S.case && S.case.result)) return;
    loadGrid(w);
  }
  $$("#whichPick button").forEach(b => b.addEventListener("click", () => setWhich(b.dataset.which)));
  async function loadGrid(which) {
    if (!S.cid || !S.case || !(S.case.stacks || []).length) { S.grid = null; S.layers = []; fillSelects(); draw(); return; }
    S.which = which;
    $$("#whichPick button").forEach(b => b.classList.toggle("on", b.dataset.which === which));
    let idx;
    try { idx = (await api("GET", `/api/cases/${S.cid}/grid?which=${which}`)).data; }
    catch (e) { toast(e.message, true); return; }
    const sameGrid = S.grid && S.grid.key === idx.key;
    S.grid = { shape: idx.grid.shape, affine: idx.grid.affine, inv: idx.inv, spacing: idx.grid.spacing, key: idx.key };
    if (!sameGrid) S.cursor = S.grid.shape.map(n => Math.floor(n / 2));
    S.layers = idx.layers.map(l => ({ ...l, data: null }));
    const ref = S.case.plan && S.case.plan.reference;
    if (which === "result") { S.base = "srr"; S.cmp = S.layers.some(l => l.id === ref) ? ref : "average"; }
    else { S.base = S.layers.some(l => l.id === ref) ? ref : (S.layers[0] && S.layers[0].id); S.cmp = null; }
    fillSelects();
    $("#viewerEmpty").hidden = true;
    await ensureData(S.base); draw();
    if (S.cmp) ensureData(S.cmp);
  }
  async function ensureData(id) {
    const l = S.layers.find(x => x.id === id); if (!l || l.data) return;
    if (l.loading) return l.loading;
    l.loading = (async () => {
      const r = await api("GET", `/api/cases/${S.cid}/layer?which=${S.which}&lid=${id}`);
      l.data = new Uint8Array(r.data);
      l.offset = parseFloat(r.headers.get("x-sw-offset") || "0"); l.scale = parseFloat(r.headers.get("x-sw-scale") || "1");
      const toWin = w => [Math.max(0, (w[0] - l.offset) / l.scale), Math.min(255, Math.max((w[1] - l.offset) / l.scale, 2))];
      const wf = (r.headers.get("x-sw-window") || "0,1").split(",").map(Number);
      const wb = (r.headers.get("x-sw-window-brain") || "").split(",").map(Number);
      l.winFull = toWin(wf); l.winBrain = wb.length === 2 && wb.every(Number.isFinite) ? toWin(wb) : l.winFull;
      l.win = (S.winMode === "full" ? l.winFull : l.winBrain).slice();
    })();
    try { await l.loading; } catch (e) { toast(e.message, true); } finally { l.loading = null; }
  }
  $("#winSel").addEventListener("change", () => {
    S.winMode = $("#winSel").value;
    S.layers.forEach(l => { if (l.winFull) l.win = (S.winMode === "full" ? l.winFull : l.winBrain).slice(); });
    draw();
  });
  function fillSelects() {
    const b = $("#baseSel"), c = $("#cmpSel"); b.textContent = ""; c.textContent = "";
    c.append(el("option", { value: "", text: "None" }));
    let g = null, gb = null, gc = null;
    for (const l of S.layers) {
      if (l.group !== g) { g = l.group; gb = el("optgroup", { label: g }); gc = el("optgroup", { label: g }); b.append(gb); c.append(gc); }
      gb.append(el("option", { value: l.id, text: l.title })); gc.append(el("option", { value: l.id, text: l.title }));
    }
    b.disabled = c.disabled = !S.layers.length;
    if (S.base) b.value = S.base; c.value = S.cmp || "";
  }
  $("#baseSel").addEventListener("change", async e => { S.base = e.target.value; await ensureData(S.base); draw(); });
  $("#cmpSel").addEventListener("change", async e => { S.cmp = e.target.value || null; if (S.cmp) await ensureData(S.cmp); });
  $("#optCross").addEventListener("change", e => { S.cross = e.target.checked; draw(); });
  $("#optRoi").addEventListener("change", e => { S.showRoi = e.target.checked; draw(); });
  async function flip() {
    if (!S.cmp) return;
    [S.base, S.cmp] = [S.cmp, S.base]; await ensureData(S.base); fillSelects(); draw();
  }

  // radiological display: patient right on screen left; anterior / superior up
  function sliceDims(view) { const [nx, ny, nz] = S.grid.shape; return view === "axial" ? [nx, ny] : view === "coronal" ? [nx, nz] : [ny, nz]; }
  function screenToVox(view, c, r) {
    const [nx, ny, nz] = S.grid.shape, [i0, j0, k0] = S.cursor;
    if (view === "axial") return [nx - 1 - c, ny - 1 - r, k0];
    if (view === "coronal") return [nx - 1 - c, j0, nz - 1 - r];
    return [i0, ny - 1 - c, nz - 1 - r];
  }
  function voxToScreen(view, i, j, k) {
    const [nx, ny, nz] = S.grid.shape;
    if (view === "axial") return [nx - 1 - i, ny - 1 - j];
    if (view === "coronal") return [nx - 1 - i, nz - 1 - k];
    return [ny - 1 - j, nz - 1 - k];
  }
  const vidx = (i, j, k) => { const [nx, ny] = S.grid.shape; return i + nx * (j + ny * k); };
  const mul = (M, p) => [0, 1, 2].map(r => M[r][0] * p[0] + M[r][1] * p[1] + M[r][2] * p[2] + M[r][3]);
  const voxToWorld = v => mul(S.grid.affine, v);
  const worldToVox = w => mul(S.grid.inv, w);
  function roiVox() {
    // the region in voxel coordinates of the current grid (only when the grid shares its axes)
    const plan = S.case && S.case.plan; if (!plan || !S.grid || S.which !== "overview") return null;
    const r = plan.roi, h = S.grid.spacing;
    return { c: worldToVox(r.center), h: r.half.map(x => x / h), shape: r.shape };
  }

  const offscreen = {};
  function draw() {
    for (const v of $$(".view")) {
      const view = v.dataset.view, cv = $("canvas", v), ctx = cv.getContext("2d");
      const W = cv.clientWidth, H = cv.clientHeight;
      if (!W || !H) continue;
      if (cv.width !== W * devicePixelRatio || cv.height !== H * devicePixelRatio) { cv.width = W * devicePixelRatio; cv.height = H * devicePixelRatio; }
      ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);
      ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
      if (!S.grid) { $(".vlabel", v).textContent = ""; continue; }
      const [w, h] = sliceDims(view);
      const oc = offscreen[view] || (offscreen[view] = document.createElement("canvas"));
      oc.width = w; oc.height = h;
      const octx = oc.getContext("2d"), img = octx.createImageData(w, h), px = img.data;
      const base = S.layers.find(l => l.id === S.base && l.data);
      if (base) {
        const [lo, hi] = base.win, d = base.data, inv = 255 / Math.max(hi - lo, 1e-3);
        for (let r = 0; r < h; r++) for (let c = 0; c < w; c++) {
          const [i, j, k] = screenToVox(view, c, r);
          let g = (d[vidx(i, j, k)] - lo) * inv; g = g < 0 ? 0 : g > 255 ? 255 : g;
          const o = (r * w + c) * 4; px[o] = px[o + 1] = px[o + 2] = g; px[o + 3] = 255;
        }
      }
      octx.putImageData(img, 0, 0);
      const sc = Math.min(W / w, H / h), dw = w * sc, dh = h * sc, dx = (W - dw) / 2, dy = (H - dh) / 2;
      v._map = { sc, dx, dy, w, h };
      ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = "high"; ctx.drawImage(oc, dx, dy, dw, dh);
      const roi = S.showRoi && roiVox();
      if (roi) drawRoi(ctx, view, roi, sc, dx, dy);
      if (S.cross) {
        const [cx, cy] = voxToScreen(view, ...S.cursor);
        ctx.strokeStyle = "rgba(224,160,48,.7)"; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(dx + (cx + .5) * sc, dy); ctx.lineTo(dx + (cx + .5) * sc, dy + dh);
        ctx.moveTo(dx, dy + (cy + .5) * sc); ctx.lineTo(dx + dw, dy + (cy + .5) * sc); ctx.stroke();
      }
      if (v._drag) {
        const [a, b] = v._drag; ctx.strokeStyle = "#7fb8ff"; ctx.setLineDash([4, 3]);
        ctx.strokeRect(dx + Math.min(a[0], b[0]) * sc, dy + Math.min(a[1], b[1]) * sc, Math.abs(b[0] - a[0]) * sc, Math.abs(b[1] - a[1]) * sc); ctx.setLineDash([]);
      }
      const ax = view === "axial" ? 2 : view === "coronal" ? 1 : 0;
      const lt = (S.layers.find(l => l.id === S.base) || {}).title || "";
      $(".vlabel", v).textContent = `${view[0].toUpperCase() + view.slice(1)} ${S.cursor[ax] + 1}/${S.grid.shape[ax]}   ${lt}`;
    }
    renderCursor();
  }
  function drawRoi(ctx, view, roi, sc, dx, dy) {
    const axes = view === "axial" ? [0, 1, 2] : view === "coronal" ? [0, 2, 1] : [1, 2, 0];
    const [a, b, n] = axes;
    const t = (S.cursor[n] - roi.c[n]) / roi.h[n];
    // bounding rectangle in screen space
    const lo = [0, 0, 0], hi = [0, 0, 0];
    for (const q of [a, b]) { lo[q] = roi.c[q] - roi.h[q]; hi[q] = roi.c[q] + roi.h[q]; }
    lo[n] = hi[n] = S.cursor[n];
    const p1 = voxToScreen(view, ...lo), p2 = voxToScreen(view, ...hi);
    const x0 = Math.min(p1[0], p2[0]) + .5, x1 = Math.max(p1[0], p2[0]) + .5, y0 = Math.min(p1[1], p2[1]) + .5, y1 = Math.max(p1[1], p2[1]) + .5;
    const inside = Math.abs(t) <= 1;
    ctx.lineWidth = 1.5; ctx.strokeStyle = inside ? "rgba(127,184,255,.9)" : "rgba(127,184,255,.35)";
    if (roi.shape === "ellipsoid") {
      ctx.setLineDash([4, 4]); ctx.strokeRect(dx + x0 * sc, dy + y0 * sc, (x1 - x0) * sc, (y1 - y0) * sc); ctx.setLineDash([]);
      if (inside) {
        const f = Math.sqrt(1 - t * t);
        ctx.beginPath(); ctx.ellipse(dx + (x0 + x1) / 2 * sc, dy + (y0 + y1) / 2 * sc, (x1 - x0) / 2 * f * sc, (y1 - y0) / 2 * f * sc, 0, 0, 2 * Math.PI); ctx.stroke();
      }
    } else {
      if (!inside) ctx.setLineDash([4, 4]);
      ctx.strokeRect(dx + x0 * sc, dy + y0 * sc, (x1 - x0) * sc, (y1 - y0) * sc); ctx.setLineDash([]);
    }
  }
  function renderCursor() {
    const box = $("#cursor"); box.textContent = ""; const line = $("#cursorLine");
    if (!S.grid) { line.textContent = ""; return; }
    const p = vidx(...S.cursor), wpos = voxToWorld(S.cursor);
    box.append(el("span", { class: "k", text: "voxel" }), el("span", { text: S.cursor.map(x => x + 1).join(", ") }),
      el("span", { class: "k", text: "RAS mm" }), el("span", { text: wpos.map(x => x.toFixed(1)).join(", ") }));
    let lineVal = "";
    for (const l of S.layers) {
      if (!l.data) continue;
      const v = l.offset + l.data[p] * l.scale;
      box.append(el("span", { class: "k", text: l.title }), el("span", { text: v.toFixed(0) }));
      if (l.id === S.base) lineVal = v.toFixed(0);
    }
    line.textContent = `RAS ${wpos.map(x => x.toFixed(1)).join(", ")}   ${lineVal}`;
  }
  function canvasPos(v, ev, clamp) {
    const m = v._map; if (!m) return null;
    const rect = $("canvas", v).getBoundingClientRect();
    let c = Math.floor((ev.clientX - rect.left - m.dx) / m.sc), r = Math.floor((ev.clientY - rect.top - m.dy) / m.sc);
    if (clamp) { c = Math.max(0, Math.min(m.w - 1, c)); r = Math.max(0, Math.min(m.h - 1, r)); }
    else if (c < 0 || r < 0 || c >= m.w || r >= m.h) return null;
    return [c, r];
  }
  for (const v of $$(".view")) {
    const view = v.dataset.view, cv = $("canvas", v);
    let drag = null;
    cv.addEventListener("contextmenu", e => e.preventDefault());
    cv.addEventListener("mousedown", ev => {
      if (!S.grid) return;
      const pos = canvasPos(v, ev); if (!pos) return;
      if (ev.button === 2) { const b = S.layers.find(l => l.id === S.base); drag = { kind: "wl", x: ev.clientX, y: ev.clientY, win: b && b.win ? b.win.slice() : null }; return; }
      if (S.tool === "roi" && S.which === "overview") { v._drag = [pos, pos]; drag = { kind: "roi" }; return; }
      S.cursor = screenToVox(view, ...pos); drag = { kind: "nav" }; draw();
    });
    window.addEventListener("mousemove", ev => {
      if (!drag) return;
      if (drag.kind === "wl" && drag.win) {
        const b = S.layers.find(l => l.id === S.base); const dxm = ev.clientX - drag.x, dym = ev.clientY - drag.y;
        const c = (drag.win[0] + drag.win[1]) / 2 - dym * 0.5, w = Math.max(2, (drag.win[1] - drag.win[0]) + dxm * 0.8);
        b.win = [c - w / 2, c + w / 2]; draw(); return;
      }
      const pos = canvasPos(v, ev, true); if (!pos) return;
      if (drag.kind === "nav") { S.cursor = screenToVox(view, ...pos); draw(); }
      else if (drag.kind === "roi") { v._drag[1] = pos; draw(); }
    });
    window.addEventListener("mouseup", async () => {
      if (drag && drag.kind === "roi" && v._drag) {
        const [p1, p2] = v._drag; v._drag = null;
        if (Math.abs(p1[0] - p2[0]) > 3 && Math.abs(p1[1] - p2[1]) > 3) await setRoiFromRect(view, p1, p2);
        else draw();
      }
      drag = null;
    });
    cv.addEventListener("wheel", ev => {
      if (!S.grid) return; ev.preventDefault();
      const d = ev.deltaY > 0 ? -1 : 1, ax = view === "axial" ? 2 : view === "coronal" ? 1 : 0;
      S.cursor[ax] = Math.max(0, Math.min(S.grid.shape[ax] - 1, S.cursor[ax] + d)); draw();
    }, { passive: false });
  }
  async function setRoiFromRect(view, p1, p2) {
    const roi = roiVox(); if (!roi) return;
    const v1 = screenToVox(view, ...p1), v2 = screenToVox(view, ...p2);
    const axes = view === "axial" ? [0, 1] : view === "coronal" ? [0, 2] : [1, 2];
    const c = roi.c.slice(), hv = roi.h.slice();
    for (const q of axes) { c[q] = (v1[q] + v2[q]) / 2; hv[q] = Math.abs(v1[q] - v2[q]) / 2 + 0.5; }
    S.tool = "nav"; $("#btnRoiTool").classList.remove("primary"); $$(".view canvas").forEach(x => x.classList.remove("roi"));
    await changePlan({ roi: { center: voxToWorld(c), half: hv.map(x => x * S.grid.spacing) } });
  }
  window.addEventListener("resize", () => { draw(); if (S.report) renderQC(); });
  window.addEventListener("keydown", ev => {
    if (!S.grid || ["INPUT", "TEXTAREA", "SELECT"].includes(document.activeElement.tagName)) return;
    if (ev.key === "ArrowUp" || ev.key === "ArrowDown") {
      const ax = S.view === "coronal" ? 1 : S.view === "sagittal" ? 0 : 2;
      S.cursor[ax] = Math.max(0, Math.min(S.grid.shape[ax] - 1, S.cursor[ax] + (ev.key === "ArrowUp" ? 1 : -1))); draw(); ev.preventDefault();
    }
    if (ev.key === "b" || ev.key === "B") flip();
    if (ev.key === "Escape" && S.tool === "roi") { S.tool = "nav"; $("#btnRoiTool").classList.remove("primary"); }
  });

  init().catch(e => toast(e.message, true));
})();
