(() => {
  const API = "/doho/research-api";
  const $ = (id) => document.getElementById(id);
  let receptor = null;
  let ligands = [];
  let pollTimer = null;
  let poseState = { report: null, vina: null, viewer: null, receptorText: "", detailById: {}, activeId: null };

  function fileToBase64(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(String(reader.result).split(",")[1] || "");
      reader.onerror = () => reject(reader.error || new Error("文件读取失败"));
      reader.readAsDataURL(file);
    });
  }

  function markReceptor(value) {
    receptor = value;
    $("receptor-name").textContent = value ? `${value.name} · ${(value.bytes / 1024).toFixed(1)} KB` : "未选择文件";
    $("receptor-status").classList.toggle("ok", Boolean(value));
    $("receptor-meta").textContent = value ? "已读取。提交时会检查 ATOM/HETATM 记录，并保存 SHA-256 输入摘要。" : "文件会在提交前做 ATOM/HETATM 记录检查，并写入任务 provenance。";
  }

  function renderLigands() {
    const host = $("ligand-list");
    if (!ligands.length) {
      host.innerHTML = '<div class="empty">还没有配体。可以加载示例，或点击“添加配体”。</div>';
      return;
    }
    host.innerHTML = ligands.map((item, index) => `<div class="ligand-row"><div><b>${escapeHtml(item.name)}</b><small>${(item.bytes / 1024).toFixed(1)} KB · PDBQT · 将独立运行</small></div><span class="engine-badge">待运行</span><button class="remove" data-remove="${index}" aria-label="移除配体">移除</button></div>`).join("");
    host.querySelectorAll("[data-remove]").forEach((button) => button.addEventListener("click", () => { ligands.splice(Number(button.dataset.remove), 1); renderLigands(); }));
  }

  function escapeHtml(value) { return String(value).replace(/[&<>"']/g, (char) => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[char])); }
  function readNumber(id) { return Number($(id).value); }
  function payload() {
    return { receptor: { name: receptor.name, content_base64: receptor.content_base64 }, ligands: ligands.map((item) => ({ id: item.id, name: item.name, content_base64: item.content_base64 })), box: { center: [readNumber("cx"), readNumber("cy"), readNumber("cz")], size: [readNumber("sx"), readNumber("sy"), readNumber("sz")] }, params: { exhaustiveness: readNumber("exhaustiveness"), num_modes: readNumber("num_modes"), energy_range: readNumber("energy_range"), seed: readNumber("seed"), cpu: readNumber("cpu") } };
  }

  function showError(message) { $("run-status").textContent = message; $("run-status").style.color = "#f2a5a5"; }
  function showProgress(text, fraction) {
    $("progress").classList.remove("hidden");
    $("progress").querySelector("i").style.width = `${Math.max(8, Math.min(100, fraction * 100))}%`;
    $("progress").querySelector("p").textContent = text;
  }
  function artifactUrl(runId, path) { return `${API}/runs/${encodeURIComponent(runId)}/artifact/${encodeURIComponent(path).replace(/%2F/g, "/")}`; }
  function resultArtifactPath(runId, detail) {
    const output = detail?.output || `results/vina-runs/${runId}/ligands/${detail?.id || "pose"}.out.pdbqt`;
    return output.replace(/^\/+/, "");
  }
  function pdbqtModel(text, mode) {
    const blocks = [];
    const pattern = /MODEL\s+(\d+)([\s\S]*?)(?:ENDMDL|$)/g;
    let match;
    while ((match = pattern.exec(text))) blocks.push({ number: Number(match[1]), body: match[0] });
    const selected = blocks.find((block) => block.number === mode) || blocks[mode - 1];
    const source = selected ? selected.body : text;
    return source.split(/\r?\n/).filter((line) => /^(ATOM  |HETATM|MODEL |ENDMDL|CONECT|TER|END)/.test(line)).join("\n");
  }
  function clearViewer() {
    if (poseState.viewer) { try { poseState.viewer.clear(); } catch (_) {} }
    poseState.viewer = null;
    $("viewer3d").innerHTML = "";
  }
  function showPoseError(message) {
    clearViewer();
    $("pose-loading").classList.add("hidden");
    $("viewer3d").innerHTML = `<div class="pose-error">${escapeHtml(message)}</div>`;
  }

  async function renderPose(id, mode) {
    const detail = poseState.detailById[id];
    const scores = detail?.scores || [];
    const score = scores.find((item) => item.mode === mode) || scores[mode - 1] || {};
    const runId = poseState.report.run_id;
    $("pose-loading").classList.remove("hidden");
    $("pose-loading").textContent = "正在加载姿态…";
    try {
      if (!poseState.receptorText) poseState.receptorText = await fetch(artifactUrl(runId, `results/vina-runs/${runId}/receptor.pdbqt`)).then((response) => { if (!response.ok) throw new Error("受体结构读取失败"); return response.text(); });
      const poseText = await fetch(artifactUrl(runId, resultArtifactPath(runId, detail))).then((response) => { if (!response.ok) throw new Error("配体姿态读取失败"); return response.text(); });
      if (!window.$3Dmol) throw new Error("三维查看器资源未加载");
      clearViewer();
      const viewer = window.$3Dmol.createViewer("viewer3d", { backgroundColor: "#0b1118" });
      viewer.addModel(pdbqtModel(poseState.receptorText, 1), "pdb").setStyle({}, { cartoon: { color: "spectrum", opacity: 0.78 } });
      viewer.addModel(pdbqtModel(poseText, mode), "pdb").setStyle({}, { stick: { colorscheme: "cyanCarbon", radius: 0.18 }, sphere: { colorscheme: "cyanCarbon", scale: 0.25 } });
      viewer.zoomTo(); viewer.render(); poseState.viewer = viewer;
      $("pose-summary").textContent = `${detail.name} · mode ${score.mode || mode} · ${score.affinity_kcal_mol ?? "--"} kcal/mol · RMSD ${score.rmsd_lb ?? "--"}–${score.rmsd_ub ?? "--"} Å`;
      $("pose-loading").classList.add("hidden");
    } catch (error) { showPoseError(error.message || "三维姿态加载失败"); }
  }

  async function openPose(id, mode = 1) {
    poseState.activeId = id;
    const detail = poseState.detailById[id];
    const select = $("pose-mode");
    select.innerHTML = (detail?.scores || []).map((score) => `<option value="${score.mode}">Mode ${score.mode} · ${score.affinity_kcal_mol} kcal/mol</option>`).join("");
    select.value = String(mode);
    $("pose-workbench").classList.remove("hidden");
    await renderPose(id, mode);
    $("pose-workbench").scrollIntoView({ behavior: "smooth", block: "center" });
  }

  function renderResult(report) {
    const vina = report.vina || report.result?.vina;
    if (!vina) { $("result").textContent = report.error || "任务没有返回结果"; return; }
    poseState = { report, vina, viewer: null, receptorText: "", activeId: null, detailById: Object.fromEntries((vina.results || []).map((item) => [item.id, item])) };
    const rows = (vina.ranking || []).map((item) => `<tr><td>${escapeHtml(item.name)}</td><td>${item.top_affinity_kcal_mol ?? "--"}</td><td>${item.status}</td><td><button type="button" class="pose-mode-button" data-pose="${escapeHtml(item.id)}" ${item.status === "succeeded" ? "" : "disabled"}>查看姿态</button> <a href="${artifactUrl(report.run_id, `results/vina-runs/${report.run_id}/ligands/${item.id}.out.pdbqt`)}" target="_blank" rel="noreferrer">PDBQT</a></td></tr>`).join("");
    $("result").classList.remove("empty");
    $("result").innerHTML = `<div class="callout neutral">${escapeHtml(vina.provenance?.scientific_limit || "Vina 结果仅用于计算排序")}</div><table><thead><tr><th>配体</th><th>最佳 affinity (kcal/mol)</th><th>状态</th><th>姿态</th></tr></thead><tbody>${rows || '<tr><td colspan="4">没有成功的 Vina pose</td></tr>'}</tbody></table><details><summary>查看完整参数和任务 ID</summary><pre>${escapeHtml(JSON.stringify(vina.configuration, null, 2))}</pre><small>${escapeHtml(report.run_id)}</small></details>`;
    $("pose-workbench").classList.add("hidden");
    $("result").querySelectorAll("[data-pose]").forEach((button) => button.addEventListener("click", () => openPose(button.dataset.pose)));
  }

  async function poll(runId) {
    try {
      const response = await fetch(`${API}/vina/runs/${encodeURIComponent(runId)}`); const record = await response.json();
      const eventsResponse = await fetch(`${API}/runs/${encodeURIComponent(runId)}/events`); const events = await eventsResponse.json();
      const stages = events.events || []; const latest = stages[stages.length - 1]?.payload || {};
      showProgress(latest.name || latest.stage || record.current_stage || "Vina 运行中", record.status === "succeeded" || record.status === "partial" ? 1 : Math.min(.9, .12 + stages.length * .08));
      if (["succeeded", "partial", "failed"].includes(record.status)) {
        clearInterval(pollTimer); pollTimer = null; $("run-vina").disabled = false;
        if (record.status === "failed") { showError(record.error || "Vina 任务失败"); return; }
        const report = await fetch(`${API}/runs/${encodeURIComponent(runId)}/report`).then((res) => res.json());
        $("run-status").textContent = record.status === "partial" ? "任务完成，但部分配体失败，请检查日志。" : "任务完成，可查看 pose、affinity 和日志。";
        $("run-status").style.color = ""; renderResult(report);
      }
    } catch (error) { showError(error.message || "无法读取任务状态"); }
  }

  async function run() {
    if (!receptor) return showError("请先选择受体 PDBQT，或加载已验证示例。");
    if (!ligands.length) return showError("请至少添加一个配体 PDBQT。");
    $("run-vina").disabled = true; $("run-status").style.color = ""; $("run-status").textContent = "正在提交输入检查…"; showProgress("提交标准 Vina 任务", .08);
    try {
      const response = await fetch(`${API}/vina/runs`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload()) }); const data = await response.json();
      if (!response.ok) throw new Error(data.detail || "Vina 输入检查失败");
      $("run-status").textContent = `任务已提交：${data.run_id}`; pollTimer = setInterval(() => poll(data.run_id), 1400); await poll(data.run_id);
    } catch (error) { $("run-vina").disabled = false; showError(error.message || "Vina 任务提交失败"); }
  }

  async function loadDemo() {
    $("load-demo").disabled = true; $("load-demo").textContent = "正在加载示例…";
    try {
      const data = await fetch(`${API}/vina/demo`).then(async (response) => { const value = await response.json(); if (!response.ok) throw new Error(value.detail || "示例不可用"); return value; });
      receptor = { ...data.receptor, bytes: Math.round(data.receptor.content_base64.length * .75) }; ligands = data.ligands.map((item) => ({ ...item, bytes: Math.round(item.content_base64.length * .75) }));
      ["cx","cy","cz"].forEach((id, i) => $(id).value = data.box.center[i]); ["sx","sy","sz"].forEach((id, i) => $(id).value = data.box.size[i]);
      Object.entries(data.params).forEach(([key, value]) => { if ($(key)) $(key).value = value; }); markReceptor(receptor); renderLigands(); $("run-status").textContent = "已加载 1IEP 实验受体和 imatinib 控制配体，可以直接运行。";
    } catch (error) { showError(error.message || "示例加载失败"); } finally { $("load-demo").disabled = false; $("load-demo").textContent = "加载已验证示例"; }
  }

  $("pose-mode").addEventListener("change", () => { if (poseState.activeId) renderPose(poseState.activeId, Number($("pose-mode").value)); });
  $("pose-reset").addEventListener("click", () => { if (poseState.viewer) { poseState.viewer.zoomTo(); poseState.viewer.render(); } });
  $("pose-fullscreen").addEventListener("click", async () => {
    const target = $("pose-workbench");
    try {
      if (document.fullscreenElement) await document.exitFullscreen();
      else if (target.requestFullscreen) await target.requestFullscreen();
      else throw new Error("fullscreen unavailable");
    } catch (_) {
      target.classList.toggle("is-expanded");
    }
    $("pose-fullscreen").textContent = document.fullscreenElement || target.classList.contains("is-expanded") ? "退出全屏" : "全屏";
    if (poseState.viewer) setTimeout(() => { poseState.viewer.resize(); poseState.viewer.render(); }, 150);
  });
  document.addEventListener("fullscreenchange", () => {
    const target = $("pose-workbench");
    $("pose-fullscreen").textContent = document.fullscreenElement || target.classList.contains("is-expanded") ? "退出全屏" : "全屏";
    if (poseState.viewer) setTimeout(() => { poseState.viewer.resize(); poseState.viewer.render(); }, 150);
  });
  $("receptor-file").addEventListener("change", async (event) => { const file = event.target.files?.[0]; if (!file) return; try { markReceptor({ name: file.name, bytes: file.size, content_base64: await fileToBase64(file) }); } catch (error) { showError(error.message); } });
  $("add-ligand").addEventListener("click", () => { const input = document.createElement("input"); input.type = "file"; input.accept = ".pdbqt,.txt"; input.onchange = async () => { const file = input.files?.[0]; if (!file) return; try { const content_base64 = await fileToBase64(file); const base = file.name.replace(/\.pdbqt$|\.txt$/i, ""); ligands.push({ id: `${base}-${Date.now()}`, name: file.name, bytes: file.size, content_base64 }); renderLigands(); } catch (error) { showError(error.message); } }; input.click(); });
  $("load-demo").addEventListener("click", loadDemo); $("run-vina").addEventListener("click", run); renderLigands();
})();
