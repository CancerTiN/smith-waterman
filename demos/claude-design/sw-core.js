/*
  Shared Smith-Waterman algorithm + animation controller.
  Each skin HTML provides the markup (with matching IDs) and visual CSS.
*/
(function () {
  function boot() {
    const POINTER = { DIAG: "diag", UP: "up", LEFT: "left", NONE: "none" };
    const ARROW = { diag: "↖", up: "↑", left: "←", none: "STOP" };
    const SPEEDS = [
      { label: "极慢", delay: 1500 },
      { label: "慢", delay: 900 },
      { label: "中", delay: 450 },
      { label: "快", delay: 150 },
      { label: "极快", delay: 50 }
    ];
    const DEFAULTS = { seqA: "ATCG", seqB: "TCG", match: 2, mismatch: -1, gap: -1 };
    const EXAMPLE_SETS = {
      dna: { alphabet: "ACGT", minLength: 8, maxLength: 14, match: 3, mismatch: -3, gap: -2, motifLength: 4 },
      protein: { alphabet: "ACDEFGHIKLMNPQRSTVWYX", minLength: 9, maxLength: 16, match: 2, mismatch: -1, gap: -2, motifLength: 3 }
    };

    const state = {
      runId: 0, phase: "idle", isPaused: false, timer: null,
      fillIndex: 0, tracebackIndex: 0, computed: null,
      previewA: "", previewMid: "", previewB: "", lastTraceCell: null
    };

    const els = {
      seqA: document.getElementById("seqA"),
      seqB: document.getElementById("seqB"),
      matchScore: document.getElementById("matchScore"),
      mismatchScore: document.getElementById("mismatchScore"),
      gapScore: document.getElementById("gapScore"),
      inputErrors: document.getElementById("inputErrors"),
      runBtn: document.getElementById("runBtn"),
      pauseBtn: document.getElementById("pauseBtn"),
      stepBtn: document.getElementById("stepBtn"),
      resetBtn: document.getElementById("resetBtn"),
      dnaExample: document.getElementById("dnaExample"),
      proteinExample: document.getElementById("proteinExample"),
      speedRange: document.getElementById("speedRange"),
      speedLabel: document.getElementById("speedLabel"),
      skipFill: document.getElementById("skipFill"),
      phaseBadge: document.getElementById("phaseBadge"),
      matrixHost: document.getElementById("matrixHost"),
      matrixWrap: document.getElementById("matrixWrap"),
      traceOverlay: document.getElementById("traceOverlay"),
      formulaPanel: document.getElementById("formulaPanel"),
      alignmentPreview: document.getElementById("alignmentPreview"),
      traceStatus: document.getElementById("traceStatus"),
      sequenceContext: document.getElementById("sequenceContext"),
      resultHost: document.getElementById("resultHost"),
      copyStatus: document.getElementById("copyStatus"),
      liveRegion: document.getElementById("liveRegion")
    };

    function sanitizeInput(value) { return value.toUpperCase().replace(/\s+/g, ""); }

    function validateInputs() {
      const seqA = sanitizeInput(els.seqA.value);
      const seqB = sanitizeInput(els.seqB.value);
      const match = Number(els.matchScore.value);
      const mismatch = Number(els.mismatchScore.value);
      const gap = Number(els.gapScore.value);
      const errors = [];
      const invalidA = [...new Set(seqA.replace(/[A-Z]/g, "").split("").filter(Boolean))];
      const invalidB = [...new Set(seqB.replace(/[A-Z]/g, "").split("").filter(Boolean))];

      els.seqA.classList.toggle("invalid", invalidA.length > 0 || seqA.length === 0 || seqA.length > 24);
      els.seqB.classList.toggle("invalid", invalidB.length > 0 || seqB.length === 0 || seqB.length > 24);
      [els.matchScore, els.mismatchScore, els.gapScore].forEach((input) => input.classList.remove("invalid"));

      if (!seqA) errors.push("Sequence A 不能为空。");
      if (!seqB) errors.push("Sequence B 不能为空。");
      if (invalidA.length) errors.push(`Sequence A 含非法字符:${invalidA.join(", ")}`);
      if (invalidB.length) errors.push(`Sequence B 含非法字符:${invalidB.join(", ")}`);
      if (seqA.length > 24) errors.push("Sequence A 超过 24 个字符。");
      if (seqB.length > 24) errors.push("Sequence B 超过 24 个字符。");
      if (!Number.isFinite(match)) { errors.push("match 必须是有效数字。"); els.matchScore.classList.add("invalid"); }
      if (!Number.isFinite(mismatch)) { errors.push("mismatch 必须是有效数字。"); els.mismatchScore.classList.add("invalid"); }
      if (!Number.isFinite(gap)) { errors.push("gap 必须是有效数字。"); els.gapScore.classList.add("invalid"); }

      els.inputErrors.textContent = errors.join(" ");
      els.runBtn.disabled = errors.length > 0;
      return { ok: errors.length === 0, seqA, seqB, match, mismatch, gap, errors };
    }

    function computeSmithWaterman(seqA, seqB, params) {
      const rows = seqA.length + 1;
      const cols = seqB.length + 1;
      const H = Array.from({ length: rows }, () => Array(cols).fill(0));
      const pointer = Array.from({ length: rows }, () => Array(cols).fill(POINTER.NONE));
      const fillSteps = [];
      let maxScore = 0;
      const maxCells = [{ i: 0, j: 0 }];

      for (let i = 1; i < rows; i += 1) {
        for (let j = 1; j < cols; j += 1) {
          const charA = seqA[i - 1];
          const charB = seqB[j - 1];
          const substitution = charA === charB ? params.match : params.mismatch;
          const diag = H[i - 1][j - 1] + substitution;
          const up = H[i - 1][j] + params.gap;
          const left = H[i][j - 1] + params.gap;
          const zero = 0;
          const value = Math.max(zero, diag, up, left);

          let ptr = POINTER.NONE;
          if (value > 0 && value === diag) ptr = POINTER.DIAG;
          else if (value > 0 && value === up) ptr = POINTER.UP;
          else if (value > 0 && value === left) ptr = POINTER.LEFT;

          H[i][j] = value;
          pointer[i][j] = ptr;

          if (value > maxScore) { maxScore = value; maxCells.length = 0; maxCells.push({ i, j }); }
          else if (value === maxScore) maxCells.push({ i, j });

          fillSteps.push({ i, j, charA, charB, substitution, diag, up, left, zero, value, pointer: ptr });
        }
      }

      const positiveMaxCells = maxCells.filter((cell) => H[cell.i][cell.j] === maxScore && maxScore > 0);
      const maxCell = positiveMaxCells[0] || { i: 0, j: 0 };
      const tracebackSteps = [];
      const alignedA = [];
      const alignedB = [];
      const matchLine = [];

      let i = maxCell.i;
      let j = maxCell.j;
      while (maxScore > 0 && i > 0 && j > 0 && H[i][j] > 0 && pointer[i][j] !== POINTER.NONE) {
        const ptr = pointer[i][j];
        let addA = "", addB = "", addMatch = " ";
        let nextI = i, nextJ = j;

        if (ptr === POINTER.DIAG) {
          addA = seqA[i - 1]; addB = seqB[j - 1];
          addMatch = addA === addB ? "|" : ".";
          nextI = i - 1; nextJ = j - 1;
        } else if (ptr === POINTER.UP) {
          addA = seqA[i - 1]; addB = "-"; nextI = i - 1;
        } else if (ptr === POINTER.LEFT) {
          addA = "-"; addB = seqB[j - 1]; nextJ = j - 1;
        }

        tracebackSteps.push({ i, j, value: H[i][j], pointer: ptr, addSeq1: addA, addSeq2: addB, addMatch, nextI, nextJ });
        alignedA.push(addA); alignedB.push(addB); matchLine.push(addMatch);
        i = nextI; j = nextJ;
      }

      const result = {
        score: maxScore,
        start: maxCell,
        end: { i, j },
        alignedA: alignedA.reverse().join(""),
        alignedB: alignedB.reverse().join(""),
        matchLine: matchLine.reverse().join("")
      };

      const length = result.alignedA.length;
      const matches = [...result.matchLine].filter((c) => c === "|").length;
      const mismatches = [...result.matchLine].filter((c) => c === ".").length;
      const gaps = [...result.alignedA + result.alignedB].filter((c) => c === "-").length;
      result.stats = { length, matches, mismatches, gaps, identity: length ? matches / length : 0, gapRate: length ? gaps / length : 0 };

      return { seqA, seqB, params, H, pointer, fillSteps, tracebackSteps, maxScore, maxCells: positiveMaxCells, maxCell, result };
    }

    function buildMatrixTable(computed, revealAll = false) {
      const { seqA, seqB, H, pointer } = computed;
      const table = document.createElement("table");
      table.className = "matrix";
      table.setAttribute("aria-label", "Smith-Waterman dynamic programming scoring matrix");

      const thead = document.createElement("thead");
      const headerRow = document.createElement("tr");
      headerRow.appendChild(createHeaderCell(""));
      headerRow.appendChild(createHeaderCell("−", "col-label-0"));
      for (let j = 0; j < seqB.length; j += 1) headerRow.appendChild(createHeaderCell(seqB[j], `col-label-${j + 1}`));
      thead.appendChild(headerRow);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      for (let i = 0; i <= seqA.length; i += 1) {
        const row = document.createElement("tr");
        row.appendChild(createHeaderCell(i === 0 ? "−" : seqA[i - 1], `row-label-${i}`));
        for (let j = 0; j <= seqB.length; j += 1) {
          const cell = document.createElement("td");
          cell.id = cellId(i, j);
          cell.dataset.i = String(i);
          cell.dataset.j = String(j);
          const initialized = i === 0 || j === 0;
          cell.textContent = initialized || revealAll ? formatScore(H[i][j]) : "·";
          cell.className = initialized ? "initialized" : revealAll ? "filled" : "pending";
          cell.title = initialized || revealAll ? tooltipFor(i, j, H[i][j], pointer[i][j]) : `H(${i},${j}) 尚未填充`;
          if (revealAll && pointer[i][j] !== POINTER.NONE) appendArrow(cell, pointer[i][j]);
          row.appendChild(cell);
        }
        tbody.appendChild(row);
      }
      table.appendChild(tbody);
      els.matrixHost.replaceChildren(table);
      resetTraceOverlay();
    }

    function createHeaderCell(text, id = "") {
      const th = document.createElement("th");
      th.scope = "col";
      th.textContent = text;
      if (!text) th.classList.add("corner");
      if (id) th.id = id;
      return th;
    }

    function renderInitialMatrix() {
      const validation = validateInputs();
      const seqA = validation.seqA || DEFAULTS.seqA;
      const seqB = validation.seqB || DEFAULTS.seqB;
      const computed = computeSmithWaterman(seqA, seqB, {
        match: Number.isFinite(validation.match) ? validation.match : DEFAULTS.match,
        mismatch: Number.isFinite(validation.mismatch) ? validation.mismatch : DEFAULTS.mismatch,
        gap: Number.isFinite(validation.gap) ? validation.gap : DEFAULTS.gap
      });
      buildMatrixTable(computed, false);
      clearHighlights();
      resetTraceOverlay();
    }

    function getInputsOrShowError() {
      const validation = validateInputs();
      if (!validation.ok) return null;
      els.seqA.value = validation.seqA;
      els.seqB.value = validation.seqB;
      return validation;
    }

    function updateFormulaPanel(step, mode = "fill") {
      if (!step) { els.formulaPanel.innerHTML = `<p class="explain">等待下一步。</p>`; return; }
      if (mode === "traceback") {
        els.formulaPanel.innerHTML = `
          <p class="formula-heading"><strong>当前阶段:</strong>回溯中</p>
          <div class="formula-list">
            <div class="formula-line chosen">当前坐标:H(${step.i}, ${step.j}) = ${formatScore(step.value)}</div>
            <div class="formula-line">pointer:${ARROW[step.pointer]} (${step.pointer})</div>
            <div class="formula-line">追加:A=${step.addSeq1}, B=${step.addSeq2}, 标识=${step.addMatch === " " ? "空格" : step.addMatch}</div>
            <div class="formula-line">下一步:H(${step.nextI}, ${step.nextJ})</div>
          </div>`;
        return;
      }
      const chosen = step.pointer;
      els.formulaPanel.innerHTML = `
        <p class="formula-heading"><strong>当前阶段:</strong>矩阵填充中</p>
        <p class="explain">当前计算:H(${step.i}, ${step.j});比较 A[${step.i - 1}] = ${step.charA},B[${step.j - 1}] = ${step.charB}</p>
        <div class="formula-list">
          <div class="formula-line ${chosen === POINTER.DIAG ? "chosen" : ""}">diag = H(${step.i - 1},${step.j - 1}) + ${step.charA === step.charB ? "match" : "mismatch"} = ${formatScore(step.diag)} ${chosen === POINTER.DIAG ? "✓ ↖" : ""}</div>
          <div class="formula-line ${chosen === POINTER.UP ? "chosen" : ""}">up = H(${step.i - 1},${step.j}) + gap = ${formatScore(step.up)} ${chosen === POINTER.UP ? "✓ ↑" : ""}</div>
          <div class="formula-line ${chosen === POINTER.LEFT ? "chosen" : ""}">left = H(${step.i},${step.j - 1}) + gap = ${formatScore(step.left)} ${chosen === POINTER.LEFT ? "✓ ←" : ""}</div>
          <div class="formula-line ${chosen === POINTER.NONE ? "chosen" : ""}">zero = 0 ${chosen === POINTER.NONE ? "✓ STOP" : ""}</div>
          <div class="formula-line chosen">H(${step.i}, ${step.j}) = ${formatScore(step.value)};方向:${ARROW[step.pointer]}</div>
        </div>`;
    }

    function updateTracebackPanel(message) { els.traceStatus.textContent = message; }

    function renderSequenceContext(computed) {
      if (!computed) { els.sequenceContext.innerHTML = "Sequence A: 等待运行\nSequence B: 等待运行"; return ""; }
      if (computed.maxScore <= 0) {
        const html = [`Sequence A: ${escapeHtml(computed.seqA)}`, `Sequence B: ${escapeHtml(computed.seqB)}`, "未找到正分局部片段"].join("\n");
        els.sequenceContext.innerHTML = html; return html;
      }
      const { result } = computed;
      const html = [
        `Sequence A: ${highlightSequenceSpan(computed.seqA, result.end.i, result.start.i)}`,
        `Sequence B: ${highlightSequenceSpan(computed.seqB, result.end.j, result.start.j)}`,
        `局部窗口:A[${result.end.i}:${result.start.i}], B[${result.end.j}:${result.start.j}]`
      ].join("\n");
      els.sequenceContext.innerHTML = html;
      return html;
    }

    function highlightSequenceSpan(sequence, start, end) {
      return [escapeHtml(sequence.slice(0, start)), `<span class="selected">${escapeHtml(sequence.slice(start, end))}</span>`, escapeHtml(sequence.slice(end))].join("");
    }

    function plainSequenceContext(computed) {
      if (!computed || computed.maxScore <= 0) return [`Sequence A: ${computed ? computed.seqA : ""}`, `Sequence B: ${computed ? computed.seqB : ""}`].join("\n");
      const { result } = computed;
      return [
        `Sequence A: ${bracketSequenceSpan(computed.seqA, result.end.i, result.start.i)}`,
        `Sequence B: ${bracketSequenceSpan(computed.seqB, result.end.j, result.start.j)}`,
        `Window: A[${result.end.i}:${result.start.i}], B[${result.end.j}:${result.start.j}]`
      ].join("\n");
    }

    function bracketSequenceSpan(sequence, start, end) {
      return `${sequence.slice(0, start)}[${sequence.slice(start, end)}]${sequence.slice(end)}`;
    }

    function updateAlignmentPreview(step) {
      if (!step) { els.alignmentPreview.textContent = "A:\n \nB:"; return; }
      state.previewA = step.addSeq1 + state.previewA;
      state.previewMid = (step.addMatch === "." ? " " : step.addMatch) + state.previewMid;
      state.previewB = step.addSeq2 + state.previewB;
      els.alignmentPreview.textContent = `A: ${state.previewA}\n   ${state.previewMid}\nB: ${state.previewB}`;
    }

    function renderFinalResult(computed) {
      if (!computed || computed.maxScore <= 0) {
        els.resultHost.innerHTML = `<p class="explain"><strong>没有找到正分局部比对,请调整序列或评分参数。</strong></p><div class="context-box">${renderSequenceContextHtml(computed)}</div>`;
        return;
      }
      const { result } = computed;
      const identity = Math.round(result.stats.identity * 1000) / 10;
      const gapRate = Math.round(result.stats.gapRate * 1000) / 10;
      const text = resultText(computed);
      els.resultHost.innerHTML = `
        <div class="result-grid">
          <div class="metric"><b>${formatScore(result.score)}</b><span>比对得分</span></div>
          <div class="metric"><b>${result.stats.length}</b><span>比对长度</span></div>
          <div class="metric"><b>${result.stats.matches}</b><span>匹配数</span></div>
          <div class="metric"><b>${result.stats.gaps}</b><span>空位数</span></div>
        </div>
        <p class="explain">起点:H(${result.start.i}, ${result.start.j});终点:H(${result.end.i}, ${result.end.j});Identity:${result.stats.matches}/${result.stats.length} (${identity}%);Gaps:${result.stats.gaps}/${result.stats.length} (${gapRate}%)</p>
        <div class="context-box">${renderSequenceContextHtml(computed)}</div>
        <div class="alignment-box">${escapeHtml(text)}</div>
        <div class="copy-row"><button type="button" id="copyResult">复制结果</button></div>`;
      document.getElementById("copyResult").addEventListener("click", async () => {
        try { await navigator.clipboard.writeText(text); els.copyStatus.textContent = "已复制"; }
        catch (e) { els.copyStatus.textContent = "复制失败"; }
        window.setTimeout(() => { els.copyStatus.textContent = ""; }, 1400);
      });
    }

    function resultText(computed) {
      const { result } = computed;
      const identity = Math.round(result.stats.identity * 1000) / 10;
      const gapRate = Math.round(result.stats.gapRate * 1000) / 10;
      return [
        `Score: ${formatScore(result.score)}`,
        `Start: H(${result.start.i}, ${result.start.j}), End: H(${result.end.i}, ${result.end.j})`,
        plainSequenceContext(computed),
        `Length: ${result.stats.length}, Identity: ${result.stats.matches}/${result.stats.length} (${identity}%), Gaps: ${result.stats.gaps}/${result.stats.length} (${gapRate}%)`,
        `A: ${result.alignedA}`, `   ${result.matchLine}`, `B: ${result.alignedB}`
      ].join("\n");
    }

    function renderSequenceContextHtml(computed) {
      if (!computed || computed.maxScore <= 0) return escapeHtml(plainSequenceContext(computed));
      const { result } = computed;
      return [
        `Sequence A: ${highlightSequenceSpan(computed.seqA, result.end.i, result.start.i)}`,
        `Sequence B: ${highlightSequenceSpan(computed.seqB, result.end.j, result.start.j)}`,
        `局部窗口:A[${result.end.i}:${result.start.i}], B[${result.end.j}:${result.start.j}]`
      ].join("\n");
    }

    function clearHighlights() {
      document.querySelectorAll(".current, .dependency, .active-label, .trace-current, .stop-cell").forEach((n) => n.classList.remove("current", "dependency", "active-label", "trace-current", "stop-cell"));
    }
    function clearTransientFillHighlights() {
      document.querySelectorAll(".current, .dependency, .active-label").forEach((n) => n.classList.remove("current", "dependency", "active-label"));
    }

    function resetState(restoreDefaults = false) {
      cancelCurrentRun();
      if (restoreDefaults) {
        els.seqA.value = DEFAULTS.seqA; els.seqB.value = DEFAULTS.seqB;
        els.matchScore.value = DEFAULTS.match; els.mismatchScore.value = DEFAULTS.mismatch; els.gapScore.value = DEFAULTS.gap;
        els.skipFill.checked = false;
      }
      Object.assign(state, { phase: "idle", isPaused: false, fillIndex: 0, tracebackIndex: 0, computed: null, previewA: "", previewMid: "", previewB: "", lastTraceCell: null });
      els.pauseBtn.textContent = "暂停"; els.pauseBtn.disabled = true;
      els.resultHost.innerHTML = `<p class="explain">完成回溯后会显示得分、坐标、统计信息和三行可复制比对结果。</p>`;
      els.formulaPanel.innerHTML = `<p class="explain">点击 “运行演示” 或 “单步” 开始。填充时显示 diag、up、left、zero 四个候选值;回溯时显示 pointer 指向和追加字符。</p>`;
      renderSequenceContext(null);
      els.alignmentPreview.textContent = "A:\n \nB:";
      updateTracebackPanel("等待回溯");
      setPhase("idle");
      renderInitialMatrix();
      updateButtons();
    }

    function cancelCurrentRun() {
      state.runId += 1;
      if (state.timer) { clearTimeout(state.timer); state.timer = null; }
      state.isPaused = false;
    }

    function sleep(runId) {
      return new Promise((resolve, reject) => {
        const tick = () => {
          if (runId !== state.runId) { reject(new Error("cancelled")); return; }
          if (state.isPaused) { state.timer = window.setTimeout(tick, 80); return; }
          state.timer = window.setTimeout(resolve, getCurrentDelay());
        };
        tick();
      });
    }

    async function animateFill(runId) {
      state.phase = "fill"; setPhase("fill"); updateButtons();
      if (els.skipFill.checked) {
        buildMatrixTable(state.computed, true);
        state.fillIndex = state.computed.fillSteps.length;
        markMaxCells(state.computed);
        await sleep(runId); return;
      }
      while (state.fillIndex < state.computed.fillSteps.length) {
        stepFill();
        await sleep(runId);
      }
      markMaxCells(state.computed);
      updateFormulaPanel(null);
      els.formulaPanel.innerHTML = `<p class="formula-heading"><strong>矩阵填充完成</strong></p><p class="explain">全矩阵最高分为 ${formatScore(state.computed.maxScore)},回溯将从行优先遇到的第一个最高分单元格开始。</p>`;
      await sleep(runId);
    }

    async function animateTraceback(runId) {
      state.phase = "traceback"; setPhase("traceback"); updateButtons();
      if (state.computed.maxScore <= 0) {
        updateTracebackPanel("没有找到正分局部比对");
        renderFinalResult(state.computed);
        state.phase = "done"; setPhase("done"); updateButtons();
        return;
      }
      while (state.tracebackIndex < state.computed.tracebackSteps.length) {
        stepTraceback();
        await sleep(runId);
      }
      const stopCell = getCell(state.computed.result.end.i, state.computed.result.end.j);
      if (stopCell) stopCell.classList.add("stop-cell");
      updateTracebackPanel("回溯完成,最优局部比对已生成");
      els.formulaPanel.innerHTML = `<p class="formula-heading"><strong>回溯完成</strong></p><p class="explain">最优局部比对已生成。终止位置为 H(${state.computed.result.end.i}, ${state.computed.result.end.j})。</p>`;
      renderFinalResult(state.computed);
      state.phase = "done"; setPhase("done"); updateButtons();
    }

    async function runDemo() {
      const validation = getInputsOrShowError();
      if (!validation) return;
      cancelCurrentRun();
      const runId = state.runId;
      initializeRun(validation);
      try { await animateFill(runId); await animateTraceback(runId); }
      catch (error) { if (error.message !== "cancelled") throw error; }
    }

    function initializeRun(validation) {
      Object.assign(state, { phase: "idle", isPaused: false, fillIndex: 0, tracebackIndex: 0, previewA: "", previewMid: "", previewB: "", lastTraceCell: null });
      state.computed = computeSmithWaterman(validation.seqA, validation.seqB, { match: validation.match, mismatch: validation.mismatch, gap: validation.gap });
      buildMatrixTable(state.computed, false);
      els.resultHost.innerHTML = `<p class="explain">动画进行中,结果将在回溯完成后生成。</p>`;
      renderSequenceContext(state.computed);
      els.alignmentPreview.textContent = "A:\n \nB:";
      updateTracebackPanel("等待回溯");
      updateButtons();
    }

    function stepFill() {
      const step = state.computed.fillSteps[state.fillIndex];
      if (!step) { markMaxCells(state.computed); state.phase = "traceback"; setPhase("traceback"); updateButtons(); return; }
      clearTransientFillHighlights();
      highlightFillStep(step);
      const cell = getCell(step.i, step.j);
      cell.textContent = formatScore(step.value);
      cell.classList.remove("pending");
      cell.classList.add("filled");
      cell.title = tooltipFor(step.i, step.j, step.value, step.pointer);
      if (step.pointer !== POINTER.NONE) appendArrow(cell, step.pointer);
      updateFormulaPanel(step, "fill");
      state.fillIndex += 1;
      announce(`已填充 H(${step.i}, ${step.j}),得分 ${formatScore(step.value)}`);
    }

    function stepTraceback() {
      const step = state.computed.tracebackSteps[state.tracebackIndex];
      if (!step) { renderFinalResult(state.computed); state.phase = "done"; setPhase("done"); updateButtons(); return; }
      const current = getCell(step.i, step.j);
      const next = getCell(step.nextI, step.nextJ);
      if (state.lastTraceCell) state.lastTraceCell.classList.remove("trace-current");
      current.classList.add("trace-cell", "trace-current");
      if (next) drawTraceLine(current, next);
      updateFormulaPanel(step, "traceback");
      updateAlignmentPreview(step);
      updateTracebackPanel(`回溯 H(${step.i}, ${step.j}),方向 ${ARROW[step.pointer]}`);
      state.lastTraceCell = next || current;
      state.tracebackIndex += 1;
      announce(`回溯 H(${step.i}, ${step.j})`);
    }

    function stepOnce() {
      if (!state.computed || state.phase === "idle") {
        const validation = getInputsOrShowError();
        if (!validation) return;
        cancelCurrentRun();
        initializeRun(validation);
        state.phase = "fill"; setPhase("fill");
      }
      state.isPaused = true; els.pauseBtn.textContent = "继续";
      if (state.phase === "fill") {
        if (els.skipFill.checked) { buildMatrixTable(state.computed, true); state.fillIndex = state.computed.fillSteps.length; }
        else if (state.fillIndex < state.computed.fillSteps.length) { stepFill(); updateButtons(); return; }
        markMaxCells(state.computed);
        state.phase = "traceback"; setPhase("traceback"); updateButtons(); return;
      }
      if (state.phase === "traceback") {
        if (state.computed.maxScore <= 0 || state.tracebackIndex >= state.computed.tracebackSteps.length) {
          renderFinalResult(state.computed); state.phase = "done"; setPhase("done");
        } else stepTraceback();
        updateButtons();
      }
    }

    function highlightFillStep(step) {
      getCell(step.i, step.j).classList.add("current");
      getCell(step.i - 1, step.j - 1).classList.add("dependency");
      getCell(step.i - 1, step.j).classList.add("dependency");
      getCell(step.i, step.j - 1).classList.add("dependency");
      document.getElementById(`row-label-${step.i}`).classList.add("active-label");
      document.getElementById(`col-label-${step.j}`).classList.add("active-label");
    }

    function markMaxCells(computed) {
      computed.maxCells.forEach(({ i, j }) => {
        const cell = getCell(i, j); if (!cell) return;
        cell.classList.add("max-cell", "flash");
        if (!cell.querySelector(".cell-star")) {
          const star = document.createElement("span");
          star.className = "cell-star"; star.textContent = "★";
          cell.appendChild(star);
        }
      });
      if (computed.maxScore <= 0) updateTracebackPanel("没有找到正分局部比对");
      else updateTracebackPanel(`最高分 ${formatScore(computed.maxScore)},共有 ${computed.maxCells.length} 个最高分单元格`);
    }

    function appendArrow(cell, pointer) {
      cell.querySelector(".cell-arrow")?.remove();
      const arrow = document.createElement("span");
      arrow.className = "cell-arrow"; arrow.textContent = ARROW[pointer];
      cell.appendChild(arrow);
    }

    function drawTraceLine(fromCell, toCell) {
      const wrapBox = els.matrixWrap.getBoundingClientRect();
      const a = fromCell.getBoundingClientRect();
      const b = toCell.getBoundingClientRect();
      const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
      line.setAttribute("x1", String(a.left + a.width / 2 - wrapBox.left));
      line.setAttribute("y1", String(a.top + a.height / 2 - wrapBox.top));
      line.setAttribute("x2", String(b.left + b.width / 2 - wrapBox.left));
      line.setAttribute("y2", String(b.top + b.height / 2 - wrapBox.top));
      els.traceOverlay.appendChild(line);
    }

    function resetTraceOverlay() {
      els.traceOverlay.replaceChildren();
      requestAnimationFrame(() => {
        const box = els.matrixHost.getBoundingClientRect();
        els.traceOverlay.setAttribute("viewBox", `0 0 ${Math.max(1, box.width)} ${Math.max(1, box.height)}`);
        els.traceOverlay.style.width = `${Math.max(1, box.width)}px`;
        els.traceOverlay.style.height = `${Math.max(1, box.height)}px`;
      });
    }

    function getCell(i, j) { return document.getElementById(cellId(i, j)); }
    function cellId(i, j) { return `cell-${i}-${j}`; }
    function tooltipFor(i, j, value, pointer) { return `H(${i},${j})=${formatScore(value)}, ptr=${ARROW[pointer]}`; }
    function formatScore(value) { return Number.isInteger(value) ? String(value) : String(Math.round(value * 1000) / 1000); }
    function getCurrentDelay() { return SPEEDS[Number(els.speedRange.value)].delay; }

    function updateSpeedLabel() {
      const speed = SPEEDS[Number(els.speedRange.value)];
      els.speedLabel.textContent = `${speed.label} · ${speed.delay}ms`;
    }

    function setPhase(phase) {
      const labels = { idle: "待运行", fill: "矩阵填充中", traceback: "回溯中", done: "完成" };
      els.phaseBadge.textContent = `${labels[phase] || phase}`;
      els.phaseBadge.dataset.phase = phase;
      announce(`状态:${labels[phase] || phase}`);
    }

    function updateButtons() {
      const isAnimating = state.phase === "fill" || state.phase === "traceback";
      els.pauseBtn.disabled = !isAnimating || state.phase === "done";
      els.pauseBtn.textContent = state.isPaused ? "继续" : "暂停";
      els.stepBtn.disabled = state.phase === "done";
    }

    function escapeHtml(text) {
      return text.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;").replace(/'/g, "&#039;");
    }

    function announce(message) { if (els.liveRegion) els.liveRegion.textContent = message; }

    function setExample(kind) {
      const example = generateRandomExample(kind);
      els.seqA.value = example.seqA; els.seqB.value = example.seqB;
      els.matchScore.value = example.match; els.mismatchScore.value = example.mismatch; els.gapScore.value = example.gap;
      cancelCurrentRun();
      state.computed = null; state.phase = "idle"; setPhase("idle");
      validateInputs();
      renderInitialMatrix();
      updateButtons();
    }

    function generateRandomExample(kind) {
      const config = EXAMPLE_SETS[kind];
      const lengthA = randomInt(config.minLength, config.maxLength);
      const lengthB = randomInt(config.minLength - 1, config.maxLength - 1);
      const motifLength = Math.min(config.motifLength + randomInt(0, 2), lengthA - 2, lengthB - 2);
      const motif = randomSequence(config.alphabet, motifLength);
      const seqA = embedMotif(randomSequence(config.alphabet, lengthA), motif);
      let seqB = embedMotif(randomSequence(config.alphabet, lengthB), maybeMutateMotif(motif, config.alphabet));
      if (seqA === seqB) seqB = mutateAt(seqB, randomInt(0, seqB.length - 1), config.alphabet);
      return { seqA, seqB, match: config.match, mismatch: config.mismatch, gap: config.gap };
    }

    function randomSequence(alphabet, length) { return Array.from({ length }, () => alphabet[randomInt(0, alphabet.length - 1)]).join(""); }
    function embedMotif(sequence, motif) { const start = randomInt(0, sequence.length - motif.length); return `${sequence.slice(0, start)}${motif}${sequence.slice(start + motif.length)}`; }
    function maybeMutateMotif(motif, alphabet) {
      const chars = motif.split("");
      if (chars.length > 2 && Math.random() < 0.55) chars[randomInt(0, chars.length - 1)] = "-";
      const mutated = chars.join("").replace(/-/g, "");
      return mutated || motif;
    }
    function mutateAt(sequence, index, alphabet) {
      const current = sequence[index];
      const choices = alphabet.split("").filter((c) => c !== current);
      const replacement = choices[randomInt(0, choices.length - 1)];
      return `${sequence.slice(0, index)}${replacement}${sequence.slice(index + 1)}`;
    }
    function randomInt(min, max) { return Math.floor(Math.random() * (max - min + 1)) + min; }

    els.runBtn.addEventListener("click", runDemo);
    els.pauseBtn.addEventListener("click", () => { state.isPaused = !state.isPaused; updateButtons(); updateTracebackPanel(state.isPaused ? "已暂停" : "继续播放"); });
    els.stepBtn.addEventListener("click", () => { if (state.phase === "fill" || state.phase === "traceback") state.isPaused = true; stepOnce(); updateButtons(); });
    els.resetBtn.addEventListener("click", () => resetState(true));
    els.dnaExample.addEventListener("click", () => setExample("dna"));
    els.proteinExample.addEventListener("click", () => setExample("protein"));
    els.speedRange.addEventListener("input", updateSpeedLabel);

    [els.seqA, els.seqB, els.matchScore, els.mismatchScore, els.gapScore].forEach((input) => {
      input.addEventListener("input", () => { validateInputs(); if (state.phase === "idle") renderInitialMatrix(); });
      input.addEventListener("blur", () => { if (input === els.seqA || input === els.seqB) input.value = sanitizeInput(input.value); });
    });

    updateSpeedLabel();
    validateInputs();
    renderInitialMatrix();
    updateButtons();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
