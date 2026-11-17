// Carquet site: tabs, copy buttons, syntax colouring, and benchmark charts built from
// bench.json (a run_benchmark.py report, copied as-is).
(function () {
  "use strict";

  const LIBS = [
    { id: "carquet", name: "carquet", color: "var(--s-carquet)" },
    { id: "arrow_cpp", name: "Arrow C++", color: "var(--s-arrow)" },
    { id: "pyarrow", name: "PyArrow", color: "var(--s-pyarrow)" },
    { id: "fastparquet", name: "fastparquet", color: "var(--s-fastparquet)" },
    { id: "polars", name: "polars", color: "var(--s-polars)" },
  ];
  const CODECS = ["none", "snappy", "zstd", "lz4"];
  const $ = (s, el = document) => el.querySelector(s);
  const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

  // ── syntax colouring (tiny, display only) ──────────────────────────────
  const RULES = {
    c: [
      ["c", /\/\*[\s\S]*?\*\/|\/\/[^\n]*/y],
      ["s", /"(?:[^"\\\n]|\\.)*"|'(?:[^'\\\n]|\\.)*'/y],
      ["pp", /#\s*\w+/y],
      ["kw", /\b(?:int|const|void|return|if|else|for|while|do|struct|sizeof|char|double|bool|true|false|NULL|int64_t|int32_t|int16_t|uint8_t|long|unsigned)\b/y],
      ["ty", /\bcarquet_\w+_t\b|\bArrow(?:Schema|Array)\b/y],
      ["mc", /\b(?:CARQUET_[A-Z0-9_]+|TRY)\b/y],
      ["fn", /\b[a-z_]\w*(?=\()/y],
      ["n", /\b\d[\d_.]*u?\b/y],
    ],
    py: [
      ["c", /#[^\n]*/y],
      ["s", /"(?:[^"\\\n]|\\.)*"|'(?:[^'\\\n]|\\.)*'/y],
      ["kw", /\b(?:import|as|with|for|in|from|def|return|print)\b/y],
      ["fn", /\b[a-z_]\w*(?=\()/y],
      ["n", /\b\d[\d_.]*\b/y],
    ],
    sh: [
      ["tp", /^(?:\$|PS>)(?= )/my],
      ["c", /(?<=\s)#[^\n]*/y],
      ["s", /"(?:[^"\\\n]|\\.)*"/y],
    ],
    cmake: [
      ["c", /#[^\n]*/y],
      ["s", /"(?:[^"\\\n]|\\.)*"/y],
      ["fn", /^\w+(?=\()/my],
      ["mc", /\b[A-Z][A-Z0-9_]{2,}\b/y],
    ],
  };
  function highlight(code, lang) {
    const rules = RULES[lang];
    const src = code.textContent;
    let out = "", plain = "", i = 0;
    outer: while (i < src.length) {
      for (const [cls, re] of rules) {
        re.lastIndex = i;
        const m = re.exec(src);
        if (m && m[0].length) {
          out += esc(plain) + `<span class="${cls}">${esc(m[0])}</span>`;
          plain = ""; i += m[0].length;
          continue outer;
        }
      }
      // advance over a whole identifier so rules only match at word starts
      const w = /\w+|[^\w]/y; w.lastIndex = i;
      const t = w.exec(src)[0];
      plain += t; i += t.length;
    }
    code.innerHTML = out + esc(plain);
  }
  document.querySelectorAll("code[class^=lang-]").forEach((c) => highlight(c, c.className.slice(5)));

  // ── SIMD lanes: same work, fewer steps (scalar 16 ticks, 8-wide 2 ticks, 16-wide 1 tick)
  const tick = 3.2 / 16;
  [[".lane.scalar", 1], [".lane.vec", 8], [".lane.wide", 16]].forEach(([sel, w]) => {
    document.querySelectorAll(sel + " i").forEach((c, i) => (c.style.animationDelay = (Math.floor(i / w) * tick).toFixed(2) + "s"));
  });

  // ── copy buttons ───────────────────────────────────────────────────────
  function copy(text, done) {
    if (navigator.clipboard) navigator.clipboard.writeText(text).then(done, () => {});
  }
  document.querySelectorAll(".panel pre").forEach((pre) => {
    const b = document.createElement("button");
    b.type = "button"; b.className = "copy-btn"; b.textContent = "Copy";
    b.addEventListener("click", () => {
      const text = pre.querySelector("code").textContent.replace(/^(\$|PS>) /gm, "");
      copy(text, () => { b.textContent = "Copied"; setTimeout(() => (b.textContent = "Copy"), 1400); });
    });
    pre.appendChild(b);
  });

  // ── tabs (each .tabs list drives the panels of its own container) ──────
  document.querySelectorAll(".tabs").forEach((list) => {
    const box = list.parentElement;
    const tabs = list.querySelectorAll("[role=tab]");
    tabs.forEach((t) => t.addEventListener("click", () => {
      tabs.forEach((x) => x.setAttribute("aria-selected", String(x === t)));
      box.querySelectorAll(":scope > .panel").forEach((p) => (p.hidden = p.dataset.panel !== t.dataset.tab));
    }));
  });

  // ── benchmark ──────────────────────────────────────────────────────────
  const fmtMs = (v) => (v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v.toFixed(2)) + " ms";
  const fmtV = (v) => (v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v.toFixed(2));
  const fmtX = (v) => (v >= 10 ? v.toFixed(0) : v.toFixed(1)) + "×";
  const fmtRows = (n) => (n >= 1e6 ? n / 1e6 + "M" : n / 1e3 + "K");
  const fmtMB = (b) => (b / 1e6).toFixed(1) + " MB";

  fetch("bench.json").then((r) => r.json()).then(init).catch(() => {
    $("#chart").textContent = "Benchmark data unavailable.";
  });

  function niceTicks(max) {
    const raw = max / 5, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 2.5, 5, 10].map((k) => k * mag).find((s) => s >= raw);
    const top = Math.ceil(max / step) * step;
    const ticks = [];
    for (let v = 0; v <= top + step / 2; v += step) ticks.push(+v.toPrecision(6));
    return { top, ticks };
  }

  function init(data) {
    const res = data.results;
    const libs = LIBS.filter((l) => res.some((r) => r.library === l.id));
    const sizes = [...new Set(res.map((r) => r.rows))].sort((a, b) => a - b);
    const get = (lib, rows, codec) => res.find((r) => r.library === lib && r.rows === rows && r.compression === codec);
    const ver = (id) => (id === "carquet" ? data.carquet_version : data[id + "_version"]);
    const state = { op: "read_ms", rows: sizes[sizes.length - 1] };

    $("#legend").innerHTML = libs.map((l) =>
      `<span><i style="background:${l.color}"></i>${l.name} <span style="color:var(--muted)">${ver(l.id) || ""}</span></span>`).join("");

    const seg = $("#size-seg");
    seg.innerHTML = sizes.map((n) => `<button type="button" role="tab" data-rows="${n}">${fmtRows(n)} rows</button>`).join("");
    const sync = () => {
      seg.querySelectorAll("button").forEach((b) => b.setAttribute("aria-selected", String(+b.dataset.rows === state.rows)));
      $("#op-seg").querySelectorAll("button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.op === state.op)));
    };
    seg.addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) { state.rows = +b.dataset.rows; sync(); render(); } });
    $("#op-seg").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) { state.op = b.dataset.op; sync(); render(); } });
    sync();

    // headline stats (largest size, against Arrow C++ when present)
    const big = sizes[sizes.length - 1];
    const rival = libs.find((l) => l.id === "arrow_cpp") || libs.find((l) => l.id !== "carquet");
    const ratio = (op, c) => { const a = get(rival.id, big, c), b = get("carquet", big, c); return a && b ? a[op] / b[op] : NaN; };
    const range = (op) => { const v = ["snappy", "zstd", "lz4"].map((c) => ratio(op, c)).filter(isFinite); return [Math.min(...v), Math.max(...v)]; };
    const [rw0, rw1] = range("write_ms"), [rr0, rr1] = range("read_ms");
    const none = get("carquet", big, "none");
    $("#stats").innerHTML = [
      [`<em>${fmtX(ratio("read_ms", "none"))}</em>`, `faster uncompressed reads than ${rival.name}`],
      [`${fmtX(rr0).slice(0, -1)}–${fmtX(rr1)}`, `faster compressed reads than ${rival.name}`],
      [`${fmtX(rw0).slice(0, -1)}–${fmtX(rw1)}`, `faster compressed writes than ${rival.name}`],
      [fmtMs(none.read_ms), `to read ${fmtRows(big)} rows (mmap zero-copy)`],
    ].map(([b, s]) => `<div class="stat"><b>${b}</b><span>${s}</span></div>`).join("");

    const sys = data.system;
    $("#bench-sub").textContent =
      `Same data, same machine: ${sys.cpu}, ${sys.ram_gb} GB RAM. INT64 + DOUBLE + INT32 columns, lower is better. Pick an operation and a size.`;
    const it = data.config.bench_iterations;
    $("#method").textContent =
      `Trimmed median of ${it.large}–${it.small} runs after ${data.config.warmup_iterations} warm-ups, page cache purged between write and read, ` +
      `measured ${data.timestamp}. Uncompressed reads use carquet's mmap zero-copy path. ` +
      `Reproduce with python3 benchmark/run_benchmark.py --skip-xlarge --fastparquet --polars.`;

    const tip = $("#tip");
    $("#chart").addEventListener("mousemove", (e) => {
      const bar = e.target.closest("[data-tip]");
      if (!bar) { tip.classList.remove("on"); return; }
      tip.innerHTML = bar.dataset.tip;
      tip.classList.add("on");
      const x = Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 8);
      tip.style.left = x + "px"; tip.style.top = e.clientY + 14 + "px";
    });
    $("#chart").addEventListener("mouseleave", () => tip.classList.remove("on"));

    function badge(op, rows, c) {
      const me = get("carquet", rows, c);
      const others = libs.filter((l) => l.id !== "carquet").map((l) => get(l.id, rows, c)).filter(Boolean);
      if (!me || !others.length) return "";
      const best = others.reduce((a, b) => (a[op] < b[op] ? a : b));
      const x = best[op] / me[op];
      const who = libs.find((l) => l.id === best.library).name;
      if (x >= 1.05) return `<span class="badge">${fmtX(x)}<i> faster</i></span><small>vs ${who}</small>`;
      if (x >= 0.95) return `<span class="badge flat">on par</span><small>with ${who}</small>`;
      return `<span class="badge flat">${fmtX(1 / x)}<i> slower</i></span><small>vs ${who}</small>`;
    }

    function render() {
      const { op, rows } = state;
      const all = CODECS.flatMap((c) => libs.map((l) => get(l.id, rows, c)).filter(Boolean));
      const { top, ticks } = niceTicks(Math.max(...all.map((r) => r[op])));
      const opName = op === "read_ms" ? "read" : "write";

      const grid = ticks.map((t) =>
        `<div class="gl" style="bottom:${(t / top) * 100}%"><span>${+t.toPrecision(4)}</span></div>`).join("");
      const groups = CODECS.map((c) => {
        const me = get("carquet", rows, c);
        const bars = libs.map((l) => {
          const r = get(l.id, rows, c);
          if (!r) return "";
          const vsMe = me && l.id !== "carquet" ? `<br>carquet is ${fmtX(r[op] / me[op])} faster` : "";
          const t = `<b>${l.name}</b> · ${c} · ${fmtRows(rows)} rows<br>${opName} ${fmtMs(r[op])}` +
            `<span class="m">${vsMe}<br>file ${fmtMB(r.file_bytes)}</span>`;
          return `<div class="col${l.id === "carquet" ? " me" : ""}" data-tip="${t.replace(/"/g, "&quot;")}">
            <div class="vbar" style="height:${Math.max(0.4, (r[op] / top) * 100)}%;background:${l.color}"><span>${fmtV(r[op])}</span></div>
          </div>`;
        }).join("");
        return `<div class="grp">${bars}</div>`;
      }).join("");
      const xl = CODECS.map((c) =>
        `<div class="xl"><b>${c === "none" ? "uncompressed" : c}</b>${badge(op, rows, c)}</div>`).join("");

      $("#chart").innerHTML =
        `<div class="axis-t">${opName} time, ms · ${fmtRows(rows)} rows · lower is better</div>
         <div class="plot">${grid}<div class="grps">${groups}</div></div>
         <div class="xls">${xl}</div>`;

      const head = `<tr><th>codec</th>${libs.map((l) => `<th>${l.name}</th>`).join("")}</tr>`;
      const body = CODECS.map((c) => `<tr><td>${c}</td>${libs.map((l) => {
        const r = get(l.id, rows, c); return `<td>${r ? fmtMs(r[op]) : "–"}</td>`;
      }).join("")}</tr>`).join("");
      $("#table").innerHTML = `<table><caption style="text-align:left;padding:4px 0;color:var(--muted)">${opName} time, ${fmtRows(rows)} rows</caption>${head}${body}</table>`;
    }
    render();
  }
})();
