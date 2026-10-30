/**
 * Photon pictures of the two BEAM process models — the physical systems the models stand for.
 *
 *   entangle(svg, data, opts): a Bell test with polarization-entangled photon pairs (Aspect 1982, Weihs 1998).
 *     A crystal sends one photon to Alice and one to Bob; each passes a polarizer set at an angle and a
 *     polarizing beam splitter with two detectors: a click on "+" is the outcome +1, on "−" it is −1.
 *     A polarizer turns half as much as the qubit's measurement axis, so the model's settings 0°/90° and
 *     45°/135° are the polarizers' 0°/45° and 22.5°/67.5°. The pairs played are the run's first rounds.
 *   teleport(svg, data, opts): teleportation of a photon's polarization (Bouwmeester et al. 1997; over
 *     1,400 km to the Micius satellite, Ren et al. 2017). Photon 1 carries ψ; a crystal makes the entangled
 *     photons 2 and 3; Alice's Bell-state analyser joins 1 and 2 and gives two bits; they go by cable to Bob,
 *     whose switchable wave plate (a Pockels cell) applies the correction to photon 3, waiting in a fibre
 *     delay. ψ is drawn as its polarization ellipse: |0⟩ = H, |1⟩ = V, ψ = cos(θ/2)|H⟩ + e^{iφ} sin(θ/2)|V⟩.
 *
 * data is the BEAM service's answer (/quantum/entangle with `pairs` and `settings`, /quantum/migrate with
 * `samples`). Both return a stop() function. No framework: the site loads it as a module and the
 * TechUPGRADE panel imports the same file.
 */
const NS = 'http://www.w3.org/2000/svg';
const DEF = { fg: '#c8d4e8', muted: '#7a879c', line: 'rgba(120,170,255,.35)', a: '#6ab0ff', b: '#a8ff3e',
  hot: '#f5c842', bg: 'rgba(120,170,255,.06)', font: 'IBM Plex Mono, monospace' };

function el(parent, tag, attrs = {}, text) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (text != null) e.textContent = text;
  if (parent) parent.appendChild(e);
  return e;
}
const lerp = (a, b, t) => a + (b - a) * t;
const ease = t => t < 0 ? 0 : t > 1 ? 1 : t * t * (3 - 2 * t);
const CLICK = .78;                                                 // when a pair's detectors fire
const seg = (p, a, b) => ease((p - a) / (b - a));                 // 0 → 1 while p goes from a to b

function loop(duration, count, frame) {
  let raf = 0, t0 = performance.now(), stopped = false;
  const tick = now => {
    if (stopped) return;
    const t = Math.max(0, now - t0) / duration, i = Math.floor(t) % count;   // a frame's time can precede t0
    frame(i, t - Math.floor(t));
    raf = requestAnimationFrame(tick);
  };
  raf = requestAnimationFrame(tick);
  return () => { stopped = true; cancelAnimationFrame(raf); };
}

// ── polarization: a Jones vector and its ellipse ─────────────────────────────────────────────────────────
// a qubit (θ, φ in degrees) as complex amplitudes [[re, im] of H, [re, im] of V]
function jones(thetaDeg, phiDeg) {
  const t = thetaDeg * Math.PI / 360, f = phiDeg * Math.PI / 180;
  return [[Math.cos(t), 0], [Math.sin(t) * Math.cos(f), Math.sin(t) * Math.sin(f)]];
}
const gateX = ([h, v]) => [v, h];                                  // H ↔ V: a half-wave plate at 45°
const gateZ = ([h, v]) => [h, [-v[0], -v[1]]];                     // V picks up π: a half-wave plate at 0°
function ellipsePath(j, cx, cy, r) {
  const pts = [];
  for (let k = 0; k <= 48; k++) {
    const w = k / 48 * 2 * Math.PI, c = Math.cos(w), s = Math.sin(w);
    const x = j[0][0] * c + j[0][1] * s, y = j[1][0] * c + j[1][1] * s;   // Re(E e^{-iωt})
    pts.push(`${(cx + x * r).toFixed(1)},${(cy - y * r).toFixed(1)}`);
  }
  return 'M' + pts.join('L') + 'Z';
}
function polGlyph(g, j, cx, cy, r, color, label, o) {
  el(g, 'circle', { cx, cy, r: r + 6, fill: o.bg, stroke: o.line });
  el(g, 'line', { x1: cx - r - 4, y1: cy, x2: cx + r + 4, y2: cy, stroke: o.line, 'stroke-dasharray': '2 3' });
  el(g, 'line', { x1: cx, y1: cy - r - 4, x2: cx, y2: cy + r + 4, stroke: o.line, 'stroke-dasharray': '2 3' });
  el(g, 'path', { d: ellipsePath(j, cx, cy, r), fill: 'none', stroke: color, 'stroke-width': 2 });
  if (label) el(g, 'text', { x: cx, y: cy + r + 20, fill: o.muted, 'font-size': 13, 'text-anchor': 'middle', 'font-family': o.font }, label);
}

// ── entanglement: a Bell test with photon pairs ──────────────────────────────────────────────────────────
export function entangle(svg, data, opts = {}) {
  const o = { ...DEF, ...opts };
  svg.setAttribute('viewBox', '0 0 760 256');
  svg.replaceChildren();
  const pol = { alice: data.settings.alice.map(x => x / 2), bob: data.settings.bob.map(x => x / 2) };
  const pairs = data.pairs ?? [];
  if (!pairs.length) return () => {};
  const txt = (x, y, s, attrs = {}) => el(svg, 'text', { x, y, fill: o.muted, 'font-size': 14, 'text-anchor': 'middle', 'font-family': o.font, ...attrs }, s);

  // source: a pump laser into a crystal
  el(svg, 'line', { x1: 380, y1: 168, x2: 380, y2: 112, stroke: '#b07cff', 'stroke-width': 3, opacity: .7 });
  txt(380, 182, 'pump laser');
  el(svg, 'rect', { x: 366, y: 88, width: 28, height: 28, fill: 'rgba(176,124,255,.18)', stroke: '#b07cff', transform: 'rotate(45 380 102)' });
  txt(380, 70, 'crystal: one pair at a time');
  el(svg, 'line', { x1: 150, y1: 102, x2: 610, y2: 102, stroke: o.line, 'stroke-dasharray': '3 5' });

  // a side: polarizer (a disk with its axis), a polarizing beam splitter, detectors + and −
  const side = (name, dir, color) => {
    const xp = 380 + dir * 170, xb = 380 + dir * 250, xd = xb + dir * 62;
    const g = el(svg, 'g');
    el(g, 'circle', { cx: xp, cy: 102, r: 26, fill: o.bg, stroke: color });
    const axis = el(g, 'line', { x1: xp - 22, y1: 102, x2: xp + 22, y2: 102, stroke: color, 'stroke-width': 2.5 });
    const ang = el(g, 'text', { x: xp, y: 150, fill: color, 'font-size': 16, 'text-anchor': 'middle', 'font-family': o.font }, '');
    el(g, 'text', { x: xp, y: 164, fill: o.muted, 'font-size': 13, 'text-anchor': 'middle', 'font-family': o.font }, 'polarizer');
    el(g, 'rect', { x: xb - 14, y: 88, width: 28, height: 28, fill: o.bg, stroke: o.fg });
    el(g, 'line', { x1: xb - 14, y1: 116, x2: xb + 14, y2: 88, stroke: o.fg });
    el(g, 'line', { x1: xb, y1: 88, x2: xb, y2: 46, stroke: o.line });
    const dPlus = el(g, 'rect', { x: xd - 12, y: 90, width: 24, height: 24, rx: 12, fill: 'none', stroke: o.fg });
    const dMinus = el(g, 'rect', { x: xb - 12, y: 22, width: 24, height: 24, rx: 12, fill: 'none', stroke: o.fg });
    el(g, 'line', { x1: xb + dir * 14, y1: 102, x2: xd - dir * 12, y2: 102, stroke: o.line });
    el(g, 'text', { x: xd, y: 132, fill: o.fg, 'font-size': 17, 'text-anchor': 'middle', 'font-family': o.font }, '+1');
    el(g, 'text', { x: xb + dir * 28, y: 38, fill: o.fg, 'font-size': 17, 'text-anchor': 'middle', 'font-family': o.font }, '−1');
    el(g, 'text', { x: xb, y: 196, fill: color, 'font-size': 17, 'text-anchor': 'middle', 'font-family': o.font, 'letter-spacing': 2 }, name);
    const photon = el(svg, 'circle', { r: 6, fill: color, opacity: 0 });
    const glyph = el(svg, 'line', { stroke: '#fff', 'stroke-width': 2, opacity: 0 });
    return { xp, xb, xd, axis, ang, dPlus, dMinus, photon, glyph, color, dir };
  };
  const A = side('ALICE', -1, o.a), B = side('BOB', 1, o.b);
  const tally2 = el(svg, 'text', { x: 380, y: 246, fill: o.muted, 'font-size': 12, 'text-anchor': 'middle', 'font-family': o.font }, '');
  const tally = el(svg, 'text', { x: 380, y: 226, fill: o.fg, 'font-size': 14, 'text-anchor': 'middle', 'font-family': o.font }, '');
  const counts = {};

  const place = (S, angle, out, p) => {
    // rotate the polarizer to this round's angle; the photon flies, is analysed, one detector clicks
    const a = angle * Math.PI / 180;
    S.axis.setAttribute('x1', S.xp - 22 * Math.cos(a)); S.axis.setAttribute('y1', 102 + 22 * Math.sin(a));
    S.axis.setAttribute('x2', S.xp + 22 * Math.cos(a)); S.axis.setAttribute('y2', 102 - 22 * Math.sin(a));
    S.ang.textContent = `${angle}°`;
    const f1 = seg(p, .08, .45), f2 = seg(p, .45, .62), f3 = seg(p, .62, .78);
    let x = lerp(380, S.xp, f1), y = 102;
    if (f2 > 0) x = lerp(S.xp, S.xb, f2);
    if (f3 > 0) { if (out > 0) x = lerp(S.xb, S.xd, f3); else y = lerp(102, 34, f3); }
    S.photon.setAttribute('cx', x); S.photon.setAttribute('cy', y);
    S.photon.setAttribute('opacity', p < .06 || p > .86 ? 0 : 1);
    // before the polarizer the photon has no polarization of its own (no glyph); after it, along or across the axis
    const pa = out > 0 ? a : a + Math.PI / 2, show = p > .45 && p < .86;
    S.glyph.setAttribute('opacity', show ? 1 : 0);
    S.glyph.setAttribute('x1', x - 9 * Math.cos(pa)); S.glyph.setAttribute('y1', y + 9 * Math.sin(pa));
    S.glyph.setAttribute('x2', x + 9 * Math.cos(pa)); S.glyph.setAttribute('y2', y - 9 * Math.sin(pa));
    const lit = p > CLICK;
    S.dPlus.setAttribute('fill', lit && out > 0 ? o.hot : 'none');
    S.dMinus.setAttribute('fill', lit && out < 0 ? o.hot : 'none');
  };

  let last = -1;
  return loop(1900, pairs.length, (i, p) => {
    const r = pairs[i], aa = pol.alice[r.a], bb = pol.bob[r.b];
    place(A, aa, r.alice, p); place(B, bb, r.bob, p);
    if (p > CLICK && last !== i) {                            // the caption changes with the clicks
      last = i;
      if (i === 0) for (const k in counts) delete counts[k];
      const k = `${aa}/${bb}`, c = counts[k] ?? (counts[k] = { same: 0, n: 0 });
      c.n++; if (r.alice === r.bob) c.same++;
      tally.textContent = `pair ${i + 1} of ${pairs.length}: Alice ${r.alice > 0 ? '+1' : '−1'}, Bob ${r.bob > 0 ? '+1' : '−1'} — ${r.alice === r.bob ? 'the same' : 'opposite'}`;
      tally2.textContent = 'the same so far: ' + Object.entries(counts).map(([s, v]) => `${s}° ${v.same}/${v.n}`).join(' · ');
    }
  });
}

// ── teleportation: a photon's polarization ───────────────────────────────────────────────────────────────
export function teleport(svg, data, opts = {}) {
  const o = { ...DEF, ...opts };
  svg.setAttribute('viewBox', '0 0 760 290');
  svg.replaceChildren();
  const runs = data.samples ?? [];
  if (!runs.length) return () => {};
  const txt = (x, y, s, attrs = {}) => el(svg, 'text', { x, y, fill: o.muted, 'font-size': 13, 'text-anchor': 'middle', 'font-family': o.font, ...attrs }, s);

  // node A and node B as two labs
  el(svg, 'rect', { x: 6, y: 18, width: 300, height: 200, rx: 8, fill: 'none', stroke: o.a, 'stroke-dasharray': '4 4', opacity: .6 });
  el(svg, 'rect', { x: 470, y: 18, width: 284, height: 200, rx: 8, fill: 'none', stroke: o.b, 'stroke-dasharray': '4 4', opacity: .6 });
  txt(156, 12, 'ALICE · node A', { fill: o.a, 'font-size': 14, 'letter-spacing': 2 });
  txt(612, 12, 'BOB · node B', { fill: o.b, 'font-size': 14, 'letter-spacing': 2 });

  // photon 1 in, the pair source, Alice's Bell-state analyser
  const P1 = [92, 70], BSA = [200, 70], SRC = [200, 175], PC = [600, 175], OUT = [720, 175];
  el(svg, 'line', { x1: P1[0], y1: P1[1], x2: BSA[0], y2: BSA[1], stroke: o.line });
  el(svg, 'line', { x1: SRC[0], y1: SRC[1], x2: BSA[0], y2: BSA[1] + 22, stroke: o.line });
  el(svg, 'path', { d: `M${SRC[0]} ${SRC[1]} L430 ${SRC[1]} m0 0 a14 14 0 1 1 0.1 0 L${PC[0] - 22} ${PC[1]}`, fill: 'none', stroke: o.line });
  txt(430, 212, 'fibre delay: photon 3', { 'font-size': 12 }); txt(430, 223, 'waits for the bits', { 'font-size': 12 });
  el(svg, 'line', { x1: PC[0] + 22, y1: PC[1], x2: OUT[0] - 20, y2: OUT[1], stroke: o.line });
  el(svg, 'rect', { x: SRC[0] - 12, y: SRC[1] - 12, width: 24, height: 24, fill: 'rgba(176,124,255,.18)', stroke: '#b07cff', transform: `rotate(45 ${SRC[0]} ${SRC[1]})` });
  txt(SRC[0], SRC[1] + 32, 'crystal: photons 2 + 3 entangled');
  el(svg, 'rect', { x: BSA[0] - 36, y: BSA[1] - 22, width: 72, height: 52, rx: 4, fill: o.bg, stroke: o.fg });
  txt(BSA[0], BSA[1] + 2, 'Bell-state', { fill: o.fg }); txt(BSA[0], BSA[1] + 14, 'analyser', { fill: o.fg });
  const dets = ['00', '01', '10', '11'].map((b, k) => {
    const x = BSA[0] + 40, y = BSA[1] - 30 + k * 18;
    const d = el(svg, 'circle', { cx: x + 22, cy: y + 6, r: 6, fill: 'none', stroke: o.fg });
    txt(x + 44, y + 10, b, { fill: o.fg, 'font-size': 13 });
    return { b, d };
  });
  // the classical channel and Bob's Pockels cell
  el(svg, 'path', { d: `M${BSA[0] + 100} 62 C 380 40, 520 40, ${PC[0]} 120 L${PC[0]} ${PC[1] - 22}`, fill: 'none', stroke: o.hot, 'stroke-dasharray': '5 4', opacity: .7 });
  txt(390, 40, 'cable: two classical bits', { fill: o.hot });
  el(svg, 'rect', { x: PC[0] - 22, y: PC[1] - 22, width: 44, height: 44, rx: 4, fill: o.bg, stroke: o.fg });
  txt(PC[0], PC[1] + 4, 'Pockels', { fill: o.fg }); txt(PC[0], PC[1] + 15, 'cell', { fill: o.fg });
  const corr = txt(PC[0], PC[1] + 40, '', { fill: o.hot, 'font-size': 14 });

  const dyn = el(svg, 'g');
  const p1 = el(svg, 'circle', { r: 6, fill: o.a }), p2 = el(svg, 'circle', { r: 5, fill: '#b07cff' }), p3 = el(svg, 'circle', { r: 5, fill: '#b07cff' });
  const bits = el(svg, 'text', { 'font-size': 18, fill: o.hot, 'text-anchor': 'middle', 'font-family': o.font, 'font-weight': 700 });
  const caption2 = el(svg, 'text', { x: 380, y: 280, fill: o.muted, 'font-size': 13, 'text-anchor': 'middle', 'font-family': o.font });
  const caption = el(svg, 'text', { x: 380, y: 258, fill: o.fg, 'font-size': 14, 'text-anchor': 'middle', 'font-family': o.font });
  const pathPt = (pts, f) => {                                    // a point along a polyline, f in 0..1
    const L = pts.slice(1).map((q, k) => Math.hypot(q[0] - pts[k][0], q[1] - pts[k][1])), tot = L.reduce((a, b) => a + b, 0);
    let d = f * tot, k = 0;
    while (k < L.length - 1 && d > L[k]) { d -= L[k]; k++; }
    const t = L[k] ? d / L[k] : 0;
    return [lerp(pts[k][0], pts[k + 1][0], t), lerp(pts[k][1], pts[k + 1][1], t)];
  };
  const put = (c, [x, y], on) => { c.setAttribute('cx', x); c.setAttribute('cy', y); c.setAttribute('opacity', on ? 1 : 0); };
  const corrName = (m1, m2) => [m2 ? 'X (H↔V)' : '', m1 ? 'Z (V sign)' : ''].filter(Boolean).join(' then ') || 'no correction';

  let last = -1;
  return loop(6400, runs.length, (i, p) => {
    const r = runs[i], psi = jones(r.input.theta_deg, r.input.phi_deg);
    const before = r.m2 ? gateX(r.m1 ? gateZ(psi) : psi) : (r.m1 ? gateZ(psi) : psi);   // = X^m2 Z^m1 ψ
    const out = jones(r.output.theta_deg, r.output.phi_deg);
    if (last !== i) {                                            // the run's static parts: ellipses and labels
      last = i;
      dyn.replaceChildren();
      polGlyph(dyn, psi, 56, 70, 22, o.a, 'photon 1: ψ', o);
      dyn.dataset.run = i;
    }
    const fIn = seg(p, .02, .25), fBits = seg(p, .32, .62), fP3 = seg(p, .05, .62), fOut = seg(p, .66, .78);
    put(p1, pathPt([P1, BSA], fIn), p < .27);
    put(p2, pathPt([SRC, [BSA[0], BSA[1] + 22]], fIn), p < .27);
    dets.forEach(({ b, d }) => d.setAttribute('fill', p > .27 && b === `${r.m1}${r.m2}` ? o.hot : 'none'));
    const bp = pathPt([[BSA[0] + 100, 62], [390, 46], [PC[0], 118]], fBits);
    bits.setAttribute('x', bp[0]); bits.setAttribute('y', bp[1] - 6); bits.textContent = p > .3 && p < .66 ? `${r.m1}${r.m2}` : '';
    corr.textContent = p > .62 ? corrName(r.m1, r.m2) : '';
    if (p > .6 && !dyn.querySelector('[data-b]')) {
      const g = el(dyn, 'g', { 'data-b': 1 });
      polGlyph(g, before, 530, 110, 18, '#b07cff', 'photon 3 before', o);
    }
    if (p > .78 && !dyn.querySelector('[data-o]')) {      // held to the end of the run
      const g = el(dyn, 'g', { 'data-o': 1 });
      polGlyph(g, out, OUT[0], 100, 22, o.b, 'after: ψ', o);
    }
    // photon 3: through the fibre delay to the Pockels cell, then out corrected
    if (p < .66) put(p3, pathPt([SRC, [416, SRC[1]], [PC[0] - 22, PC[1]]], fP3), p > .02);
    else put(p3, pathPt([[PC[0] + 22, PC[1]], [OUT[0] - 20, OUT[1]]], fOut), p < .8);
    p3.setAttribute('fill', p > .66 ? o.b : '#b07cff');
    caption.textContent = `run ${i + 1} of ${runs.length}: in θ = ${r.input.theta_deg}°, φ = ${r.input.phi_deg}° → out θ = ${r.output.theta_deg}°, φ = ${r.output.phi_deg}°`;
    caption2.textContent = `bits ${r.m1}${r.m2} · ${corrName(r.m1, r.m2)} · fidelity ${r.fidelity}`;
  });
}
