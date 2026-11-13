// Tests a "roll-over shape" in the raw Fukuchi 2017 data (subjects 1-28, 2.5 / 3.5 / 4.5 m/s): does the pressure point under the foot follow a fixed
// curve in a frame attached to the shank (position set by shank lean, not by time)? Also: foot creep relative to the belt, foot flat / heel-off timing,
// and the forward displacement of the hip, mid-pelvis and plate-derived COM in the belt frame.
// Writes tests/data/fukuchi2017-rollover.json. Needs data/ from  node tools/fetch-fukuchi.js --raw .   Run: node tools/rollover-fukuchi.js
// Shares the marker reconstruction, contact detection and force-filter helpers of tools/landing-fukuchi.js (it exports them).
//
// Conventions: x forward, y up (plate = 0), metres; 21 stance fractions 0, 5 ... 100 % between the raw 20 N strike and toe-off; COP-based quantities are NaN where
// Fy < 0.2 BW; Fx, Fy, COPx low-passed at 30 Hz (zero-lag) per contact as in the forceLine block; markers 15 Hz zero-lag as in landing-fukuchi.js.
'use strict';
const fs = require('fs');
const path = require('path');
const R = require('./landing-fukuchi.js');
const { readTable, markers, rigidFit, cond2, cubic, loadStatic, prepTrial, segFilter, THIGH, SHANK, mean, sd, median, fin, r4,
        FS_M, FS_F, DIR, SPEEDS, SPEED_LABEL, NSUBJ, REAR } = R;

const flag = name => { const i = process.argv.indexOf(name); return i > 0 ? process.argv[i + 1] : null; };
const OUT = flag('--roll-out') || path.join(__dirname, '..', 'tests', 'data', 'fukuchi2017-rollover.json');
const G = 9.81, DEG = 180 / Math.PI, FL_LP = 30, MINFY = 0.2;
const PCT = Array.from({ length: 21 }, (_, i) => i * 5);
const FIT_IDX = PCT.map((p, i) => i).filter(i => PCT[i] >= 15 && PCT[i] <= 85);
const VARS = ['lat', 'jc'];   // lat = lateral malleolus / lateral epicondyle markers; jc = ankle and knee joint centres (mid-malleoli, mid-epicondyles)

const lin = (a, u) => { if (!(u >= 0)) return a[0]; const j = Math.min(a.length - 2, Math.floor(u)); return a[j] + (a[j + 1] - a[j]) * (u - j); };
const scale = (p, s) => ({ x: Float64Array.from(p.x, v => v * s), y: Float64Array.from(p.y, v => v * s) });
const cumInt = (a, dt) => { const o = new Float64Array(a.length); for (let i = 1; i < a.length; i++) o[i] = o[i - 1] + 0.5 * (a[i] + a[i - 1]) * dt; return o; };
const angV = (dx, dz) => Math.atan2(dx, dz) * DEG;      // angle from vertical of a vector, forward positive

// ordinary least squares with intercept; r2
function olsFull(x, y) {
  const p = []; for (let i = 0; i < x.length; i++) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) p.push([x[i], y[i]]);
  if (p.length < 6) return null;
  const mx = mean(p.map(q => q[0])), my = mean(p.map(q => q[1])); let sxx = 0, sxy = 0, syy = 0;
  for (const [a, b] of p) { sxx += (a - mx) ** 2; sxy += (a - mx) * (b - my); syy += (b - my) ** 2; }
  const slope = sxy / sxx; return { slope, intercept: my - slope * mx, r2: sxy * sxy / (sxx * syy), n: p.length };
}
// algebraic (Kasa) least-squares circle through points; returns centre, radius and the RMS of the radial residuals
function circleFit(pts) {
  if (pts.length < 8) return null;
  let sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sz = 0, sxz = 0, syz = 0; const n = pts.length;
  for (const [x, y] of pts) { const z = x * x + y * y; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; sz += z; sxz += x * z; syz += y * z; }
  const A = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]], b = [-sxz, -syz, -sz];
  // solve 3x3 by Cramer's rule
  const det = M => M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  const d = det(A); if (Math.abs(d) < 1e-30) return null;
  const col = k => A.map((row, i) => row.map((v, j) => (j === k ? b[i] : v)));
  const D = det(col(0)) / d, E = det(col(1)) / d, F = det(col(2)) / d;
  const cx = -D / 2, cy = -E / 2, r2 = cx * cx + cy * cy - F; if (!(r2 > 0)) return null;
  const r = Math.sqrt(r2); let e2 = 0; for (const [x, y] of pts) e2 += (Math.hypot(x - cx, y - cy) - r) ** 2;
  return { cx, cy, r, rms: Math.sqrt(e2 / n) };
}

function contactForce(T, a, BW) {      // as in momentum-fukuchi.js / the forceLine block
  const i0 = Math.ceil(a.t0 * FS_F), i1 = Math.floor(a.t1 * FS_F), len = i1 - i0 + 1;
  if (len < 10) return null;
  const raw = arr => Array.from(arr.slice(i0, i1 + 1)), fyR = raw(T.force.fy), fxR = raw(T.force.fx), copR = raw(T.force.cop);
  const good = fyR.map(v => v >= MINFY * BW), g0 = good.indexOf(true);
  if (g0 < 0) return null;
  for (let j = g0 - 1; j >= 0; j--) copR[j] = copR[j + 1];
  for (let j = g0 + 1; j < len; j++) if (!good[j]) copR[j] = copR[j - 1];
  return { i0, len, fy: segFilter(fyR, FL_LP), fx: segFilter(fxR, FL_LP), cop: Float64Array.from(segFilter(copR, FL_LP), v => v / 1000) };
}

// ---------------------------------------------------------------- one trial
function prepRoll(id, sp, st) {
  const T = prepTrial(id, sp, st, { rms: { pelvis: [], thigh: [], shank: [] } });
  const mt = readTable(path.join(DIR, `RBDS${id}runT${sp}markers.txt`)), M = markers(mt), n = mt.n, P = {};
  const mid = (A, B) => ({ x: Float64Array.from(A.x, (v, k) => (v + B.x[k]) / 2), y: Float64Array.from(A.y, (v, k) => (v + B.y[k]) / 2) });
  for (const s of ['R', 'L']) {
    const th = rigidFit(M, st.ref, THIGH(s), [`${s}.Knee`, `${s}.Knee.Medial`], n), sk = rigidFit(M, st.ref, SHANK(s), [`${s}.Ankle`, `${s}.Ankle.Medial`], n), L = T.legs[s];
    const mtm = { x: new Float64Array(n).fill(NaN), y: new Float64Array(n).fill(NaN) }, mode = new Int8Array(n).fill(-1);
    for (let k = 0; k < n; k++) {
      const a = L.mt1.x[k], b = L.mt5.x[k];
      if (Number.isFinite(a) && Number.isFinite(b)) { mtm.x[k] = (a + b) / 2; mtm.y[k] = (L.mt1.y[k] + L.mt5.y[k]) / 2; mode[k] = 0; }
      else if (Number.isFinite(a)) { mtm.x[k] = a; mtm.y[k] = L.mt1.y[k]; mode[k] = 1; }
      else if (Number.isFinite(b)) { mtm.x[k] = b; mtm.y[k] = L.mt5.y[k]; mode[k] = 2; }
    }
    const foreOff = st.legs[s].foreOff, off = [foreOff.both, foreOff.mt1, foreOff.mt5];
    const foot = new Float64Array(n);
    for (let k = 0; k < n; k++) foot[k] = mode[k] < 0 ? NaN : Math.atan2(mtm.y[k] - L.heel.y[k], mtm.x[k] - L.heel.x[k]) * DEG - off[mode[k]];
    P[s] = {
      hip: scale(L.hjc, 0.001), kneeLat: scale(L.knee, 0.001), ankleLat: scale(L.ankle, 0.001),
      kneeJc: scale(cond2(mid(th.out[`${s}.Knee`], th.out[`${s}.Knee.Medial`])), 0.001), ankleJc: scale(cond2(mid(sk.out[`${s}.Ankle`], sk.out[`${s}.Ankle.Medial`])), 0.001),
      heel: scale(L.heel, 0.001), mtm: scale(mtm, 0.001), foot
    };
  }
  // belt speed of the trial: median over contacts of the backward speed of the lowest foot marker over 40-60 % of the contact
  const sp_ = [];
  for (const c of T.contacts) {
    const L = T.legs[c.leg], w0 = Math.round((c.t0 + 0.4 * c.dur) * FS_M), w1 = Math.round((c.t0 + 0.6 * c.dur) * FS_M); let best = null;
    for (const nm of ['heel', 'mt1', 'mt5']) {
      const ys = []; for (let k = w0; k <= w1; k++) if (Number.isFinite(L[nm].y[k]) && Number.isFinite(L[nm].x[k])) ys.push(L[nm].y[k]);
      if (ys.length >= 4) { const my = mean(ys); if (!best || my < best.my) best = { nm, my }; }
    }
    if (!best) continue;
    const ts = [], xs = []; for (let k = w0; k <= w1; k++) if (Number.isFinite(L[best.nm].x[k])) { ts.push(k / FS_M); xs.push(L[best.nm].x[k]); }
    const mt_ = mean(ts), mx = mean(xs); let num = 0, den = 0; for (let q = 0; q < ts.length; q++) { num += (ts[q] - mt_) * (xs[q] - mx); den += (ts[q] - mt_) ** 2; }
    sp_.push(-num / den / 1000);
  }
  return { T, n, P, vBelt: median(sp_), fxAll: segFilter(T.force.fx, FL_LP) };
}

// standing shank angle of both variants (static trial means)
function standShank(st, s) {
  const r = st.ref, mid = (a, b) => [(r[a][0] + r[b][0]) / 2, (r[a][1] + r[b][1]) / 2];
  const kl = r[`${s}.Knee`], al = r[`${s}.Ankle`], kj = mid(`${s}.Knee`, `${s}.Knee.Medial`), aj = mid(`${s}.Ankle`, `${s}.Ankle.Medial`);
  return { lat: angV(al[0] - kl[0], kl[1] - al[1]), jc: angV(aj[0] - kj[0], kj[1] - aj[1]) };
}

// ---------------------------------------------------------------- one stride
function strideRec(TR, a, c, leg, subj, ss) {
  const T = TR.T, BW = subj.mass * G, P = TR.P[leg], oth = leg === 'R' ? 'L' : 'R', Po = TR.P[oth], tN = c.t0, dur = a.dur, r = {};
  const cf = contactForce(T, a, BW); if (!cf) return null;
  const at = (arr, t) => lin(arr, Math.min(cf.len - 1, Math.max(0, t * FS_F - cf.i0)));
  const fi = (arr, t) => cubic(arr, t * FS_M), vb = TR.vBelt;
  // plate-derived periodic COM vx over the stride -> displacement
  const jS = Math.ceil(a.t0 * FS_F), jE = Math.floor(tN * FS_F), nj = jE - jS + 1, ax = new Float64Array(nj);
  for (let j = 0; j < nj; j++) ax[j] = TR.fxAll[jS + j] / subj.mass;
  const mAx = mean(Array.from(ax)), vx = cumInt(Float64Array.from(ax, v => v - mAx), 1 / FS_F), mv = mean(Array.from(vx)); for (let j = 0; j < nj; j++) vx[j] -= mv;
  const comDisp = cumInt(vx, 1 / FS_F);
  // belt-frame displacement of a lab x series relative to straight travel at its stride-mean velocity
  const disp = (arr, t) => { const x0 = fi(arr, a.t0), xN = fi(arr, tN), vbar = (xN - x0) / (tN - a.t0); return fi(arr, t) - x0 - vbar * (t - a.t0); };
  const midx = Float64Array.from(P.hip.x, (v, k) => (v + Po.hip.x[k]) / 2);
  // level of heel / forefoot height in mid-stance (40-60 %)
  const lev = arr => median(fin([0.4, 0.45, 0.5, 0.55, 0.6].map(f => fi(arr, a.t0 + f * dur))));
  const heelLev = lev(P.heel.y), foreLev = lev(P.mtm.y);
  // foot flat: first stance fraction (1 % grid) within 2 deg of the 40-50 % footDeg
  const ref = median(fin([0.4, 0.45, 0.5].map(f => fi(P.foot, a.t0 + f * dur)))); r.footFlatPct = NaN;
  for (let p = 0; p <= 100; p++) { const v = fi(P.foot, a.t0 + p / 100 * dur); if (Number.isFinite(v) && Math.abs(v - ref) <= 2) { r.footFlatPct = p; break; } }
  // heel-off: Heel.Bottom first 10 mm above its mid-stance level
  r.heelOffPct = NaN;
  for (let k = Math.ceil((a.t0 + 0.5 * dur) * FS_M); k / FS_M <= a.t1 + 0.02; k++) {
    const y1 = P.heel.y[k], y0 = P.heel.y[k - 1], lv = heelLev + 0.010;
    if (Number.isFinite(y1) && Number.isFinite(y0) && y0 <= lv && y1 > lv) { r.heelOffPct = 100 * ((k - 1 + (lv - y0) / (y1 - y0)) / FS_M - a.t0) / dur; break; }
  }
  const crp = (arr, t) => fi(arr, t) + vb * (t - a.t0);       // marker x in the belt frame (belt moves backward at vb)
  const c30 = { ank: crp(P.ankleLat.x, a.t0 + 0.3 * dur), heel: crp(P.heel.x, a.t0 + 0.3 * dur), fore: crp(P.mtm.x, a.t0 + 0.3 * dur) };
  PCT.forEach((p, i) => {
    const t = a.t0 + p / 100 * dur, fy = at(cf.fy, t), cop = at(cf.cop, t), valid = fy >= MINFY * BW;
    const hx = fi(P.hip.x, t), hy = fi(P.hip.y, t);
    r[`c_footDeg_${i}`] = fi(P.foot, t);
    const thighJc = angV(fi(P.kneeJc.x, t) - hx, hy - fi(P.kneeJc.y, t));
    for (const v of VARS) {
      const kx = fi(v === 'lat' ? P.kneeLat.x : P.kneeJc.x, t), ky = fi(v === 'lat' ? P.kneeLat.y : P.kneeJc.y, t), axx = fi(v === 'lat' ? P.ankleLat.x : P.ankleJc.x, t), ay = fi(v === 'lat' ? P.ankleLat.y : P.ankleJc.y, t);
      const shank = angV(axx - kx, ky - ay);
      r[`c_shankDeg_${v}_${i}`] = shank; r[`c_ankleDF_${v}_${i}`] = -(shank - ss[v]) + r[`c_footDeg_${i}`];
      const L_ = Math.hypot(kx - axx, ky - ay), ux = (kx - axx) / L_, uy = (ky - ay) / L_, px = uy, py = -ux, dx = cop - axx, dy = 0 - ay;
      r[`c_copAheadAnkle_${v}_${i}`] = valid ? cop - axx : NaN;
      r[`c_perp_${v}_${i}`] = valid ? dx * px + dy * py : NaN; r[`c_along_${v}_${i}`] = valid ? dx * ux + dy * uy : NaN;
    }
    r[`c_thighDeg_${i}`] = thighJc; r[`c_kneeDeg_${i}`] = thighJc - r[`c_shankDeg_jc_${i}`];
    r[`c_creepAnkle_${i}`] = crp(P.ankleLat.x, t) - c30.ank; r[`c_creepHeel_${i}`] = crp(P.heel.x, t) - c30.heel; r[`c_creepFore_${i}`] = crp(P.mtm.x, t) - c30.fore;
    r[`c_heelH_${i}`] = fi(P.heel.y, t) - heelLev; r[`c_foreH_${i}`] = fi(P.mtm.y, t) - foreLev;
    const dh = disp(P.hip.x, t), dm = disp(midx, t), dc = lin(comDisp, Math.min(nj - 1, Math.max(0, t * FS_F - jS)));
    r[`c_hipDisp_${i}`] = dh; r[`c_midDisp_${i}`] = dm; r[`c_comDisp_${i}`] = dc;
    r[`c_hipMinusMid_${i}`] = dh - dm; r[`c_midMinusCom_${i}`] = dm - dc; r[`c_hipMinusCom_${i}`] = dh - dc;
  });
  // creep scalars (m): strike -> 30 % and heel-off -> toe-off
  for (const [k, nm] of [['Ankle', 'creepAnkle'], ['Heel', 'creepHeel'], ['Fore', 'creepFore']]) {
    r[`creep_0_30_${k}`] = r[`c_${nm}_6`] - r[`c_${nm}_0`];
    const arr = k === 'Ankle' ? P.ankleLat.x : k === 'Heel' ? P.heel.x : P.mtm.x;
    r[`creep_ho_to_${k}`] = Number.isFinite(r.heelOffPct) ? crp(arr, a.t1) - crp(arr, a.t0 + r.heelOffPct / 100 * dur) : NaN;
  }
  r.heelH0 = r.c_heelH_0; r.foreH0 = r.c_foreH_0; r.vBelt = vb;
  return r;
}

// ---------------------------------------------------------------- main
function main() {
  const infoLines = fs.readFileSync(path.join(DIR, 'RBDSinfo.txt'), 'utf8').split(/\r?\n/).filter(Boolean);
  const head = infoLines[0].split('\t'), ix = n => head.indexOf(n), info = {};
  for (const l of infoLines.slice(1)) { const f = l.split('\t'), id = parseInt(f[ix('Subject')], 10); if (id >= 1 && id <= NSUBJ && !info[id]) info[id] = f; }
  const ids = Object.keys(info).map(Number).sort((a, b) => a - b);
  const heights = ids.map(i => parseFloat(info[i][ix('Height')]) / 100), masses = ids.map(i => parseFloat(info[i][ix('Mass')]));
  const men = ids.filter(i => info[i][ix('Gender')].trim() === 'M').length;

  const rows = [], rej = { timing: 0, nostride: 0, few: 0 }, legLen = {};
  for (const sid of ids) {
    const id = String(sid).padStart(3, '0'), k = ids.indexOf(sid), subj = { height: heights[k], mass: masses[k] }, st = loadStatic(id);
    legLen[sid] = mean(['R', 'L'].map(s => { const r = st.ref, d = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]); return (d(r[`${s}.GTR`], r[`${s}.Knee`]) + d(r[`${s}.Knee`], r[`${s}.Ankle`])) / 1000; }));
    for (const sp of SPEEDS) {
      const TR = prepRoll(id, sp, st), cs = TR.T.contacts, legRecs = {}, fsOf = {};
      const steps = []; for (let i = 0; i + 1 < cs.length; i++) steps.push(cs[i + 1].t0 - cs[i].t0);
      const medStep = median(steps.filter(x => x > 0.2 && x < 1.0));
      for (const s of ['R', 'L']) {
        fsOf[s] = (info[sid][ix(s + 'FSI' + sp)] || '').trim();
        const recs = [], ss = standShank(st, s);
        for (let i = 0; i + 2 < cs.length; i++) {
          const a = cs[i], b = cs[i + 1], c = cs[i + 2];
          if (a.leg !== s || b.leg === s || c.leg !== s) continue;
          if (Math.abs((b.t0 - a.t0) / medStep - 1) > 0.3 || Math.abs((c.t0 - b.t0) / medStep - 1) > 0.3) { rej.timing++; continue; }
          const r = strideRec(TR, a, c, s, subj, ss); if (r) recs.push(r); else rej.nostride++;
        }
        if (recs.length < 5) { rej.few++; legRecs[s] = null; continue; }
        const rec = { nStrides: recs.length };
        for (const key in recs[0]) { const v = fin(recs.map(q => q[key])); rec[key] = v.length ? mean(v) : NaN; }
        legRecs[s] = rec;
      }
      rows.push({ id: sid, sp, fsOf, legRecs, height: subj.height });
    }
    process.stderr.write(`\rsubject ${sid}/${NSUBJ}`);
  }
  process.stderr.write('\n');

  // ------------------------------------------------ cohorts and subject-level fits
  const avgLegs = legs => { const out = { nLegs: legs.length }; for (const k in legs[0]) { const v = fin(legs.map(l => l[k])); out[k] = v.length ? mean(v) : NaN; } return out; };
  const cohort = { all: { 25: [], 35: [], 45: [] }, rearfoot: { 25: [], 35: [], 45: [] } };
  for (const r of rows) {
    const all = ['R', 'L'].filter(s => r.legRecs[s]).map(s => r.legRecs[s]), rear = ['R', 'L'].filter(s => r.legRecs[s] && r.fsOf[s] === REAR).map(s => r.legRecs[s]);
    if (all.length) cohort.all[r.sp].push({ id: r.id, height: r.height, v: avgLegs(all) });
    if (rear.length) cohort.rearfoot[r.sp].push({ id: r.id, height: r.height, v: avgLegs(rear) });
  }
  const sinOf = d => Math.sin(-d / DEG), radOf = d => -d / DEG;
  const pts = (v, vr, idx) => idx.map(i => [v[`c_perp_${vr}_${i}`], v[`c_along_${vr}_${i}`]]).filter(p => Number.isFinite(p[0]) && Number.isFinite(p[1]));
  for (const c in cohort) {
    for (const sp of SPEEDS) for (const row of cohort[c][sp]) {
      const v = row.v;
      for (const vr of VARS) {
        const f = circleFit(pts(v, vr, FIT_IDX));
        v[`circ_cx_${vr}`] = f ? f.cx : NaN; v[`circ_cy_${vr}`] = f ? f.cy : NaN; v[`circ_r_${vr}`] = f ? f.r : NaN; v[`circ_rms_${vr}`] = f ? f.rms : NaN;
        v[`circ_rOverHeight_${vr}`] = f ? f.r / row.height : NaN; v[`circ_rOverLeg_${vr}`] = f ? f.r / legLen[row.id] : NaN;
        for (const [nm, fn] of [['sin', sinOf], ['rad', radOf]]) {
          const g = olsFull(FIT_IDX.map(i => fn(v[`c_shankDeg_${vr}_${i}`])), FIT_IDX.map(i => v[`c_copAheadAnkle_${vr}_${i}`]));
          v[`reg_${nm}_slope_${vr}`] = g ? g.slope : NaN; v[`reg_${nm}_int_${vr}`] = g ? g.intercept : NaN; v[`reg_${nm}_r2_${vr}`] = g ? g.r2 : NaN;
        }
      }
    }
  }
  // pooled over speeds per subject (subjects present at all three speeds in the cohort)
  const pooled = { all: [], rearfoot: [] };
  for (const c in cohort) {
    const byId = {}; for (const sp of SPEEDS) for (const row of cohort[c][sp]) (byId[row.id] = byId[row.id] || {})[sp] = row;
    for (const id in byId) {
      if (!SPEEDS.every(sp => byId[id][sp])) continue;
      const o = { id: +id };
      for (const vr of VARS) for (const [nm, fn] of [['sin', sinOf], ['rad', radOf]]) {
        const x = [], y = [];
        for (const sp of SPEEDS) for (const i of FIT_IDX) { x.push(fn(byId[id][sp].v[`c_shankDeg_${vr}_${i}`])); y.push(byId[id][sp].v[`c_copAheadAnkle_${vr}_${i}`]); }
        const g = olsFull(x, y); o[`${nm}_${vr}`] = g;
      }
      pooled[c].push(o);
    }
  }

  const KEYS = Object.keys(cohort.all[25][0].v).filter(k => k !== 'nLegs');
  const scalarKeys = KEYS.filter(k => !/^c_/.test(k)), curveNames = [...new Set(KEYS.filter(k => /^c_/.test(k)).map(k => k.replace(/^c_(.*)_\d+$/, '$1')))];
  const stat = vals => { const v = fin(vals); return { mean: r4(mean(v)), sd: r4(sd(v)), n: v.length }; };
  const summary = {
    source: { citation: 'Fukuchi RK, Fukuchi CA, Duarte M (2017). A public dataset of running biomechanics and the effects of running speed on lower extremity kinematics and kinetics. PeerJ 5:e3298',
              doi: '10.7717/peerj.3298', figshare: 'https://doi.org/10.6084/m9.figshare.4543435 (article 4543435)', licence: 'CC BY 4.0', note: 'derived summary statistics from the raw marker and force files' },
    generated: new Date().toISOString(),
    method: {
      frame: 'x forward, y up (plate = 0), metres; stance fractions 0, 5 ... 100 % between the raw 20 N strike and toe-off; COP-based values masked where Fy < 0.2 BW; Fx, Fy, COPx 30 Hz zero-lag per contact (COP held at low force before filtering), markers 15 Hz zero-lag (as landing-fukuchi.js / the forceLine block)',
      variants: 'lat: lateral malleolus (static Ankle carried by the shank-cluster fit) and lateral epicondyle (Knee); jc: ankle joint centre = midpoint of Ankle and Ankle.Medial, knee joint centre = midpoint of Knee and Knee.Medial (all carried by the cluster fits)',
      shankDeg: 'angle from vertical of the knee-to-ankle line, positive when the ankle is ahead of the knee (shank leaning backward); shank lean forward = -shankDeg; thighDeg = hip joint centre (Harrington) to knee joint centre; kneeDeg = thighDeg - shankDeg_jc',
      footDeg: 'Heel.Bottom to forefoot (MT1/MT5 midpoint, or the single marker available) line above horizontal minus the standing value with the same marker choice',
      ankleDF: 'dorsiflexion change from standing = (-(shankDeg - standing shankDeg)) + footDeg: forward shank lean plus toe-up foot rotation (shank vertical and foot flat in standing is the zero)',
      shankFrame: 'origin at the ankle point; u = unit vector from ankle to knee (along), p = u rotated by -90 deg (perpendicular; forward when the shank is vertical); the pressure point is (COPx, 0); perp = d.p, along = d.u',
      circle: 'algebraic (Kasa) least-squares circle through the leg-averaged subject trace (perp, along) at 15-85 % of stance (13 points); rms = rms radial residual; radius also / subject height and / leg length (|GTR-Knee| + |Knee-Ankle| standing, 3D)',
      regression: 'per subject, copAheadAnkle on sin(-shankDeg) (slope = effective radius, m) and on -shankDeg in radians, over 15-85 % of stance, on the leg-averaged curves; pooled = the 3 speeds together per subject',
      creep: 'x of a marker + v_belt (t - t_strike) relative to its value at 30 % of stance; v_belt = trial median of the backward speed of the lowest foot marker in mid-stance; creep_0_30 = change from strike to 30 % (forward positive); creep_ho_to = change from heel-off to toe-off',
      footFlat: 'first stance fraction (1 % grid) at which footDeg is within 2 deg of its median at 40, 45, 50 %; heelOff = Heel.Bottom first 10 mm above its 40-60 % median, after 50 %',
      heights: 'heelH / foreH = Heel.Bottom / forefoot-midpoint height above its median over 40-60 % of stance (m)',
      displacement: 'hipDisp / midDisp / comDisp: forward displacement since foot strike relative to straight travel at the stride-mean velocity of the same point (hip = stance-side hip joint centre, mid = midpoint of the two hip joint centres, com = plate-derived COM: integral of the periodic velocity from the stride Fx/m, zero mean velocity); hipMinusMid, midMinusCom, hipMinusCom are differences',
      cohorts: 'strides -> legs -> subjects as in landing-fukuchi.js; rearfoot = legs whose RFSI/LFSI is "Rearfoot"; subject-level fits are computed on the leg-averaged curves, then averaged across subjects'
    },
    subjects: { n: ids.length, men, women: ids.length - men }, cohorts: { all: {}, rearfoot: {} }, pooledOverSpeeds: {}
  };
  for (const c in cohort) for (const sp of SPEEDS) {
    const rs = cohort[c][sp], e = { n: rs.length, legs: rs.reduce((s, r) => s + r.v.nLegs, 0), scalars: {}, stanceCurves: {} };
    for (const k of scalarKeys) e.scalars[k] = stat(rs.map(r => r.v[k]));
    for (const nm of curveNames) { const o = { mean: [], sd: [], n: [] }; PCT.forEach((p, i) => { const s = stat(rs.map(r => r.v[`c_${nm}_${i}`])); o.mean.push(s.mean); o.sd.push(s.sd); o.n.push(s.n); }); e.stanceCurves[nm] = o; }
    summary.cohorts[c][SPEED_LABEL[sp]] = e;
  }
  for (const c in pooled) { summary.pooledOverSpeeds[c] = { n: pooled[c].length }; for (const vr of VARS) for (const nm of ['sin', 'rad']) for (const q of ['slope', 'intercept', 'r2']) summary.pooledOverSpeeds[c][`${nm}_${q}_${vr}`] = stat(pooled[c].map(o => o[`${nm}_${vr}`] ? o[`${nm}_${vr}`][q] : NaN)); }
  summary.rejected = rej;
  fs.mkdirSync(path.dirname(OUT), { recursive: true });
  fs.writeFileSync(OUT, JSON.stringify(summary));

  // ------------------------------------------------ report
  const out = [], P = s => out.push(s), E = (c, sp) => summary.cohorts[c][SPEED_LABEL[sp]];
  const sc = (c, sp, k, d = 3, m = 1) => { const s = E(c, sp).scalars[k]; return s && s.mean !== null ? `${(s.mean * m).toFixed(d)} ± ${s.sd === null ? '-' : (s.sd * m).toFixed(d)}` : '-'; };
  P(`subjects ${ids.length}; strides rejected ${JSON.stringify(rej)}; measured belt speed (trial median, m/s): ` + SPEEDS.map(sp => `${SPEED_LABEL[sp]}: ${sc('all', sp, 'vBelt', 3)}`).join(', '));
  for (const c of ['all', 'rearfoot']) {
    P(`\n===== cohort ${c}: scalars (mean ± SD across subjects; n = ${SPEEDS.map(sp => E(c, sp).n).join('/')}) =====`);
    P('metric'.padEnd(46) + SPEEDS.map(sp => (SPEED_LABEL[sp] + ' m/s').padStart(20)).join(''));
    const row = (lab, k, d = 3, m = 1) => P(lab.padEnd(46) + SPEEDS.map(sp => sc(c, sp, k, d, m).padStart(20)).join(''));
    for (const vr of VARS) {
      P(`-- circle fit in the shank frame, ${vr}, 15-85 % of stance`);
      row('  centre perpendicular (m)', `circ_cx_${vr}`); row('  centre along shank axis (m)', `circ_cy_${vr}`); row('  radius (m)', `circ_r_${vr}`); row('  radius / height', `circ_rOverHeight_${vr}`); row('  radius / leg length', `circ_rOverLeg_${vr}`); row('  rms residual (mm)', `circ_rms_${vr}`, 1, 1000);
      P(`-- copAheadAnkle regression, ${vr} (per subject, 15-85 %)`);
      row('  sin(-shankDeg): slope = radius (m)', `reg_sin_slope_${vr}`); row('  sin(-shankDeg): intercept (m)', `reg_sin_int_${vr}`); row('  sin(-shankDeg): r2', `reg_sin_r2_${vr}`, 3);
      row('  -shankDeg (rad): slope (m/rad)', `reg_rad_slope_${vr}`); row('  -shankDeg (rad): intercept (m)', `reg_rad_int_${vr}`); row('  -shankDeg (rad): r2', `reg_rad_r2_${vr}`, 3);
    }
    P(`-- pooled over the three speeds (per subject, n = ${summary.pooledOverSpeeds[c].n}); mean ± SD`);
    for (const vr of VARS) for (const nm of ['sin', 'rad']) { const g = q => { const s = summary.pooledOverSpeeds[c][`${nm}_${q}_${vr}`]; return `${s.mean.toFixed(3)} ± ${s.sd.toFixed(3)}`; }; P(`  ${vr} ${nm}: slope ${g('slope')}, intercept ${g('intercept')}, r2 ${g('r2')}`); }
    P('-- creep relative to the belt (mm)');
    for (const nm of ['Ankle', 'Heel', 'Fore']) { row(`  ${nm} strike -> 30 % of stance`, `creep_0_30_${nm}`, 1, 1000); row(`  ${nm} heel-off -> toe-off`, `creep_ho_to_${nm}`, 1, 1000); }
    row('heel height above mid-stance at strike (mm)', 'heelH0', 1, 1000); row('forefoot height above mid-stance at strike (mm)', 'foreH0', 1, 1000);
    P('-- timing (% of stance)'); row('  foot flat (footDeg within 2 deg of 40-50 %)', 'footFlatPct', 1); row('  heel-off (heel +10 mm)', 'heelOffPct', 1);
  }
  const ct = (c, sp, list, every) => {
    const e = E(c, sp).stanceCurves;
    P(''.padEnd(34) + PCT.filter((p, i) => i % every === 0).map(p => String(p).padStart(12)).join(''));
    for (const [k, lab, d, m] of list) { const o = e[k]; P(lab.padEnd(34) + o.mean.map((v, i) => i % every === 0 ? ((v === null ? '-' : `${(v * m).toFixed(d)}±${o.sd[i] === null ? '-' : (o.sd[i] * m).toFixed(d)}`).padStart(12)) : '').join('')); }
  };
  for (const c of ['all', 'rearfoot']) for (const sp of SPEEDS) {
    P(`\n===== curves, cohort ${c}, ${SPEED_LABEL[sp]} m/s (% of stance, every 10 %) =====`);
    ct(c, sp, [['shankDeg_jc', 'shankDeg jc (deg)', 1, 1], ['shankDeg_lat', 'shankDeg lat (deg)', 1, 1], ['footDeg', 'footDeg (deg)', 1, 1], ['thighDeg', 'thighDeg (deg)', 1, 1], ['kneeDeg', 'kneeDeg (deg)', 1, 1], ['ankleDF_jc', 'ankle DF jc (deg)', 1, 1], ['ankleDF_lat', 'ankle DF lat (deg)', 1, 1],
      ['copAheadAnkle_lat', 'copAheadAnkle lat (m)', 3, 1], ['copAheadAnkle_jc', 'copAheadAnkle jc (m)', 3, 1], ['perp_lat', 'trace perp lat (m)', 3, 1], ['along_lat', 'trace along lat (m)', 3, 1], ['perp_jc', 'trace perp jc (m)', 3, 1], ['along_jc', 'trace along jc (m)', 3, 1],
      ['creepAnkle', 'ankle creep vs belt (mm)', 1, 1000], ['creepHeel', 'heel creep vs belt (mm)', 1, 1000], ['creepFore', 'forefoot creep vs belt (mm)', 1, 1000], ['heelH', 'heel height vs mid-stance (mm)', 1, 1000], ['foreH', 'forefoot height vs mid-stance (mm)', 1, 1000],
      ['hipDisp', 'hip disp. in belt frame (m)', 3, 1], ['midDisp', 'mid-pelvis disp. (m)', 3, 1], ['comDisp', 'plate COM disp. (m)', 3, 1], ['hipMinusMid', 'hip - mid-pelvis (m)', 3, 1], ['midMinusCom', 'mid-pelvis - COM (m)', 3, 1], ['hipMinusCom', 'hip - COM (m)', 3, 1]], 2);
  }
  console.log(out.join('\n'));
  console.log('\nwrote ' + path.relative(process.cwd(), OUT));
}

main();
