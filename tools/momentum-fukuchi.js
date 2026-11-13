// Tests whole-body mechanics against the measured ground force in the raw Fukuchi 2017 data (subjects 1-28, 2.5 / 3.5 / 4.5 m/s):
// a seven-mass sagittal segment model (thighs, shanks, feet, head-arms-trunk point mass) is compared with the force plate.
// Writes tests/data/fukuchi2017-momentum.json. Needs data/ from  node tools/fetch-fukuchi.js --raw .   Run: node tools/momentum-fukuchi.js
// Shares the marker reconstruction, contact detection and force-filter helpers of tools/landing-fukuchi.js (it exports them).
//
// Conventions (sagittal plane): x forward, y up (plate = 0), metres; angles anticlockwise positive.
//   Angular momentum H (per kg of body mass) is the z component of  sum_i f_i [ (r_i - r_com) x (v_i - v_com) ] + sum_i f_i k_i^2 L_i^2 w_i
//   with w = d(angle)/dt, angle = atan2(dy, dx) of the segment axis (anticlockwise positive).
//   Torque of the ground force about the COM = (x_cop - x_com) Fy + y_com Fx = m dH/dt  (checked numerically below).
'use strict';
const fs = require('fs');
const path = require('path');
const R = require('./landing-fukuchi.js');
const { readTable, markers, rigidFit, cond2, conditioned, cubic, loadStatic, prepTrial, segFilter, THIGH, SHANK, PELVIS, mean, sd, median, fin, r4, ols,
        FS_M, FS_F, DIR, SPEEDS, SPEED_LABEL, SPEED_MS, NSUBJ, REAR } = R;

const flag = name => { const i = process.argv.indexOf(name); return i > 0 ? process.argv[i + 1] : null; };
const OUT = flag('--mom-out') || path.join(__dirname, '..', 'tests', 'data', 'fukuchi2017-momentum.json');
const G = 9.81, DEG = 180 / Math.PI, FL_LP = 30, MINFY = 0.2;
const PCT = Array.from({ length: 21 }, (_, i) => i * 5);
// de Leva 1996 (male): mass fraction, COM from proximal joint (fraction of length), radius of gyration (fraction of length)
const SEG = { thigh: { f: 0.1416, com: 0.4095, k: 0.329 }, shank: { f: 0.0433, com: 0.4459, k: 0.255 }, foot: { f: 0.0137, com: 0.44, k: 0.257 } };
const F_HAT = 1 - 2 * (SEG.thigh.f + SEG.shank.f + SEG.foot.f), HAT_H = 0.180;
const FIX_LEANS = [0, 5, 10];

// ---------------------------------------------------------------- helpers
const diff1 = a => { const o = new Float64Array(a.length).fill(NaN); for (let k = 1; k < a.length - 1; k++) o[k] = (a[k + 1] - a[k - 1]) * FS_M / 2; return o; };
const diff2 = a => { const o = new Float64Array(a.length).fill(NaN); for (let k = 1; k < a.length - 1; k++) o[k] = (a[k + 1] - 2 * a[k] + a[k - 1]) * FS_M * FS_M; return o; };
const cumInt = (a, dt) => { const o = new Float64Array(a.length); for (let i = 1; i < a.length; i++) o[i] = o[i - 1] + 0.5 * (a[i] + a[i - 1]) * dt; return o; };
const rmsOf = d => { const v = fin(d); return v.length ? Math.sqrt(mean(v.map(x => x * x))) : NaN; };
function corr(x, y) {
  const p = []; for (let i = 0; i < x.length; i++) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) p.push([x[i], y[i]]);
  if (p.length < 5) return NaN;
  const mx = mean(p.map(q => q[0])), my = mean(p.map(q => q[1])); let a = 0, b = 0, c = 0;
  for (const [u, v] of p) { a += (u - mx) * (v - my); b += (u - mx) ** 2; c += (v - my) ** 2; }
  return a / Math.sqrt(b * c);
}
const rangeOf = a => { const v = fin(Array.from(a)); return v.length ? Math.max(...v) - Math.min(...v) : NaN; };
const lin = (a, u) => { if (!(u >= 0)) return a[0]; const j = Math.min(a.length - 2, Math.floor(u)); return a[j] + (a[j + 1] - a[j]) * (u - j); };
const scale = (p, s) => ({ x: Float64Array.from(p.x, v => v * s), y: Float64Array.from(p.y, v => v * s) });

// ---------------------------------------------------------------- one trial: segment arrays (metres), forces, pelvis angles
function prepMomentum(id, sp, st, log) {
  const T = prepTrial(id, sp, st, { rms: { pelvis: [], thigh: [], shank: [] } });
  const mt = readTable(path.join(DIR, `RBDS${id}runT${sp}markers.txt`)), M = markers(mt), n = mt.n;
  const P = {}, S = {};
  for (const s of ['R', 'L']) {
    for (const nm of ['Knee', 'Knee.Medial', 'Ankle', 'Ankle.Medial']) if (!st.ref[`${s}.${nm}`]) throw new Error(`static ${s}.${nm} missing`);
    const th = rigidFit(M, st.ref, THIGH(s), [`${s}.Knee`, `${s}.Knee.Medial`], n), sk = rigidFit(M, st.ref, SHANK(s), [`${s}.Ankle`, `${s}.Ankle.Medial`], n);
    const mid = (A, B) => ({ x: Float64Array.from(A.x, (v, k) => (v + B.x[k]) / 2), y: Float64Array.from(A.y, (v, k) => (v + B.y[k]) / 2) });
    const L = T.legs[s];
    const mtm = { x: new Float64Array(n).fill(NaN), y: new Float64Array(n).fill(NaN) };
    for (let k = 0; k < n; k++) {      // forefoot point: MT1/MT5 midpoint, or the single one available
      const xs = fin([L.mt1.x[k], L.mt5.x[k]]), ys = fin([L.mt1.y[k], L.mt5.y[k]]);
      if (xs.length) { mtm.x[k] = mean(xs); mtm.y[k] = mean(ys); }
    }
    P[s] = {
      hip: scale(L.hjc, 0.001), knee: scale(cond2(mid(th.out[`${s}.Knee`], th.out[`${s}.Knee.Medial`])), 0.001),
      ankle: scale(cond2(mid(sk.out[`${s}.Ankle`], sk.out[`${s}.Ankle.Medial`])), 0.001), heel: scale(L.heel, 0.001), mtm: scale(mtm, 0.001)
    };
    // segments: centre, axis angle, length
    const seg = (A, B, fc, key) => {
      const cx = new Float64Array(n), cy = new Float64Array(n), ang = new Float64Array(n), len = new Float64Array(n);
      for (let k = 0; k < n; k++) {
        const dx = B.x[k] - A.x[k], dy = B.y[k] - A.y[k];
        cx[k] = A.x[k] + fc * dx; cy[k] = A.y[k] + fc * dy; ang[k] = Math.atan2(dy, dx); len[k] = Math.hypot(dx, dy);
      }
      const I = Float64Array.from(len, l => SEG[key].f * (SEG[key].k * l) ** 2);
      return { f: SEG[key].f, cx, cy, w: diff1(ang), I, vx: diff1(cx), vy: diff1(cy) };
    };
    const p = P[s];
    // foot: toe tip = heel + 1.5 (forefoot midpoint - heel); COM 0.44 of the way heel -> toe; length heel -> toe
    const toe = { x: Float64Array.from(p.heel.x, (v, k) => v + 1.5 * (p.mtm.x[k] - v)), y: Float64Array.from(p.heel.y, (v, k) => v + 1.5 * (p.mtm.y[k] - v)) };
    S[s] = { thigh: seg(p.hip, p.knee, SEG.thigh.com, 'thigh'), shank: seg(p.knee, p.ankle, SEG.shank.com, 'shank'), foot: seg(p.heel, toe, SEG.foot.com, 'foot') };
  }
  // pelvis angles from the rigidly fitted pelvis markers
  const pel = rigidFit(M, st.ref, PELVIS, ['R.ASIS', 'L.ASIS', 'R.PSIS', 'L.PSIS'], n);
  const tilt = new Float64Array(n).fill(NaN), rot = new Float64Array(n).fill(NaN);
  for (let k = 0; k < n; k++) {
    const g = (nm, c) => pel.out[nm][c][k];
    const ax = (g('R.ASIS', 'x') + g('L.ASIS', 'x')) / 2, ay = (g('R.ASIS', 'y') + g('L.ASIS', 'y')) / 2, px = (g('R.PSIS', 'x') + g('L.PSIS', 'x')) / 2, py = (g('R.PSIS', 'y') + g('L.PSIS', 'y')) / 2;
    tilt[k] = Math.atan2(-(ay - py), ax - px) * DEG;                                                    // + = anterior tilt (ASIS below PSIS)
    rot[k] = Math.atan2(g('R.ASIS', 'x') - g('L.ASIS', 'x'), g('R.ASIS', 'z') - g('L.ASIS', 'z')) * DEG;   // + = right ASIS ahead of left (right hip forward)
  }
  return { T, n, P, S, tilt: conditioned(tilt), rot: conditioned(rot), fyAll: segFilter(T.force.fy, FL_LP), fxAll: segFilter(T.force.fx, FL_LP) };
}

// whole-body COM, angular momentum and its parts for a given trunk lean (per kg of body mass)
function bodyModel(TR, height, leanDeg) {
  const n = TR.n, dy = HAT_H * height, dx = dy * Math.tan(leanDeg / DEG);
  const midx = new Float64Array(n), midy = new Float64Array(n);
  for (let k = 0; k < n; k++) { midx[k] = (TR.P.R.hip.x[k] + TR.P.L.hip.x[k]) / 2; midy[k] = (TR.P.R.hip.y[k] + TR.P.L.hip.y[k]) / 2; }
  const hx = Float64Array.from(midx, v => v + dx), hy = Float64Array.from(midy, v => v + dy), hvx = diff1(hx), hvy = diff1(hy);
  const comx = new Float64Array(n), comy = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    let X = F_HAT * hx[k], Y = F_HAT * hy[k];
    for (const s of ['R', 'L']) for (const g of ['thigh', 'shank', 'foot']) { const q = TR.S[s][g]; X += q.f * q.cx[k]; Y += q.f * q.cy[k]; }
    comx[k] = X; comy[k] = Y;
  }
  const vcx = diff1(comx), vcy = diff1(comy), acx = diff2(comx), acy = diff2(comy);
  const part = { R: new Float64Array(n), L: new Float64Array(n), H: new Float64Array(n) };
  for (let k = 0; k < n; k++) {
    for (const s of ['R', 'L']) { let h = 0; for (const g of ['thigh', 'shank', 'foot']) { const q = TR.S[s][g]; h += q.f * ((q.cx[k] - comx[k]) * (q.vy[k] - vcy[k]) - (q.cy[k] - comy[k]) * (q.vx[k] - vcx[k])) + q.I[k] * q.w[k]; } part[s][k] = h; }
    part.H[k] = F_HAT * ((hx[k] - comx[k]) * (hvy[k] - vcy[k]) - (hy[k] - comy[k]) * (hvx[k] - vcx[k]));
  }
  const H = Float64Array.from(part.R, (v, k) => v + part.L[k] + part.H[k]);
  return { comx, comy, vcx, vcy, acx, acy, H, dH: diff1(H), dpart: { R: diff1(part.R), L: diff1(part.L), H: diff1(part.H) }, midx, midy, dyHat: dy };
}

// filtered force of one contact (as in the forceLine block) -> {i0, len, fy, fx, cop (m)} ; COP held where Fy < 0.2 BW before filtering
function contactForce(T, a, BW) {
  const i0 = Math.ceil(a.t0 * FS_F), i1 = Math.floor(a.t1 * FS_F), len = i1 - i0 + 1;
  if (len < 10) return null;
  const raw = arr => Array.from(arr.slice(i0, i1 + 1)), fyR = raw(T.force.fy), fxR = raw(T.force.fx), copR = raw(T.force.cop);
  const good = fyR.map(v => v >= MINFY * BW), g0 = good.indexOf(true);
  if (g0 < 0) return null;
  for (let j = g0 - 1; j >= 0; j--) copR[j] = copR[j + 1];
  for (let j = g0 + 1; j < len; j++) if (!good[j]) copR[j] = copR[j - 1];
  return { i0, len, fy: segFilter(fyR, FL_LP), fx: segFilter(fxR, FL_LP), cop: Float64Array.from(segFilter(copR, FL_LP), v => v / 1000) };
}

// ---------------------------------------------------------------- one stride
function strideRec(TR, models, variants, a, c, leg, subj, acc) {
  const m = subj.mass, BW = m * G, T = TR.T, r = {};
  const oth = leg === 'R' ? 'L' : 'R', dur = a.dur, tN = c.t0;
  const cf = contactForce(T, a, BW);
  if (!cf) return null;
  const at = (arr, t) => lin(arr, Math.min(cf.len - 1, Math.max(0, t * FS_F - cf.i0)));
  const fi = (arr, t) => cubic(arr, t * FS_M);
  const kS0 = Math.ceil(a.t0 * FS_M), kS1 = Math.floor(a.t1 * FS_M), k0 = kS0, k1 = Math.floor(tN * FS_M);
  const hip = TR.P[leg].hip, m0 = models[0];
  const hipvx = diff1(hip.x);
  // ---- plate-derived series over the stride (300 Hz, 30 Hz zero-lag)
  const jS = Math.ceil(a.t0 * FS_F), jE = Math.floor(tN * FS_F), nj = jE - jS + 1;
  const axP = new Float64Array(nj), ayP = new Float64Array(nj);
  for (let j = 0; j < nj; j++) { axP[j] = TR.fxAll[jS + j] / m; ayP[j] = TR.fyAll[jS + j] / m - G; }
  const mAx = mean(Array.from(axP)), mAy = mean(Array.from(ayP));
  const vxP = cumInt(Float64Array.from(axP, v => v - mAx), 1 / FS_F); const mvx = mean(Array.from(vxP)); for (let j = 0; j < nj; j++) vxP[j] -= mvx;
  const vyP = cumInt(Float64Array.from(ayP, v => v - mAy), 1 / FS_F); const mvy = mean(Array.from(vyP)); for (let j = 0; j < nj; j++) vyP[j] -= mvy;
  const yP = cumInt(vyP, 1 / FS_F);
  r.rangeComYPlate = rangeOf(yP); r.plateMeanAyOffset = mAy;          // mean of Fy/m - g over the stride (should be ~0)
  const plateAt = (arr, t) => lin(arr, Math.min(nj - 1, Math.max(0, t * FS_F - jS)));
  // ---- check 1: vertical
  const ayM = [], ayPl = [], axM = [], axPl = [], flight = [];
  for (let k = k0; k <= k1; k++) {
    const t = k / FS_M, pl = lin(TR.fyAll, t * FS_F) / m - G;
    ayM.push(m0.acy[k]); ayPl.push(pl);
    if (lin(TR.fyAll, t * FS_F) < 20 && Number.isFinite(m0.acy[k])) flight.push(m0.acy[k]);
  }
  r.rmsAy_stride = rmsOf(ayM.map((v, i) => v - ayPl[i])); r.corrAy_stride = corr(ayM, ayPl); r.flightAy = flight.length ? mean(flight) : NaN; r.flightN = flight.length;
  const st_ayM = [], st_ayP = [], st_axM = [], st_axP = [];
  for (let k = kS0; k <= kS1; k++) { const t = k / FS_M; st_ayM.push(m0.acy[k]); st_ayP.push(at(cf.fy, t) / m - G); st_axM.push(m0.acx[k]); st_axP.push(at(cf.fx, t) / m); }
  r.rmsAy_stance = rmsOf(st_ayM.map((v, i) => v - st_ayP[i])); r.corrAy_stance = corr(st_ayM, st_ayP);
  r.rmsAx_stance = rmsOf(st_axM.map((v, i) => v - st_axP[i])); r.corrAx_stance = corr(st_axM, st_axP);
  r.rangeComY = rangeOf(m0.comy.slice(k0, k1 + 1)); r.rangeMidHipY = rangeOf(m0.midy.slice(k0, k1 + 1)); r.rangeStanceHipY = rangeOf(hip.y.slice(k0, k1 + 1));
  const hmc = Float64Array.from(hip.y, (v, k) => v - m0.comy[k]);
  r.rangeHipMinusCom = rangeOf(hmc.slice(k0, k1 + 1));
  // ---- stride curves (21 points over the stride, % of stride)
  const comyMean = mean(fin(Array.from(m0.comy.slice(k0, k1 + 1)))), midMean = mean(fin(Array.from(m0.midy.slice(k0, k1 + 1)))), hipMean = mean(fin(Array.from(hip.y.slice(k0, k1 + 1))));
  PCT.forEach((p, i) => {
    const t = a.t0 + p / 100 * (tN - a.t0);
    r[`s_hipMinusComY_${i}`] = fi(hmc, t); r[`s_comYrel_${i}`] = fi(m0.comy, t) - comyMean; r[`s_midHipYrel_${i}`] = fi(m0.midy, t) - midMean; r[`s_stanceHipYrel_${i}`] = fi(hip.y, t) - hipMean;
    r[`s_aComYModel_${i}`] = fi(m0.acy, t); r[`s_aComYPlate_${i}`] = plateAt(ayP, t) ;
    r[`s_pelvisTilt_${i}`] = fi(TR.tilt, t); r[`s_pelvisRotStance_${i}`] = (leg === 'R' ? 1 : -1) * fi(TR.rot, t);
  });
  const tl = Array.from(TR.tilt.slice(k0, k1 + 1)), rt = Array.from(TR.rot.slice(k0, k1 + 1)).map(v => (leg === 'R' ? 1 : -1) * v);
  r.tiltMean = mean(fin(tl)); r.tiltRange = rangeOf(tl); r.rotStanceMean = mean(fin(rt)); r.rotStanceRange = rangeOf(rt);
  // ---- stance curves, lean independent
  // plate vx relative to strike from the stance force alone, and cumulative integral of Fx/m over the contact
  const fxm = Float64Array.from(cf.fx, v => v / m), cumFx = cumInt(fxm, 1 / FS_F);
  const vModel0 = fi(m0.vcx, a.t0);
  PCT.forEach((p, i) => {
    const t = a.t0 + p / 100 * dur, fy = at(cf.fy, t), fx = at(cf.fx, t), cop = at(cf.cop, t), valid = fy >= MINFY * BW;
    r[`c_aComYModel_${i}`] = fi(m0.acy, t); r[`c_aComYPlate_${i}`] = fy / m - G;
    r[`c_aComXModel_${i}`] = fi(m0.acx, t); r[`c_aComXPlate_${i}`] = fx / m;
    r[`c_vxModelRel_${i}`] = fi(m0.vcx, t) - vModel0; r[`c_vxPlateRel_${i}`] = at(cumFx, t);
    const hvx = fi(hipvx, t);
    r[`c_lagModel_${i}`] = hvx - fi(m0.vcx, t); r[`c_lagPlate_${i}`] = hvx - plateAt(vxP, t);
    r[`c_hipVx_${i}`] = hvx; r[`c_vxComModel_${i}`] = fi(m0.vcx, t); r[`c_vxComPlate_${i}`] = plateAt(vxP, t);
    r[`c_hipY_${i}`] = fi(hip.y, t); r[`c_comY_${i}`] = fi(m0.comy, t);
    r[`c_fxBW_${i}`] = fx / BW; r[`c_fyBW_${i}`] = fy / BW;
    const hx = fi(hip.x, t), hy = fi(hip.y, t);
    r[`c_copAheadHip_${i}`] = valid ? cop - hx : NaN;
    r[`c_fxPredA_${i}`] = valid ? -(cop - hx) * (fy / BW) / hy : NaN;
  });
  r.dvxModel = fi(m0.vcx, a.t1) - vModel0; r.dvxPlate = at(cumFx, a.t1);
  { const v = []; for (let k = kS0; k <= kS1; k++) v.push(m0.vcx[k] - vModel0); r.dipModel = Math.min(...fin(v)); const w = []; for (let k = kS0; k <= kS1; k++) w.push(at(cumFx, k / FS_M)); r.dipPlate = Math.min(...w); }
  // ---- impulses and RMS over valid points, prediction A (hip) and measured
  const frames = []; for (let k = kS0; k <= kS1; k++) { const t = k / FS_M, fy = at(cf.fy, t); if (fy >= MINFY * BW) frames.push({ k, t, fy, fx: at(cf.fx, t), cop: at(cf.cop, t) }); }
  const imp = f => { let b = 0, p = 0; for (const x of frames) { const v = f(x); if (!Number.isFinite(v)) continue; if (v < 0) b -= v; else p += v; } return [b / m / FS_M, p / m / FS_M]; };
  [r.brakeMeasV, r.propMeasV] = imp(x => x.fx);
  [r.brakeA, r.propA] = imp(x => -(x.cop - hip.x[x.k]) * x.fy / hip.y[x.k]);
  { let b = 0, p = 0; for (let j = 0; j < cf.len - 1; j++) { const u = 0.5 * (cf.fx[j] + cf.fx[j + 1]); if (u < 0) b -= u; else p += u; } r.brakeMeas = b / m / FS_F; r.propMeas = p / m / FS_F; }
  { const d = []; for (let i = 2; i <= 18; i++) { const v = r[`c_fxPredA_${i}`]; if (Number.isFinite(v)) d.push(v - r[`c_fxBW_${i}`]); } r.rmsPredA = rmsOf(d); }
  // ---- lean dependent
  variants.forEach((V, vi) => {
    const M = models[vi], sf = '_' + V.name, dyHat = M.dyHat;
    const stIdx = leg === 'R' ? 'R' : 'L', swIdx = oth;
    PCT.forEach((p, i) => {
      const t = a.t0 + p / 100 * dur, fy = at(cf.fy, t), fx = at(cf.fx, t), cop = at(cf.cop, t), valid = fy >= MINFY * BW;
      const cx = fi(M.comx, t), cy = fi(M.comy, t), dH = fi(M.dH, t);
      r[`c_hipMinusComX${sf}_${i}`] = fi(hip.x, t) - cx;
      r[`c_dH${sf}_${i}`] = dH; r[`c_dHstance${sf}_${i}`] = fi(M.dpart[stIdx], t); r[`c_dHswing${sf}_${i}`] = fi(M.dpart[swIdx], t); r[`c_dHhat${sf}_${i}`] = fi(M.dpart.H, t);
      r[`c_offset${sf}_${i}`] = valid ? (cop - cx) + (fx / fy) * cy : NaN;
      r[`c_torque${sf}_${i}`] = valid ? ((cop - cx) * fy + cy * fx) / m : NaN;
      r[`c_fxPredB${sf}_${i}`] = valid ? -(cop - cx) * (fy / BW) / cy : NaN;
      r[`c_fxPredC${sf}_${i}`] = valid ? (dH / G - (cop - cx) * (fy / BW)) / cy : NaN;
    });
    const dB = [], dC = [];
    for (let i = 2; i <= 18; i++) { const b = r[`c_fxPredB${sf}_${i}`], cc = r[`c_fxPredC${sf}_${i}`]; if (Number.isFinite(b)) { dB.push(b - r[`c_fxBW_${i}`]); dC.push(cc - r[`c_fxBW_${i}`]); } }
    r['rmsPredB' + sf] = rmsOf(dB); r['rmsPredC' + sf] = rmsOf(dC);
    [r['brakeB' + sf], r['propB' + sf]] = imp(x => -(x.cop - M.comx[x.k]) * x.fy / M.comy[x.k]);
    [r['brakeC' + sf], r['propC' + sf]] = imp(x => (m * M.dH[x.k] - (x.cop - M.comx[x.k]) * x.fy) / M.comy[x.k]);
    const off = frames.map(x => (x.cop - M.comx[x.k]) + (x.fx / x.fy) * M.comy[x.k]), tq = frames.map(x => ((x.cop - M.comx[x.k]) * x.fy + M.comy[x.k] * x.fx) / m), dh = frames.map(x => M.dH[x.k]);
    r['offsetMean' + sf] = off.length ? mean(off) : NaN; r['torqueMean' + sf] = tq.length ? mean(tq) : NaN; r['dHmeanValid' + sf] = dh.length ? mean(dh) : NaN;
    r['corrTorque' + sf] = corr(tq, dh);
    if (vi === 0) for (let q = 0; q < tq.length; q++) if (Number.isFinite(tq[q] + dh[q])) { acc.n++; acc.sx += dh[q]; acc.sy += tq[q]; acc.sxy += dh[q] * tq[q]; acc.sxx += dh[q] * dh[q]; acc.syy += tq[q] * tq[q]; }
    const avg = arr => { const v = []; for (let k = kS0; k <= kS1; k++) v.push(arr[k]); const f = fin(v); return f.length ? mean(f) : NaN; };
    r['dHmeanStance' + sf] = avg(M.dH); r['dHstanceLegMean' + sf] = avg(M.dpart[stIdx]); r['dHswingLegMean' + sf] = avg(M.dpart[swIdx]); r['dHhatMean' + sf] = avg(M.dpart.H);
    { const v = []; for (let k = k0; k <= k1; k++) v.push(M.dH[k]); r['dHmeanStride' + sf] = mean(fin(v)); }
    { const sd_ = (arr) => PCT.forEach((p, i) => { const t = a.t0 + p / 100 * (tN - a.t0); r[`s_${arr[0]}${sf}_${i}`] = fi(arr[1], t); });
      sd_(['dH', M.dH]); sd_(['dHstance', M.dpart[stIdx]]); sd_(['dHswing', M.dpart[swIdx]]); sd_(['dHhat', M.dpart.H]); }
  });
  return r;
}

// ---------------------------------------------------------------- pass over all subjects
function runPass(variantsFor, info, ix, ids, heights, masses, label) {
  const rows = [], acc = { n: 0, sx: 0, sy: 0, sxy: 0, sxx: 0, syy: 0 }, rej = { timing: 0, nostride: 0, few: 0 };
  for (const sid of ids) {
    const id = String(sid).padStart(3, '0'), k = ids.indexOf(sid), subj = { height: heights[k], mass: masses[k] }, st = loadStatic(id);
    for (const sp of SPEEDS) {
      const TR = prepMomentum(id, sp, st), variants = variantsFor(sp), models = variants.map(V => bodyModel(TR, subj.height, V.deg(subj, sid)));
      const cs = TR.T.contacts, legRecs = {}, fsOf = {};
      const steps = []; for (let i = 0; i + 1 < cs.length; i++) steps.push(cs[i + 1].t0 - cs[i].t0);
      const medStep = median(steps.filter(x => x > 0.2 && x < 1.0));
      for (const s of ['R', 'L']) {
        fsOf[s] = (info[sid][ix(s + 'FSI' + sp)] || '').trim();
        const recs = [];
        for (let i = 0; i + 2 < cs.length; i++) {
          const a = cs[i], b = cs[i + 1], c = cs[i + 2];
          if (a.leg !== s || b.leg === s || c.leg !== s) continue;
          if (Math.abs((b.t0 - a.t0) / medStep - 1) > 0.3 || Math.abs((c.t0 - b.t0) / medStep - 1) > 0.3) { rej.timing++; continue; }
          const r = strideRec(TR, models, variants, a, c, s, subj, acc);
          if (r) recs.push(r); else rej.nostride++;
        }
        if (recs.length < 5) { rej.few++; legRecs[s] = null; continue; }
        const rec = { nStrides: recs.length };
        for (const key in recs[0]) { const v = fin(recs.map(q => q[key])); rec[key] = v.length ? mean(v) : NaN; }
        legRecs[s] = rec;
      }
      rows.push({ id: sid, sp, fsOf, legRecs, height: subj.height });
    }
    process.stderr.write(`\r${label} subject ${sid}/${NSUBJ}`);
  }
  process.stderr.write('\n');
  return { rows, acc, rej };
}

function main() {
  const infoLines = fs.readFileSync(path.join(DIR, 'RBDSinfo.txt'), 'utf8').split(/\r?\n/).filter(Boolean);
  const head = infoLines[0].split('\t'), ix = n => head.indexOf(n), info = {};
  for (const l of infoLines.slice(1)) { const f = l.split('\t'), id = parseInt(f[ix('Subject')], 10); if (id >= 1 && id <= NSUBJ && !info[id]) info[id] = f; }
  const ids = Object.keys(info).map(Number).sort((a, b) => a - b);
  const heights = ids.map(i => parseFloat(info[i][ix('Height')]) / 100), masses = ids.map(i => parseFloat(info[i][ix('Mass')]));
  const men = ids.filter(i => info[i][ix('Gender')].trim() === 'M').length;

  // pass 1: lean 0 only, to find the lean that puts the stance-average force line through the COM
  const p1 = runPass(() => [{ name: 'L0', deg: () => 0 }], info, ix, ids, heights, masses, 'pass 1');
  const bestBy = {};    // speed -> {perSubject: [deg], meanDeg}
  for (const sp of SPEEDS) {
    const per = [];
    for (const r of p1.rows.filter(q => q.sp === sp)) {
      const legs = ['R', 'L'].filter(s => r.legRecs[s]).map(s => r.legRecs[s].offsetMean_L0);
      if (!legs.length) continue;
      per.push(Math.atan(mean(legs) / (F_HAT * HAT_H * r.height)) * DEG);
    }
    bestBy[sp] = { per, meanDeg: mean(per), sdDeg: sd(per), n: per.length };
  }
  const bestDeg = sp => bestBy[sp].meanDeg;
  // pass 2: leans 0, 5, 10 and the per-speed best fit
  const p2 = runPass(sp => [...FIX_LEANS.map(d => ({ name: 'L' + d, deg: () => d })), { name: 'Lbest', deg: () => bestDeg(sp) }], info, ix, ids, heights, masses, 'pass 2');
  const rows = p2.rows;

  // ------------------------------------------------ cohorts
  const avgLegs = legs => { const out = { nLegs: legs.length }; for (const k in legs[0]) { const v = fin(legs.map(l => l[k])); out[k] = v.length ? mean(v) : NaN; } return out; };
  const cohort = { all: { 25: [], 35: [], 45: [] }, rearfoot: { 25: [], 35: [], 45: [] } };
  for (const r of rows) {
    const all = ['R', 'L'].filter(s => r.legRecs[s]).map(s => r.legRecs[s]), rear = ['R', 'L'].filter(s => r.legRecs[s] && r.fsOf[s] === REAR).map(s => r.legRecs[s]);
    if (all.length) cohort.all[r.sp].push({ id: r.id, v: avgLegs(all) });
    if (rear.length) cohort.rearfoot[r.sp].push({ id: r.id, v: avgLegs(rear) });
  }
  const KEYS = Object.keys(cohort.all[25][0].v).filter(k => k !== 'nLegs');
  const scalarKeys = KEYS.filter(k => !/^[cs]_/.test(k)), curveKeys = {};
  for (const k of KEYS) { const mm = /^([cs])_(.*)_(\d+)$/.exec(k); if (mm) (curveKeys[mm[1] + '_' + mm[2]] = curveKeys[mm[1] + '_' + mm[2]] || []).push(+mm[3]); }
  const stat = vals => { const v = fin(vals); return { mean: r4(mean(v)), sd: r4(sd(v)), n: v.length }; };
  const summary = {
    source: { citation: 'Fukuchi RK, Fukuchi CA, Duarte M (2017). A public dataset of running biomechanics and the effects of running speed on lower extremity kinematics and kinetics. PeerJ 5:e3298',
              doi: '10.7717/peerj.3298', figshare: 'https://doi.org/10.6084/m9.figshare.4543435 (article 4543435)', licence: 'CC BY 4.0', note: 'derived summary statistics from the raw marker and force files' },
    generated: new Date().toISOString(),
    method: {
      model: 'sagittal seven-mass model, x forward, y up (plate = 0), metres: thigh, shank, foot of both legs plus one head-arms-trunk point mass; angular momentum H and COM per kg of body mass',
      joints: 'hip = Harrington 2007 pelvis-regression hip joint centre; knee = midpoint of Knee and Knee.Medial (static markers carried by the thigh-cluster rigid fit); ankle = midpoint of Ankle and Ankle.Medial (shank-cluster fit)',
      foot: 'toe tip = Heel.Bottom + 1.5 (MT1/MT5 midpoint - Heel.Bottom); foot centre 0.44 of the way heel to toe tip (0.66 of heel to MT midpoint); foot length heel to toe tip; a single MT marker is used when the other is missing',
      inertia: 'de Leva 1996 male set for all subjects (the one woman too): thigh 0.1416 / 0.4095 / 0.329, shank 0.0433 / 0.4459 / 0.255, foot 0.0137 / 0.44 / 0.257 (mass fraction / COM from proximal / gyration radius, fractions of length); segment lengths per frame in the sagittal plane',
      hat: `point mass ${F_HAT.toFixed(4)} of body mass at the midpoint of the two hip joint centres + (dx, dy), dy = ${HAT_H} x height, dx = dy tan(lean); no rotation; leans 0, 5, 10 deg and a per-speed best fit`,
      bestLean: 'tan(lean) = mean stance offset of the force line from the lean-0 COM / (HAT fraction x dy), per subject and speed; the speed value is the mean over the 28 subjects; offset = (x_cop - x_com) + (Fx/Fy) y_com, i.e. the horizontal distance of the force line from the COM at COM height',
      sign: 'H z component anticlockwise positive; torque of the ground force about the COM = (x_cop - x_com) Fy + y_com Fx = m dH/dt (numerical check in checks.torqueSign)',
      predictions: 'A: Fx = -(x_cop - x_hip) Fy / y_hip (stance hip joint centre); B: same with the COM; C: Fx = (m dH/dt - (x_cop - x_com) Fy) / y_com; compared in body weights over 10-90 % of stance where Fy >= 0.2 BW',
      filtering: 'markers 15 Hz zero-lag 2nd-order Butterworth (as landing-fukuchi.js); velocities and accelerations by central differences at 150 Hz (a = x[k+1] - 2x[k] + x[k-1]); forces 30 Hz zero-lag, per contact for stance curves (COP held where Fy < 0.2 BW before filtering) and over the whole record for stride-wide plate series',
      plate: 'aComY_plate = Fy/m - g; aComX_plate = Fx/m; plate-derived COM vx and y are periodic (mean acceleration removed, mean velocity zero) over the stride; vxPlateRel is the integral of Fx/m from strike',
      stances: 'stance = raw 20 N strike to toe-off; stance curves at 0, 5 ... 100 % (cubic interpolation of 150 Hz series, linear for forces); stride curves at 0 ... 100 % of the stride (strike to next strike of the same leg)',
      pelvis: 'tilt = angle of the mid-ASIS - mid-PSIS line in the sagittal plane, positive = anterior tilt (ASIS below PSIS); rotation = angle of the left-to-right ASIS line seen from above, positive when the right ASIS is ahead; rotStance flips the sign for left-stance strides so that positive = stance-side hip forward',
      cohorts: 'strides -> legs -> subjects as in landing-fukuchi.js; rearfoot = legs whose RFSI/LFSI is "Rearfoot"'
    },
    subjects: { n: ids.length, men, women: ids.length - men },
    lean: {}, cohorts: { all: {}, rearfoot: {} }, checks: {}
  };
  for (const sp of SPEEDS) summary.lean[SPEED_LABEL[sp]] = { bestFitDeg: r4(bestBy[sp].meanDeg), perSubjectSdDeg: r4(bestBy[sp].sdDeg), n: bestBy[sp].n };
  for (const c in cohort) for (const sp of SPEEDS) {
    const rs = cohort[c][sp], e = { n: rs.length, legs: rs.reduce((s, r) => s + r.v.nLegs, 0), scalars: {}, stanceCurves: {}, strideCurves: {} };
    for (const k of scalarKeys) e.scalars[k] = stat(rs.map(r => r.v[k]));
    for (const ck in curveKeys) {
      const o = { mean: [], sd: [], n: [] };
      for (const i of curveKeys[ck].sort((a, b) => a - b)) { const s = stat(rs.map(r => r.v[`${ck.slice(0, 1)}_${ck.slice(2)}_${i}`])); o.mean.push(s.mean); o.sd.push(s.sd); o.n.push(s.n); }
      (ck[0] === 'c' ? e.stanceCurves : e.strideCurves)[ck.slice(2)] = o;
    }
    summary.cohorts[c][SPEED_LABEL[sp]] = e;
  }
  const a = p1.acc, slope = (a.n * a.sxy - a.sx * a.sy) / (a.n * a.sxx - a.sx * a.sx), rr = (a.n * a.sxy - a.sx * a.sy) / Math.sqrt((a.n * a.sxx - a.sx * a.sx) * (a.n * a.syy - a.sy * a.sy));
  summary.checks.torqueSign = { note: 'pooled valid stance frames, lean 0: measured torque of the ground force about the COM = (x_cop - x_com) Fy + y_com Fx (per kg) regressed on model dH/dt; sign is right if the slope is positive and near 1', slope: r4(slope), r: r4(rr), frames: a.n };
  summary.checks.rejected = p2.rej;
  fs.mkdirSync(path.dirname(OUT), { recursive: true });
  fs.writeFileSync(OUT, JSON.stringify(summary));

  // ------------------------------------------------ report
  const out = [], P = s => out.push(s);
  const E = (c, sp) => summary.cohorts[c][SPEED_LABEL[sp]];
  const sc = (c, sp, k, d = 3) => { const s = E(c, sp).scalars[k]; return s && s.mean !== null ? `${s.mean.toFixed(d)} ± ${s.sd === null ? '-' : s.sd.toFixed(d)}` : '-'; };
  P(`subjects ${ids.length}; strides rejected ${JSON.stringify(p2.rej)}; torque check: slope of measured torque on model m dH/dt ${slope.toFixed(2)}, r ${rr.toFixed(2)} (n ${a.n} frames)`);
  P('best-fit HAT lean (deg, mean ± SD over subjects): ' + SPEEDS.map(sp => `${SPEED_LABEL[sp]}: ${bestBy[sp].meanDeg.toFixed(1)} ± ${bestBy[sp].sdDeg.toFixed(1)}`).join(', '));
  const LV = ['L0', 'L5', 'L10', 'Lbest'];
  const HEAD = ['all', 'rearfoot'];
  for (const c of HEAD) {
    P(`\n===== cohort ${c}: scalars (mean ± SD across subjects; n = ${SPEEDS.map(sp => E(c, sp).n).join('/')}) =====`);
    P('metric'.padEnd(30) + SPEEDS.map(sp => (SPEED_LABEL[sp] + ' m/s').padStart(20)).join(''));
    const rowS = (lab, k, d = 3) => P(lab.padEnd(30) + SPEEDS.map(sp => sc(c, sp, k, d).padStart(20)).join(''));
    P('-- check 1: vertical');
    rowS('rmsAy stance (m/s2)', 'rmsAy_stance', 2); rowS('corrAy stance', 'corrAy_stance', 2); rowS('rmsAy stride (m/s2)', 'rmsAy_stride', 2); rowS('corrAy stride', 'corrAy_stride', 2); rowS('model aY in flight (m/s2)', 'flightAy', 2);
    rowS('range COM y (m)', 'rangeComY'); rowS('range COM y from plate (m)', 'rangeComYPlate'); rowS('range mid-hip y (m)', 'rangeMidHipY'); rowS('range stance hip y (m)', 'rangeStanceHipY'); rowS('range hip - COM y (m)', 'rangeHipMinusCom');
    P('-- check 2: fore-aft');
    rowS('rmsAx stance (m/s2)', 'rmsAx_stance', 2); rowS('corrAx stance', 'corrAx_stance', 2); rowS('dvx model strike->TO (m/s)', 'dvxModel'); rowS('dvx plate integral (m/s)', 'dvxPlate'); rowS('min vx dip model (m/s)', 'dipModel'); rowS('min vx dip plate (m/s)', 'dipPlate');
    P('-- check 3: prediction RMS (BW, 10-90 %) and impulses (m/s); measured on the same valid points');
    rowS('measured braking/propulsive', 'brakeMeasV'); rowS('  (propulsive)', 'propMeasV'); rowS('A hip: rms', 'rmsPredA'); rowS('A hip: braking', 'brakeA'); rowS('A hip: propulsive', 'propA');
    for (const v of LV) { rowS(`B COM ${v}: rms`, 'rmsPredB_' + v); rowS(`B COM ${v}: braking`, 'brakeB_' + v); rowS(`B COM ${v}: propulsive`, 'propB_' + v); rowS(`C dH/dt ${v}: rms`, 'rmsPredC_' + v); rowS(`C dH/dt ${v}: braking`, 'brakeC_' + v); rowS(`C dH/dt ${v}: propulsive`, 'propC_' + v); }
    P('-- torque about the COM and dH/dt (N m/kg)');
    for (const v of LV) { rowS(`stance-mean torque ${v}`, 'torqueMean_' + v); rowS(`stance-mean dH/dt (valid) ${v}`, 'dHmeanValid_' + v); rowS(`stride-mean dH/dt ${v}`, 'dHmeanStride_' + v); rowS(`stance-leg part ${v}`, 'dHstanceLegMean_' + v); rowS(`swing-leg part ${v}`, 'dHswingLegMean_' + v); rowS(`HAT part ${v}`, 'dHhatMean_' + v); rowS(`stance-mean force-line offset ${v}`, 'offsetMean_' + v); rowS(`corr torque vs m dH/dt ${v}`, 'corrTorque_' + v, 2); }
    P('-- check 4: pelvis (deg)');
    rowS('tilt mean', 'tiltMean', 1); rowS('tilt range', 'tiltRange', 1); rowS('rotation mean (stance side fwd +)', 'rotStanceMean', 1); rowS('rotation range', 'rotStanceRange', 1);
  }
  const curveTable = (c, sp, list, stride) => {
    const e = E(c, sp), cv = stride ? e.strideCurves : e.stanceCurves;
    P(''.padEnd(34) + PCT.filter((p, i) => i % 2 === 0).map(p => String(p).padStart(12)).join(''));
    for (const [k, lab, d] of list) { const o = cv[k]; if (!o) continue; P(lab.padEnd(34) + o.mean.map((m, i) => i % 2 === 0 ? ((m === null ? '-' : `${m.toFixed(d)}±${o.sd[i] === null ? '-' : o.sd[i].toFixed(d)}`).padStart(12)) : '').join('')); }
  };
  for (const c of HEAD) for (const sp of SPEEDS) {
    P(`\n===== curves, cohort ${c}, ${SPEED_LABEL[sp]} m/s: stance (% of stance, every 10 %) =====`);
    curveTable(c, sp, [['aComYModel', 'aComY model (m/s2)', 1], ['aComYPlate', 'aComY plate Fy/m-g (m/s2)', 1], ['aComXModel', 'aComX model (m/s2)', 1], ['aComXPlate', 'aComX plate Fx/m (m/s2)', 1],
      ['vxModelRel', 'COM vx model - strike (m/s)', 3], ['vxPlateRel', 'int Fx/m dt (m/s)', 3], ['hipVx', 'stance hip vx lab (m/s)', 2], ['vxComModel', 'COM vx model (m/s)', 2], ['vxComPlate', 'COM vx plate periodic (m/s)', 2],
      ['lagModel', 'hip vx - COM vx model (m/s)', 2], ['lagPlate', 'hip vx - COM vx plate (m/s)', 2], ['hipY', 'stance hip y (m)', 3], ['comY', 'COM y (m)', 3], ['fxBW', 'Fx measured (BW)', 3], ['fyBW', 'Fy (BW)', 2], ['copAheadHip', 'cop - hip x (m)', 3], ['fxPredA', 'Fx pred A hip (BW)', 3]], false);
    for (const v of LV) curveTable(c, sp, [['hipMinusComX_' + v, `hip x - COM x ${v} (m)`, 3], ['offset_' + v, `force-line offset ${v} (m)`, 3], ['fxPredB_' + v, `Fx pred B COM ${v} (BW)`, 3], ['fxPredC_' + v, `Fx pred C dH/dt ${v} (BW)`, 3],
      ['torque_' + v, `GRF torque about COM ${v} (Nm/kg)`, 2], ['dH_' + v, `dH/dt ${v} (Nm/kg)`, 2], ['dHstance_' + v, `  stance leg ${v}`, 2], ['dHswing_' + v, `  swing leg ${v}`, 2], ['dHhat_' + v, `  HAT ${v}`, 2]], false);
    P(`--- stride curves, cohort ${c}, ${SPEED_LABEL[sp]} m/s (% of stride)`);
    curveTable(c, sp, [['hipMinusComY', 'stance hip y - COM y (m)', 3], ['comYrel', 'COM y - mean (m)', 3], ['midHipYrel', 'mid-hip y - mean (m)', 3], ['stanceHipYrel', 'stance hip y - mean (m)', 3], ['aComYModel', 'aComY model (m/s2)', 1], ['aComYPlate', 'aComY plate (m/s2)', 1],
      ['pelvisTilt', 'pelvis tilt (deg, + anterior)', 1], ['pelvisRotStance', 'pelvis rotation (deg, + stance hip fwd)', 1], ['dH_Lbest', 'dH/dt Lbest (Nm/kg)', 2], ['dHstance_Lbest', '  stance leg', 2], ['dHswing_Lbest', '  swing leg', 2], ['dHhat_Lbest', '  HAT', 2]], true);
  }
  P('\nlargest subject-level RMS errors of prediction C (lean 0, BW) per speed, "subject: value" (median in brackets) - outlier diagnostic');
  for (const sp of SPEEDS) { const v = cohort.all[sp].map(r => [r.id, r.v.rmsPredC_L0]).sort((x, y) => y[1] - x[1]); P(`  ${SPEED_LABEL[sp]} m/s: ` + v.slice(0, 5).map(q => `${q[0]}: ${q[1].toFixed(2)}`).join(', ') + `  (median ${median(v.map(q => q[1])).toFixed(2)})`); }
  console.log(out.join('\n'));
  console.log('\nwrote ' + path.relative(process.cwd(), OUT));
}

main();
