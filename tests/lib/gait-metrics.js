// Gait metrics shared by the dataset reducer and the model validator, so both are measured with the same definitions.
//
// Input: a stride record in dataset conventions
//   { T,                         stride time, s
//     n,                         samples on a uniform stride-fraction grid from foot strike (both ends included)
//     kTO,                       index of toe-off on that grid (stance = 0..kTO, swing = kTO+1..n-1)
//     grfX, grfY,                N/kg, 0 outside stance
//     hipAng, kneeAng, ankleAng, deg (ankleAng may be NaN outside stance)
//     hipMom, kneeMom, ankleMom, N m/kg: hip flexor +, knee extensor +, ankle plantarflexor +
//     hipPow, kneePow, anklePow }  W/kg
// Output: { scalars, stanceCurves (51 points over 0-100% of stance), strideCurves (n points over the stride) }
'use strict';

const G = 9.81;
const NSTANCE = 51;
const STANCE_CURVES = ['grfX', 'grfY', 'kneeAng', 'ankleAng', 'hipMom', 'kneeMom', 'ankleMom', 'hipPow', 'kneePow', 'anklePow'];
const STRIDE_CURVES = ['kneeAng', 'hipMom', 'kneeMom', 'hipPow', 'kneePow'];

// linear interpolation at a fractional index
function at(a, x) {
  const n = a.length;
  if (x <= 0) return a[0];
  if (x >= n - 1) return a[n - 1];
  const i = Math.floor(x), f = x - i;
  return f === 0 ? a[i] : a[i] + (a[i + 1] - a[i]) * f;
}
// max / argmax over indices lo..hi (inclusive); NaN samples are ignored
function maxOver(a, lo, hi) {
  let m = -Infinity, im = -1;
  for (let k = lo; k <= hi; k++) if (a[k] > m) { m = a[k]; im = k; }
  return { v: im < 0 ? NaN : m, i: im };
}
function minOver(a, lo, hi) {
  let m = Infinity, im = -1;
  for (let k = lo; k <= hi; k++) if (a[k] < m) { m = a[k]; im = k; }
  return { v: im < 0 ? NaN : m, i: im };
}
// trapezoid work of a power curve between samples lo and hi, positive and negative parts apart (J/kg per stride)
function work(p, lo, hi, T, n) {
  const h = T / (n - 1);
  let pos = 0, neg = 0;
  for (let k = lo; k < hi; k++) {
    const a = p[k], b = p[k + 1];
    pos += 0.5 * (Math.max(a, 0) + Math.max(b, 0)) * h;
    neg += 0.5 * (Math.min(a, 0) + Math.min(b, 0)) * h;
  }
  return [pos, neg];
}
function area(a, T, n) {   // trapezoid over the whole stride fraction, times T
  let s = 0;
  for (let k = 0; k < n - 1; k++) s += 0.5 * (a[k] + a[k + 1]);
  return s * T / (n - 1);
}

function computeMetrics(r) {
  const { T, n, kTO } = r, last = n - 1, s = {};
  s.contactMs = kTO / last * T * 1000;
  s.stancePct = 100 * kTO / last;

  const pk = maxOver(r.grfY, 0, kTO);
  s.peakFz = pk.v / G;
  // time of the peak: parabola through the three samples around it (the grid is coarse), as % of stance
  let pkI = pk.i;
  if (pkI > 0 && pkI < kTO) {
    const a = r.grfY[pkI - 1], b = r.grfY[pkI], c = r.grfY[pkI + 1], den = a - 2 * b + c;
    if (den < 0) pkI += 0.5 * (a - c) / den;
  }
  s.peakFzAt = 100 * pkI / kTO;

  const mn = minOver(r.grfX, 0, kTO), mx = maxOver(r.grfX, 0, kTO);
  s.peakBrake = mn.v / G;
  s.peakProp = mx.v / G;
  s.brakeImpulse = area(r.grfX.map(x => Math.min(x, 0)), T, n);
  s.propImpulse = area(r.grfX.map(x => Math.max(x, 0)), T, n);
  s.brakeEndsAt = NaN;
  for (let k = mn.i; k < kTO; k++) {
    if (r.grfX[k] < 0 && r.grfX[k + 1] >= 0) { s.brakeEndsAt = 100 * (k + (-r.grfX[k]) / (r.grfX[k + 1] - r.grfX[k])) / kTO; break; }
  }
  s.vertImpulseCheck = (area(r.grfY, 1, n)) / (G / 2);

  s.kneeAtStrike = r.kneeAng[0];
  s.kneePeakStance = maxOver(r.kneeAng, 0, kTO).v;
  s.kneeAtToeOff = r.kneeAng[kTO];
  s.kneePeakSwing = maxOver(r.kneeAng, kTO + 1, last).v;
  s.ankleAtStrike = r.ankleAng[0];
  s.anklePeakDF = maxOver(r.ankleAng, 0, kTO).v;
  s.ankleAtToeOff = r.ankleAng[kTO];

  // ankle as a spring: ankle moment against ankle angle over stance samples lo..hi (no interpolation; at least 4 samples, else NaN)
  const fitMomAng = (lo, hi) => {
    const m = hi - lo + 1;
    if (m < 4) return null;
    let sx = 0, sy = 0;
    for (let k = lo; k <= hi; k++) { if (!Number.isFinite(r.ankleAng[k]) || !Number.isFinite(r.ankleMom[k])) return null; sx += r.ankleAng[k]; sy += r.ankleMom[k]; }
    const mx = sx / m, my = sy / m;
    let sxx = 0, sxy = 0, syy = 0;
    for (let k = lo; k <= hi; k++) { const dx = r.ankleAng[k] - mx, dy = r.ankleMom[k] - my; sxx += dx * dx; sxy += dx * dy; syy += dy * dy; }
    if (!(sxx > 0)) return null;
    const slope = sxy / sxx;
    return { slope, mx, my, r2: syy > 0 ? sxy * sxy / (sxx * syy) : NaN };
  };
  const angPk = maxOver(r.ankleAng, 0, kTO), momPk = maxOver(r.ankleMom, 0, kTO);
  const load = angPk.i < 0 ? null : fitMomAng(Math.ceil(0.2 * kTO - 1e-9), angPk.i);
  const unload = momPk.i < 0 ? null : fitMomAng(momPk.i, Math.floor(0.9 * kTO + 1e-9));
  s.ankleStiffLoad = load ? load.slope : NaN;
  s.ankleEngageDeg = load && load.slope !== 0 ? load.mx - load.my / load.slope : NaN;   // angle where the fitted line crosses zero moment
  s.ankleStiffR2 = load ? load.r2 : NaN;
  s.ankleStiffUnload = unload ? unload.slope : NaN;
  s.anklePeakDFAt = angPk.i < 0 ? NaN : 100 * angPk.i / kTO;

  const hmax = maxOver(r.hipAng, 0, last).v, hmin = minOver(r.hipAng, 0, last).v;
  s.hipROM = hmax - hmin;

  s.kneeExtMom = maxOver(r.kneeMom, 0, kTO).v;
  s.anklePFMom = maxOver(r.ankleMom, 0, kTO).v;
  s.hipExtMomStance = -minOver(r.hipMom, 0, kTO).v;
  s.hipFlxMomStance = maxOver(r.hipMom, 0, kTO).v;
  s.hipFlxMomSwing = maxOver(r.hipMom, kTO + 1, last).v;
  s.hipExtMomSwing = -minOver(r.hipMom, kTO + 1, last).v;
  s.kneeFlxMomSwing = -minOver(r.kneeMom, kTO + 1, last).v;
  // extra, whole-stride extremes (used to compare with published peak hip moments)
  s.hipExtMomAll = -minOver(r.hipMom, 0, last).v;
  s.hipFlxMomAll = maxOver(r.hipMom, 0, last).v;

  let sumPos = 0;
  for (const j of ['hip', 'knee', 'ankle']) {
    const p = r[j + 'Pow'];
    const [ps, ns] = work(p, 0, kTO, T, n), [pw, nw] = work(p, kTO, last, T, n);
    s[j + 'PosStance'] = ps; s[j + 'NegStance'] = ns; s[j + 'PosSwing'] = pw; s[j + 'NegSwing'] = nw;
    sumPos += ps + pw;
  }
  s.sumPos = sumPos;

  for (let pct = 10; pct <= 90; pct += 10) {
    const x = pct / 100 * kTO, fx = at(r.grfX, x), fy = at(r.grfY, x), mh = at(r.hipMom, x);
    s['fxOverFz_' + pct] = fx / fy;
    s['hipLever_' + pct] = -mh / Math.sqrt(fx * fx + fy * fy);
  }

  const stanceCurves = {}, strideCurves = {};
  for (const c of STANCE_CURVES) {
    const out = new Array(NSTANCE);
    for (let i = 0; i < NSTANCE; i++) out[i] = at(r[c], i / (NSTANCE - 1) * kTO);
    stanceCurves[c] = out;
  }
  for (const c of STRIDE_CURVES) strideCurves[c] = Array.from(r[c]);
  return { scalars: s, stanceCurves, strideCurves };
}

module.exports = { computeMetrics, G, NSTANCE, STANCE_CURVES, STRIDE_CURVES };
