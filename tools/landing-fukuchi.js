// Measures where the foot lands relative to the hip in the raw Fukuchi 2017 marker and force files (subjects 1-28,
// 2.5 / 3.5 / 4.5 m/s), and reduces it to across-subject statistics: tests/data/fukuchi2017-landing.json.
// Needs data/ from  node tools/fetch-fukuchi.js --raw .   Run: node tools/landing-fukuchi.js [--csv <file>]
//
// Method in short (details are written into the JSON "method" block):
//  * Axes (checked on the data): X is the direction of travel, Y is vertical, Z is mediolateral; mm; markers 150 Hz, forces 300 Hz.
//    The sagittal plane is (X, Y).
//  * The running trials have no hip, knee or ankle markers, so the greater trochanter, lateral femoral epicondyle and lateral
//    malleolus (static-trial markers R.GTR, R.Knee, R.Ankle, ...) are reconstructed in every frame by fitting the static position of
//    the thigh / shank marker cluster to the cluster as measured (rigid-body fit, Horn 1987). Pelvis: ASIS/PSIS/iliac crests.
//  * Foot strike / toe-off: vertical force (Fy) crossing 20 N; the leg is the one whose lowest foot marker is lower.
//  * Variant "hjc": hip = Harrington et al. (2007) pelvis-regression hip joint centre instead of the trochanter marker.
'use strict';
const fs = require('fs');
const path = require('path');

// --force-lp <Hz>: optional zero-lag low-pass of the force before thresholding (sensitivity check only); --out / --csv <file>: output paths
const flagVal = name => { const i = process.argv.indexOf(name); return i > 0 ? process.argv[i + 1] : null; };
const DIR = path.join(__dirname, '..', 'data', 'fukuchi2017');
const OUT = flagVal('--out') || path.join(__dirname, '..', 'tests', 'data', 'fukuchi2017-landing.json');
const CSV = flagVal('--csv') || path.join(DIR, 'landing_subjects.csv');
const SPEEDS = ['25', '35', '45'], SPEED_LABEL = { 25: '2.5', 35: '3.5', 45: '4.5' }, SPEED_MS = { 25: 2.5, 35: 3.5, 45: 4.5 };
const NSUBJ = 28, REAR = 'Rearfoot';
const FS_M = 150, FS_F = 300, THR = 20;     // Hz, Hz, N (threshold used by the paper)
const FC = 15, MAXGAP = 20, MINRUN = 25;    // marker low-pass (Hz), longest interpolated gap and shortest usable run (frames)
const DEG = 180 / Math.PI;
const TAU = [-30, -20, -10, 0, 10, 20, 30];   // ms around a force event for the force / marker alignment check
const FORCE_LP = flagVal('--force-lp') ? parseFloat(flagVal('--force-lp')) : 0;

const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
const sd = a => { if (a.length < 2) return NaN; const m = mean(a); return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / (a.length - 1)); };
const median = a => { const s = [...a].sort((x, y) => x - y), n = s.length; return n ? (n % 2 ? s[n >> 1] : 0.5 * (s[n / 2 - 1] + s[n / 2])) : NaN; };
const r4 = x => Number.isFinite(x) ? Number(x.toPrecision(4)) : null;
const fin = a => a.filter(Number.isFinite);

// ---------------------------------------------------------------- file reading
function readTable(file) {
  const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(l => l.trim() !== '');
  const head = lines[0].split('\t').map(s => s.trim());
  const n = lines.length - 1, cols = {};
  head.forEach(h => { cols[h] = new Float64Array(n); });
  for (let i = 0; i < n; i++) {
    const f = lines[i + 1].split('\t');
    for (let j = 0; j < head.length; j++) cols[head[j]][i] = parseFloat(f[j]);
  }
  return { head, cols, n };
}
// marker "base" -> {x,y,z} (all three columns must exist)
function markers(tab) {
  const out = {};
  for (const h of tab.head) {
    if (!/X$/.test(h)) continue;
    const b = h.slice(0, -1);
    if (tab.cols[b + 'Y'] && tab.cols[b + 'Z']) out[b] = { x: tab.cols[b + 'X'], y: tab.cols[b + 'Y'], z: tab.cols[b + 'Z'] };
  }
  return out;
}
function staticMeans(M) {
  const out = {};
  for (const b in M) {
    const m = ['x', 'y', 'z'].map(c => fin(Array.from(M[b][c])));
    if (m.every(a => a.length > 10)) out[b] = m.map(mean);
  }
  return out;
}

// ---------------------------------------------------------------- rigid-body fit (Horn 1987, quaternion form)
function topEigvec4(A) {
  const V = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
  for (let sweep = 0; sweep < 30; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) off += A[p * 4 + q] ** 2;
    if (off < 1e-20) break;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) {
      const apq = A[p * 4 + q];
      if (Math.abs(apq) < 1e-30) continue;
      const th = (A[q * 4 + q] - A[p * 4 + p]) / (2 * apq);
      const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 4; k++) { const a = A[k * 4 + p], b = A[k * 4 + q]; A[k * 4 + p] = c * a - s * b; A[k * 4 + q] = s * a + c * b; }
      for (let k = 0; k < 4; k++) { const a = A[p * 4 + k], b = A[q * 4 + k]; A[p * 4 + k] = c * a - s * b; A[q * 4 + k] = s * a + c * b; }
      for (let k = 0; k < 4; k++) { const a = V[k * 4 + p], b = V[k * 4 + q]; V[k * 4 + p] = c * a - s * b; V[k * 4 + q] = s * a + c * b; }
    }
  }
  let im = 0;
  for (let i = 1; i < 4; i++) if (A[i * 5] > A[im * 5]) im = i;
  return [V[im], V[4 + im], V[8 + im], V[12 + im]];
}
// refPts: {name: [x,y,z]} static positions; trackNames: cluster markers; virt: names of static points to carry into each frame.
// Returns { out: {name: {x,y,z}}, rms } (NaN where fewer than 3 cluster markers are visible; rms = mean fit residual, mm)
function rigidFit(M, refPts, trackNames, virt, n) {
  const tn = trackNames.filter(nm => refPts[nm] && M[nm]);
  const out = {};
  for (const v of virt) out[v] = { x: new Float64Array(n).fill(NaN), y: new Float64Array(n).fill(NaN), z: new Float64Array(n).fill(NaN) };
  let rs = 0, rn = 0;
  const idx = [];
  for (let k = 0; k < n; k++) {
    idx.length = 0;
    for (let m = 0; m < tn.length; m++) { const q = M[tn[m]]; if (Number.isFinite(q.x[k]) && Number.isFinite(q.y[k]) && Number.isFinite(q.z[k])) idx.push(m); }
    if (idx.length < 3) continue;
    const pc = [0, 0, 0], qc = [0, 0, 0];
    for (const m of idx) { const p = refPts[tn[m]], q = M[tn[m]]; for (let a = 0; a < 3; a++) pc[a] += p[a]; qc[0] += q.x[k]; qc[1] += q.y[k]; qc[2] += q.z[k]; }
    for (let a = 0; a < 3; a++) { pc[a] /= idx.length; qc[a] /= idx.length; }
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const m of idx) {
      const p = refPts[tn[m]], q = M[tn[m]], pp = [p[0] - pc[0], p[1] - pc[1], p[2] - pc[2]], qq = [q.x[k] - qc[0], q.y[k] - qc[1], q.z[k] - qc[2]];
      for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) S[a][b] += pp[a] * qq[b];
    }
    const [Sxx, Sxy, Sxz] = S[0], [Syx, Syy, Syz] = S[1], [Szx, Szy, Szz] = S[2];
    const q4 = topEigvec4([
      Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx,
      Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz,
      Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy,
      Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz]);
    const [w, x, y, z] = q4, nq = Math.hypot(w, x, y, z);
    const W = w / nq, X = x / nq, Y = y / nq, Z = z / nq;
    const R = [[W * W + X * X - Y * Y - Z * Z, 2 * (X * Y - W * Z), 2 * (X * Z + W * Y)],
               [2 * (X * Y + W * Z), W * W - X * X + Y * Y - Z * Z, 2 * (Y * Z - W * X)],
               [2 * (X * Z - W * Y), 2 * (Y * Z + W * X), W * W - X * X - Y * Y + Z * Z]];
    const apply = p => { const d = [p[0] - pc[0], p[1] - pc[1], p[2] - pc[2]]; return [0, 1, 2].map(a => R[a][0] * d[0] + R[a][1] * d[1] + R[a][2] * d[2] + qc[a]); };
    for (const v of virt) { const o = apply(refPts[v]); out[v].x[k] = o[0]; out[v].y[k] = o[1]; out[v].z[k] = o[2]; }
    let e2 = 0;
    for (const m of idx) { const o = apply(refPts[tn[m]]), q = M[tn[m]]; e2 += (o[0] - q.x[k]) ** 2 + (o[1] - q.y[k]) ** 2 + (o[2] - q.z[k]) ** 2; }
    rs += Math.sqrt(e2 / idx.length); rn++;
  }
  return { out, rms: rn ? rs / rn : NaN };
}

// ---------------------------------------------------------------- signal helpers
// fill gaps up to MAXGAP frames by linear interpolation, drop runs shorter than MINRUN, zero-lag 2nd-order Butterworth
// (forward and backward) with odd-reflection padding on each remaining run
function conditioned(src) {
  const a = Float64Array.from(src), n = a.length;
  for (let i = 0; i < n;) {
    if (Number.isFinite(a[i])) { i++; continue; }
    let j = i; while (j < n && !Number.isFinite(a[j])) j++;
    if (i > 0 && j < n && j - i <= MAXGAP) for (let k = i; k < j; k++) a[k] = a[i - 1] + (a[j] - a[i - 1]) * (k - i + 1) / (j - i + 1);
    i = j;
  }
  const K = Math.tan(Math.PI * FC / FS_M), nrm = 1 / (1 + Math.SQRT2 * K + K * K);
  const b0 = K * K * nrm, b1 = 2 * b0, b2 = b0, a1 = 2 * (K * K - 1) * nrm, a2 = (1 - Math.SQRT2 * K + K * K) * nrm;
  const pass = x => { let z1 = 0, z2 = 0; const y = new Float64Array(x.length);
    z1 = x[0] * (1 - b0) ; z2 = x[0] * (b2 - a2) + 0; // steady-state start for a constant input
    for (let i = 0; i < x.length; i++) { const o = b0 * x[i] + z1; z1 = b1 * x[i] - a1 * o + z2; z2 = b2 * x[i] - a2 * o; y[i] = o; }
    return y; };
  const out = new Float64Array(n).fill(NaN);
  for (let i = 0; i < n;) {
    if (!Number.isFinite(a[i])) { i++; continue; }
    let j = i; while (j < n && Number.isFinite(a[j])) j++;
    const len = j - i;
    if (len >= MINRUN) {
      const pad = Math.min(15, len - 1), x = new Float64Array(len + 2 * pad);
      for (let k = 0; k < pad; k++) { x[pad - 1 - k] = 2 * a[i] - a[i + 1 + k]; x[pad + len + k] = 2 * a[j - 1] - a[j - 2 - k]; }
      for (let k = 0; k < len; k++) x[pad + k] = a[i + k];
      const f = pass(x).reverse(), g = pass(f).reverse();
      for (let k = 0; k < len; k++) out[i + k] = g[pad + k];
    }
    i = j;
  }
  return out;
}
const cond2 = p => ({ x: conditioned(p.x), y: conditioned(p.y) });

function cubic(a, u) {      // Catmull-Rom interpolation at fractional index u
  const i = Math.floor(u), f = u - i;
  if (i < 1 || i + 2 >= a.length) return NaN;
  const p0 = a[i - 1], p1 = a[i], p2 = a[i + 1], p3 = a[i + 2];
  return p1 + 0.5 * f * (p2 - p0 + f * (2 * p0 - 5 * p1 + 4 * p2 - p3 + f * (3 * (p1 - p2) + p3 - p0)));
}
const ptAt = (s, t) => ({ x: cubic(s.x, t * FS_M), y: cubic(s.y, t * FS_M) });

// ---------------------------------------------------------------- pelvis hip joint centre (Harrington et al. 2007)
// pelvis frame: origin mid-ASIS, Z = left ASIS -> right ASIS, X anterior (mid-PSIS -> mid-ASIS, made orthogonal to Z), Y = Z x X.
//   x = -0.24 PD - 9.9,  y = -0.30 PW - 10.9,  z = +-(0.33 PW + 7.3)   (mm; PW = inter-ASIS distance, PD = mid-PSIS to mid-ASIS distance)
function hjcFrom(RA, LA, RP, LP, side) {
  const mA = [0, 1, 2].map(i => (RA[i] + LA[i]) / 2), mP = [0, 1, 2].map(i => (RP[i] + LP[i]) / 2);
  const zv = [RA[0] - LA[0], RA[1] - LA[1], RA[2] - LA[2]], PW = Math.hypot(...zv), Z = zv.map(v => v / PW);
  const A = [mA[0] - mP[0], mA[1] - mP[1], mA[2] - mP[2]], dot = A[0] * Z[0] + A[1] * Z[1] + A[2] * Z[2];
  const xv = [A[0] - dot * Z[0], A[1] - dot * Z[1], A[2] - dot * Z[2]], xn = Math.hypot(...xv), X = xv.map(v => v / xn);
  const Y = [Z[1] * X[2] - Z[2] * X[1], Z[2] * X[0] - Z[0] * X[2], Z[0] * X[1] - Z[1] * X[0]];
  const PD = Math.hypot(...A), cx = -0.24 * PD - 9.9, cy = -0.30 * PW - 10.9, cz = (side === 'R' ? 1 : -1) * (0.33 * PW + 7.3);
  return [0, 1, 2].map(i => mA[i] + X[i] * cx + Y[i] * cy + Z[i] * cz);
}

// ---------------------------------------------------------------- static trial
const THIGH = s => [`${s}.Thigh.Top.Lateral`, `${s}.Thigh.Bottom.Lateral`, `${s}.Thigh.Top.Medial`, `${s}.Thigh.Bottom.Medial`];
const SHANK = s => [`${s}.Shank.Top.Lateral`, `${s}.Shank.Bottom.Lateral`, `${s}.Shank.Top.Medial`, `${s}.Shank.Bottom.Medial`];
const PELVIS = ['R.ASIS', 'L.ASIS', 'R.PSIS', 'L.PSIS', 'R.Iliac.Crest', 'L.Iliac.Crest'];

const ang = (dx, dz) => Math.atan2(dx, dz) * DEG;   // angle from vertical of a vector (dx forward, dz upward-pointing)
// sagittal snapshot of one leg from points {hip, knee, ankle, heel, fore} (each {x, y}, mm) -> degrees / metres
function snapshot(P, ground, foreOff) {
  const ankleAhead = (P.ankle.x - P.hip.x) / 1000, hipHeight = (P.hip.y - ground) / 1000;
  const thigh = ang(P.knee.x - P.hip.x, P.hip.y - P.knee.y), shank = ang(P.ankle.x - P.knee.x, P.knee.y - P.ankle.y);
  return {
    ankleAhead, heelAhead: (P.heel.x - P.hip.x) / 1000, hipHeight,
    landDeg: Math.atan(ankleAhead / hipHeight) * DEG,
    thighDeg: thigh, shankDeg: shank, kneeDeg: thigh - shank,      // thigh - shank = 180 deg minus the hip-knee-ankle included angle (signed: flexion +)
    footDeg: P.fore ? Math.atan2(P.fore.y - P.heel.y, P.fore.x - P.heel.x) * DEG - foreOff[P.foreMode] : NaN
  };
}

function loadStatic(id) {
  const tab = readTable(path.join(DIR, `RBDS${id}static.txt`)), M = markers(tab), ref = staticMeans(M);
  const res = { ref, legs: {} };
  for (const s of ['R', 'L']) {
    res.legs[s] = {};
    const heel = ref[`${s}.Heel.Bottom`], mt1 = ref[`${s}.MT1`], mt5 = ref[`${s}.MT5`];
    const foreOff = { both: NaN, mt1: NaN, mt5: NaN };
    const fa = (h, f) => Math.atan2(f[1] - h[1], f[0] - h[0]) * DEG;
    if (heel && mt1 && mt5) foreOff.both = fa(heel, [(mt1[0] + mt5[0]) / 2, (mt1[1] + mt5[1]) / 2]);
    if (heel && mt1) foreOff.mt1 = fa(heel, mt1);
    if (heel && mt5) foreOff.mt5 = fa(heel, mt5);
    res.legs[s].foreOff = foreOff;
    res.legs[s].footLowY = Math.min(...[heel, mt1, mt5].filter(Boolean).map(p => p[1]));
  }
  // pelvis hip joint centres from the static pelvis markers
  const need = ['R.ASIS', 'L.ASIS', 'R.PSIS', 'L.PSIS'];
  if (need.every(n => ref[n])) for (const s of ['R', 'L']) res.legs[s].hjc = hjcFrom(ref['R.ASIS'], ref['L.ASIS'], ref['R.PSIS'], ref['L.PSIS'], s);
  res.ground = mean([res.legs.R.footLowY, res.legs.L.footLowY]);
  return res;
}

// ---------------------------------------------------------------- one running trial: reconstruct, condition, find contacts
function lowpassForce(x, fc) {      // zero-lag 2nd-order Butterworth, forward and backward (sensitivity option only)
  const K = Math.tan(Math.PI * fc / FS_F), nrm = 1 / (1 + Math.SQRT2 * K + K * K);
  const b0 = K * K * nrm, b1 = 2 * b0, a1 = 2 * (K * K - 1) * nrm, a2 = (1 - Math.SQRT2 * K + K * K) * nrm;
  const pass = v => { let z1 = 0, z2 = 0; const y = new Float64Array(v.length); for (let i = 0; i < v.length; i++) { const o = b0 * v[i] + z1; z1 = b1 * v[i] - a1 * o + z2; z2 = b0 * v[i] - a2 * o; y[i] = o; } return y; };
  return pass(pass(Float64Array.from(x)).reverse()).reverse();
}
function detectContacts(fy0) {
  const fy = FORCE_LP ? lowpassForce(fy0, FORCE_LP) : fy0;
  const n = fy.length, s = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const w = [];
    for (let d = -2; d <= 2; d++) w.push(fy[Math.min(n - 1, Math.max(0, i + d))]);
    w.sort((a, b) => a - b); s[i] = w[2];
  }
  const out = [];
  for (let i = 0; i < n;) {
    if (!(s[i] > THR)) { i++; continue; }
    let j = i, pk = s[i];
    while (j + 1 < n && s[j + 1] > THR) { j++; if (s[j] > pk) pk = s[j]; }
    if (i > 0 && j < n - 1) {
      const t0 = (i - 1 + (THR - s[i - 1]) / (s[i] - s[i - 1])) / FS_F, t1 = (j + (s[j] - THR) / (s[j] - s[j + 1])) / FS_F;
      out.push({ t0, t1, dur: t1 - t0, peak: pk });
    }
    i = j + 1;
  }
  return out;
}

function prepTrial(id, sp, st, log) {
  const mt = readTable(path.join(DIR, `RBDS${id}runT${sp}markers.txt`)), ft = readTable(path.join(DIR, `RBDS${id}runT${sp}forces.txt`));
  const M = markers(mt), n = mt.n, T = { id, sp, n, legs: {}, notes: [] };
  const pel = rigidFit(M, st.ref, PELVIS, ['R.ASIS', 'L.ASIS', 'R.PSIS', 'L.PSIS'], n);
  log.rms.pelvis.push(pel.rms);
  const hjc = {};
  for (const s of ['R', 'L']) {
    hjc[s] = { x: new Float64Array(n).fill(NaN), y: new Float64Array(n).fill(NaN) };
    for (let k = 0; k < n; k++) {
      const P = nm => [pel.out[nm].x[k], pel.out[nm].y[k], pel.out[nm].z[k]];
      const RA = P('R.ASIS');
      if (!Number.isFinite(RA[0])) continue;
      const h = hjcFrom(RA, P('L.ASIS'), P('R.PSIS'), P('L.PSIS'), s);
      hjc[s].x[k] = h[0]; hjc[s].y[k] = h[1];
    }
  }
  for (const s of ['R', 'L']) {
    const th = rigidFit(M, st.ref, THIGH(s), [`${s}.GTR`, `${s}.Knee`], n), sk = rigidFit(M, st.ref, SHANK(s), [`${s}.Ankle`], n);
    log.rms.thigh.push(th.rms); log.rms.shank.push(sk.rms);
    const raw = nm => M[nm] ? { x: M[nm].x, y: M[nm].y } : { x: new Float64Array(n).fill(NaN), y: new Float64Array(n).fill(NaN) };
    const L = {
      gtr: cond2(th.out[`${s}.GTR`]), hjc: cond2(hjc[s]), knee: cond2(th.out[`${s}.Knee`]), ankle: cond2(sk.out[`${s}.Ankle`]),
      heel: cond2(raw(`${s}.Heel.Bottom`)), mt1: cond2(raw(`${s}.MT1`)), mt5: cond2(raw(`${s}.MT5`))
    };
    L.lowY = new Float64Array(n);       // lowest foot marker height per frame
    for (let k = 0; k < n; k++) { const v = fin([L.heel.y[k], L.mt1.y[k], L.mt5.y[k]]); L.lowY[k] = v.length ? Math.min(...v) : NaN; }
    // per-frame angle series for the two hip variants
    for (const v of ['gtr', 'hjc']) {
      const thigh = new Float64Array(n), knee = new Float64Array(n);
      for (let k = 0; k < n; k++) {
        const hx = L[v].x[k], hy = L[v].y[k];
        thigh[k] = ang(L.knee.x[k] - hx, hy - L.knee.y[k]);
        knee[k] = thigh[k] - ang(L.ankle.x[k] - L.knee.x[k], L.knee.y[k] - L.ankle.y[k]);
      }
      L[v].thigh = thigh; L[v].kneeSeries = knee;
    }
    T.legs[s] = L;
  }
  // contacts from the vertical force (Fy), then clean ones: duration within 0.6-1.5 x the median of the real contacts
  const all = detectContacts(ft.cols.Fy), real = all.filter(c => c.dur >= 0.08), med = median(real.map(c => c.dur));
  T.nRaw = all.length;
  T.contacts = real.filter(c => c.dur >= 0.6 * med && c.dur <= 1.5 * med);
  T.nDropped = all.length - T.contacts.length;
  // leg of each contact: the foot whose lowest marker is lower over the middle 40 % of the contact
  const ambiguous = [];
  for (const c of T.contacts) {
    const k0 = Math.round((c.t0 + 0.3 * c.dur) * FS_M), k1 = Math.round((c.t0 + 0.7 * c.dur) * FS_M), m = {};
    for (const s of ['R', 'L']) m[s] = mean(fin(Array.from(T.legs[s].lowY.slice(k0, k1 + 1))));
    if (!Number.isFinite(m.R) || !Number.isFinite(m.L) || Math.abs(m.R - m.L) < 20) { c.leg = null; ambiguous.push(c); continue; }
    c.leg = m.R < m.L ? 'R' : 'L'; c.low = m[c.leg];
    // COP cross-check: the stance foot is the one whose mean x (heel + forefoot) is closer to the centre of pressure in x
    const kc = Math.round((c.t0 + 0.5 * c.dur) * FS_M), ci = Math.round((c.t0 + 0.5 * c.dur) * FS_F), cop = ft.cols.COPx[ci], dx = {};
    for (const s of ['R', 'L']) { const l = T.legs[s]; dx[s] = Math.abs(mean(fin([l.heel.x[kc], l.mt1.x[kc], l.mt5.x[kc]])) - cop); }
    c.copAgrees = Number.isFinite(dx.R + dx.L) ? ((dx.R < dx.L) === (c.leg === 'R')) : null;
  }
  T.contacts = T.contacts.filter(c => c.leg);
  T.nAmbiguous = ambiguous.length;
  return T;
}

// ---------------------------------------------------------------- strides of a trial for one leg
const HIPDEP = ['ankleAhead', 'heelAhead', 'hipHeight', 'landDeg', 'thighDeg', 'kneeDeg'];
function strideMetrics(T, st, subj, ground, leg, rej) {
  const L = T.legs[leg], cs = T.contacts, recs = [];
  const stepTimes = []; for (let i = 0; i + 1 < cs.length; i++) stepTimes.push(cs[i + 1].t0 - cs[i].t0);
  const medStep = median(stepTimes.filter(x => x > 0.2 && x < 1.0));
  const foreOff = st.legs[leg].foreOff;
  for (let i = 0; i + 2 < cs.length; i++) {
    const a = cs[i], b = cs[i + 1], c = cs[i + 2];
    if (a.leg !== leg || b.leg === leg || c.leg !== leg) continue;
    const s1 = b.t0 - a.t0, s2 = c.t0 - b.t0;
    if (Math.abs(s1 / medStep - 1) > 0.3 || Math.abs(s2 / medStep - 1) > 0.3) { rej.timing++; continue; }
    const tFS = a.t0, tTO = a.t1, tN = c.t0, Tstride = tN - tFS;
    const k0 = Math.ceil(tFS * FS_M), k1 = Math.floor(tN * FS_M), kt = Math.floor(tTO * FS_M);
    const r = { strideTime: Tstride, stepRate: 120 / Tstride, contactMs: a.dur * 1000, stancePct: 100 * a.dur / Tstride };
    // belt speed: slope of x of the lowest foot marker over the middle 20 % of the contact
    {
      const w0 = Math.round((tFS + 0.4 * a.dur) * FS_M), w1 = Math.round((tFS + 0.6 * a.dur) * FS_M);
      let best = null;
      for (const nm of ['heel', 'mt1', 'mt5']) {
        const xs = [], ys = [];
        for (let k = w0; k <= w1; k++) if (Number.isFinite(L[nm].x[k]) && Number.isFinite(L[nm].y[k])) { xs.push(k / FS_M); ys.push(L[nm].y[k]); }
        if (xs.length >= 4) { const my = mean(ys); if (!best || my < best.my) best = { nm, my }; }
      }
      if (best) {
        const ts = [], xs = [];
        for (let k = w0; k <= w1; k++) if (Number.isFinite(L[best.nm].x[k])) { ts.push(k / FS_M); xs.push(L[best.nm].x[k]); }
        const mt = mean(ts), mx = mean(xs); let num = 0, den = 0;
        for (let q = 0; q < ts.length; q++) { num += (ts[q] - mt) * (xs[q] - mx); den += (ts[q] - mt) ** 2; }
        r.beltSpeed = -num / den / 1000;
      }
    }
    let ok = true;
    for (const v of ['gtr', 'hjc']) {
      const pre = v === 'gtr' ? '' : 'hjc_';
      const ev = (t, tag) => {
        const P = { hip: ptAt(L[v], t), knee: ptAt(L.knee, t), ankle: ptAt(L.ankle, t), heel: ptAt(L.heel, t) };
        const m1 = ptAt(L.mt1, t), m5 = ptAt(L.mt5, t);
        if (Number.isFinite(m1.x) && Number.isFinite(m5.x)) { P.fore = { x: (m1.x + m5.x) / 2, y: (m1.y + m5.y) / 2 }; P.foreMode = 'both'; }
        else if (Number.isFinite(m5.x)) { P.fore = m5; P.foreMode = 'mt5'; }
        else if (Number.isFinite(m1.x)) { P.fore = m1; P.foreMode = 'mt1'; }
        const sn = snapshot(P, ground, foreOff);
        if (![P.hip.x, P.hip.y, P.knee.x, P.knee.y, P.ankle.x, P.ankle.y, P.heel.x, P.heel.y].every(Number.isFinite)) return null;
        const o = {};
        for (const k in sn) if (!(v === 'hjc' && !HIPDEP.includes(k))) o[`${pre}${tag}_${k}`] = sn[k];
        o[`${pre}${tag}_ankleAheadNorm`] = sn.ankleAhead / subj.height;
        return { o, hipY: P.hip.y };
      };
      const fsv = ev(tFS, 'fs'), tov = ev(tTO, 'to');
      if (!fsv || !tov) { ok = false; break; }
      Object.assign(r, fsv.o, tov.o);
      // extremes over the swing (toe-off -> next strike) and stance (strike -> toe-off), from the 150 Hz series
      const rng = (a0, a1, arr) => { let mx = -Infinity, mn = Infinity; for (let k = a0; k <= a1; k++) { const x = arr[k]; if (!Number.isFinite(x)) return null; if (x > mx) mx = x; if (x < mn) mn = x; } return { mx, mn }; };
      const sw = rng(kt + 1, k1, L[v].thigh), swK = rng(kt + 1, k1, L[v].kneeSeries), stK = rng(k0, kt, L[v].kneeSeries);
      const hipY = L[v].y, stH = rng(k0, kt, hipY), strH = rng(k0, k1, hipY), strT = rng(k0, k1, L[v].thigh);
      if (!sw || !swK || !stK || !stH || !strH || !strT) { ok = false; break; }
      r[pre + 'peakThighFwd'] = sw.mx; r[pre + 'peakThighBack'] = sw.mn; r[pre + 'peakKneeSwing'] = swK.mx;
      r[pre + 'thighMinStride'] = strT.mn; r[pre + 'peakKneeStance'] = stK.mx;
      r[pre + 'minHipHeight'] = (Math.min(stH.mn, fsv.hipY) - ground) / 1000;
      r[pre + 'hipDrop'] = (fsv.hipY - Math.min(stH.mn, fsv.hipY)) / 1000;
      r[pre + 'hipRange'] = (strH.mx - strH.mn) / 1000;
    }
    if (!ok) { rej.nan++; continue; }
    recs.push(r);
  }
  return recs;
}

// ---------------------------------------------------------------- processed-file comparison (same legs and strides logic as reduce-fukuchi)
function processedKnee(ptab, side, sp) {
  const c = nm => ptab.cols[side + nm + sp];
  const knee = c('kneeAngZ'), gy = c('grfY');
  if (!knee || !gy || knee.length !== 101 || Array.from(knee).some(x => !Number.isFinite(x) || Math.abs(x) > 1000)) return null;
  let kTO = -1; while (kTO + 1 < 101 && Number.isFinite(gy[kTO + 1])) kTO++;
  if (kTO < 5 || kTO > 80) return null;
  const kn = Array.from(knee);
  return { proc_kneeFS: kn[0], proc_kneePeakStance: Math.max(...kn.slice(0, kTO + 1)), proc_kneePeakSwing: Math.max(...kn.slice(kTO + 1)) };
}

// ---------------------------------------------------------------- regression
function ols(x, y) {
  const p = []; for (let i = 0; i < x.length; i++) if (Number.isFinite(x[i]) && Number.isFinite(y[i])) p.push([x[i], y[i]]);
  const n = p.length; if (n < 4) return null;
  const mx = mean(p.map(q => q[0])), my = mean(p.map(q => q[1]));
  let sxx = 0, sxy = 0, syy = 0;
  for (const [a, b] of p) { sxx += (a - mx) ** 2; sxy += (a - mx) * (b - my); syy += (b - my) ** 2; }
  const slope = sxy / sxx, rss = syy - slope * sxy, se = Math.sqrt(rss / (n - 2) / sxx);
  return { slope, se, r: sxy / Math.sqrt(sxx * syy), n };
}
const fmtReg = g => g ? `${g.slope.toPrecision(3)} ± ${g.se.toPrecision(2)}  r=${g.r.toFixed(2)}  n=${g.n}` : '-';

// ---------------------------------------------------------------- main
function main() {
  const infoLines = fs.readFileSync(path.join(DIR, 'RBDSinfo.txt'), 'utf8').split(/\r?\n/).filter(Boolean);
  const head = infoLines[0].split('\t'), ix = n => head.indexOf(n);
  const info = {};
  for (const l of infoLines.slice(1)) { const f = l.split('\t'), id = parseInt(f[ix('Subject')], 10); if (id >= 1 && id <= NSUBJ && !info[id]) info[id] = f; }
  const ids = Object.keys(info).map(Number).sort((a, b) => a - b);
  const heights = ids.map(i => parseFloat(info[i][ix('Height')]) / 100), masses = ids.map(i => parseFloat(info[i][ix('Mass')]));
  const men = ids.filter(i => info[i][ix('Gender')].trim() === 'M').length;

  const log = { rms: { pelvis: [], thigh: [], shank: [] }, rej: { timing: 0, nan: 0, fewStrides: 0 }, trials: [], prof: { 25: { fs: TAU.map(() => []), to: TAU.map(() => []) }, 35: { fs: TAU.map(() => []), to: TAU.map(() => []) }, 45: { fs: TAU.map(() => []), to: TAU.map(() => []) } },
                contactsKept: 0, contactsRaw: 0, ambiguous: 0, copAgree: 0, copTotal: 0, alternation: { ok: 0, total: 0 }, hipSign: { neg: 0, tot: 0 },
                staticOffsets: [], groundMm: [], strides: { 25: 0, 35: 0, 45: 0 }, hjcMinusGtr: [] };
  const subjRows = [];                                   // one per subject and speed
  const legVals = { R: { 25: [], 35: [], 45: [] }, L: { 25: [], 35: [], 45: [] } };   // for the left/right comparison: {id, rec}
  const fsNames = new Set();

  for (const sid of ids) {
    const id = String(sid).padStart(3, '0'), height = heights[ids.indexOf(sid)];
    const st = loadStatic(id), ptab = readTable(path.join(DIR, `RBDS${id}processed.txt`)), subj = { height };
    const trials = SPEEDS.map(sp => prepTrial(id, sp, st, log));
    // ground: median lowest foot-marker height at mid-stance over every clean contact of the subject
    const gs = []; for (const T of trials) for (const c of T.contacts) gs.push(c.low);
    const ground = median(gs); log.groundMm.push(ground);
    // standing values (static trial, same definitions; mean marker positions)
    const stand = {};
    {
      const legs = [];
      for (const s of ['R', 'L']) {
        const r = st.ref, g = r[`${s}.GTR`], k = r[`${s}.Knee`], a = r[`${s}.Ankle`];
        if (!g || !k || !a) continue;
        const o = { gtr: snapshot({ hip: { x: g[0], y: g[1] }, knee: { x: k[0], y: k[1] }, ankle: { x: a[0], y: a[1] }, heel: { x: 0, y: 0 } }, st.ground, {}) };
        const h = st.legs[s].hjc;
        if (h) { o.hjc = snapshot({ hip: { x: h[0], y: h[1] }, knee: { x: k[0], y: k[1] }, ankle: { x: a[0], y: a[1] }, heel: { x: 0, y: 0 } }, st.ground, {}); log.hjcMinusGtr.push([h[0] - g[0], h[1] - g[1]]); }
        o.leg = Math.hypot(g[0] - k[0], g[1] - k[1], g[2] - k[2]) + Math.hypot(k[0] - a[0], k[1] - a[1], k[2] - a[2]);
        legs.push(o);
      }
      stand.stand_kneeDeg = mean(legs.map(o => o.gtr.kneeDeg)); stand.stand_thighDeg = mean(legs.map(o => o.gtr.thighDeg));
      stand.stand_hipHeight = mean(legs.map(o => o.gtr.hipHeight)); stand.stand_hipHeightOverHeight = stand.stand_hipHeight / height;
      stand.stand_legLength = mean(legs.map(o => o.leg)) / 1000;
      const hj = legs.filter(o => o.hjc);
      if (hj.length) { stand.hjc_stand_kneeDeg = mean(hj.map(o => o.hjc.kneeDeg)); stand.hjc_stand_thighDeg = mean(hj.map(o => o.hjc.thighDeg)); stand.hjc_stand_hipHeight = mean(hj.map(o => o.hjc.hipHeight)); }
      log.staticOffsets.push(stand);
    }
    for (let si = 0; si < 3; si++) {
      const sp = SPEEDS[si], T = trials[si];
      log.contactsRaw += T.nRaw; log.contactsKept += T.contacts.length; log.ambiguous += T.nAmbiguous;
      for (const c of T.contacts) if (c.copAgrees !== null) { log.copTotal++; if (c.copAgrees) log.copAgree++; }
      for (let i = 0; i + 1 < T.contacts.length; i++) { log.alternation.total++; if (T.contacts[i].leg !== T.contacts[i + 1].leg) log.alternation.ok++; }
      // force / marker alignment: height of the lowest foot marker above the ground level around the force strike and toe-off, rearfoot legs
      for (const c of T.contacts) {
        if ((info[sid][ix(c.leg + 'FSI' + sp)] || '').trim() !== REAR) continue;
        const lowY = T.legs[c.leg].lowY;
        TAU.forEach((tau, j) => {
          const a = cubic(lowY, (c.t0 + tau / 1000) * FS_M) - ground, b = cubic(lowY, (c.t1 + tau / 1000) * FS_M) - ground;
          if (Number.isFinite(a)) { log.prof[sp].fs[j].push(a); }
          if (Number.isFinite(b)) { log.prof[sp].to[j].push(b); }
        });
      }
      const rej = log.rej;
      const legRecs = {}, fsOf = {};
      for (const s of ['R', 'L']) {
        fsOf[s] = (info[sid][ix(s + 'FSI' + sp)] || '').trim(); fsNames.add(fsOf[s]);
        const recs = strideMetrics(T, st, subj, ground, s, rej);
        if (recs.length < 5) { if (recs.length) rej.fewStrides++; legRecs[s] = null; continue; }
        const rec = { nStrides: recs.length };
        for (const k in recs[0]) { const v = fin(recs.map(r => r[k])); rec[k] = v.length ? mean(v) : NaN; }
        Object.assign(rec, stand);
        Object.assign(rec, processedKnee(ptab, s, sp) || { proc_kneeFS: NaN, proc_kneePeakStance: NaN, proc_kneePeakSwing: NaN });
        rec.diff_kneeFS = rec.fs_kneeDeg - rec.proc_kneeFS; rec.diff_kneePeakStance = rec.peakKneeStance - rec.proc_kneePeakStance; rec.diff_kneePeakSwing = rec.peakKneeSwing - rec.proc_kneePeakSwing;
        rec.contactOverHeight = rec.contactMs / 1000 * SPEED_MS[sp] / height;
        rec.stridesPerMin = 60 / rec.strideTime;
        legRecs[s] = rec; log.strides[sp] += recs.length;
        legVals[s][sp].push({ id: sid, rec });
      }
      subjRows.push({ id: sid, sp, height, mass: masses[ids.indexOf(sid)], fsOf, legRecs, ground });
    }
    process.stdout.write(`\rsubject ${sid}/${NSUBJ}`);
  }
  process.stdout.write('\n');

  // ------------------------------------------------ cohorts: legs averaged within subject, then across subjects
  const avgLegs = legs => {
    const out = { nLegs: legs.length };
    for (const k in legs[0]) { const v = fin(legs.map(l => l[k])); out[k] = v.length ? mean(v) : NaN; }
    return out;
  };
  const cohortRows = { all: { 25: [], 35: [], 45: [] }, rearfoot: { 25: [], 35: [], 45: [] } };
  for (const r of subjRows) {
    const all = ['R', 'L'].filter(s => r.legRecs[s]).map(s => r.legRecs[s]);
    const rear = ['R', 'L'].filter(s => r.legRecs[s] && r.fsOf[s] === REAR).map(s => r.legRecs[s]);
    if (all.length) cohortRows.all[r.sp].push({ id: r.id, height: r.height, v: avgLegs(all) });
    if (rear.length) cohortRows.rearfoot[r.sp].push({ id: r.id, height: r.height, v: avgLegs(rear) });
  }
  const KEYS = Object.keys(cohortRows.all[25][0].v).filter(k => k !== 'nLegs');
  const stat = vals => { const v = fin(vals); return { mean: r4(mean(v)), sd: r4(sd(v)), n: v.length }; };
  const summary = {
    source: { citation: 'Fukuchi RK, Fukuchi CA, Duarte M (2017). A public dataset of running biomechanics and the effects of running speed on lower extremity kinematics and kinetics. PeerJ 5:e3298',
              doi: '10.7717/peerj.3298', figshare: 'https://doi.org/10.6084/m9.figshare.4543435 (article 4543435)', licence: 'CC BY 4.0', note: 'derived summary statistics from the raw marker and force files' },
    generated: new Date().toISOString(),
    method: {
      axes: 'X = direction of travel (forward +), Y = vertical, Z = mediolateral; mm; markers 150 Hz, forces 300 Hz (Fy vertical, single plate under the treadmill). Sagittal plane = (X, Y).',
      events: `foot strike / toe-off = first / last sample of Fy (5-sample median${FORCE_LP ? `, plus ${FORCE_LP} Hz zero-lag low-pass` : ''}) above 20 N, linearly interpolated; contacts lasting 0.6-1.5 x the median contact; stride = strike to next strike of the same leg with alternating legs and step times within 30 % of the median`,
      leg: 'contact assigned to the foot whose lowest marker (heel bottom, MT1, MT5) is lower over the middle 40 % of the contact',
      hip: 'gtr = greater trochanter reconstructed from the thigh marker cluster (static R.GTR / L.GTR carried by a rigid-body fit); hjc = Harrington 2007 pelvis regression on ASIS/PSIS (pelvis marker group fitted rigidly)',
      knee: 'lateral femoral epicondyle (static R.Knee) carried by the thigh-cluster fit', ankle: 'lateral malleolus (static R.Ankle) carried by the shank-cluster fit',
      heel: 'R.Heel.Bottom marker', forefoot: 'midpoint of MT1 and MT5 (a single one if the other is missing; the static angle is taken with the same choice)',
      filter: `short gaps (<= ${MAXGAP} frames) interpolated; zero-lag 2nd-order Butterworth forward and backward, nominal ${FC} Hz; events by cubic interpolation`,
      ground: 'per subject: median over all clean contacts of the lowest foot-marker height at mid-stance (marker height, not the belt surface)',
      kneeDeg: 'thighDeg - shankDeg (equals 180 deg minus the hip-knee-ankle included angle when the knee is flexed)',
      footDeg: 'heel-bottom to forefoot line above horizontal minus the same angle in the static trial (toe up positive)',
      swing: 'toe-off to next strike; stance: strike to toe-off', cohorts: 'legs averaged within subject; rearfoot = legs whose RFSI/LFSI value is "Rearfoot" (as tools/reduce-fukuchi.js)'
    },
    subjects: { n: ids.length, heightM: r4(mean(heights)), massKg: r4(mean(masses)), men, women: ids.length - men },
    cohorts: { all: {}, rearfoot: {} }
  };
  for (const c in cohortRows) for (const sp of SPEEDS) {
    const rows = cohortRows[c][sp], e = { n: rows.length, legs: rows.reduce((s, r) => s + r.v.nLegs, 0), strides: Math.round(rows.reduce((s, r) => s + r.v.nStrides * r.v.nLegs, 0)), scalars: {} };
    for (const k of KEYS) e.scalars[k] = stat(rows.map(r => r.v[k]));
    summary.cohorts[c][SPEED_LABEL[sp]] = e;
  }

  // ------------------------------------------------ regressions
  const regs = { vsStepRate: {}, vsSpeed: {} };
  const RATE_Y = [['fs_ankleAheadNorm', 'ankleAhead/height'], ['fs_landDeg', 'landDeg'], ['peakThighFwd', 'peak fwd thighDeg'], ['peakKneeSwing', 'peak swing kneeDeg'], ['contactMs', 'contactMs'], ['fs_thighDeg', 'thighDeg at strike']];
  const SPEED_Y = [['fs_landDeg', 'landDeg'], ['fs_thighDeg', 'thighDeg at strike'], ['fs_kneeDeg', 'kneeDeg at strike'], ['fs_ankleAheadNorm', 'ankleAhead/height'], ['peakThighFwd', 'peak fwd thighDeg'],
                    ['peakKneeSwing', 'peak swing kneeDeg'], ['contactOverHeight', 'contactMs*speed/height (-)'], ['stepRate', 'step rate (/min)']];
  const reg = (c, pre, kx, ky) => {
    const out = {};
    for (const sp of SPEEDS) { const rows = cohortRows[c][sp]; out[SPEED_LABEL[sp]] = ols(rows.map(r => r.v[kx]), rows.map(r => r.v[ky])); }
    return out;
  };
  for (const c of ['all', 'rearfoot']) {
    regs.vsStepRate[c] = {}; regs.vsSpeed[c] = {};
    for (const [k] of RATE_Y) regs.vsStepRate[c][k] = reg(c, '', 'stepRate', k);
    for (const [k] of SPEED_Y) {
      const xs = [], ys = [], bySubj = {};
      for (const sp of SPEEDS) for (const r of cohortRows[c][sp]) { xs.push(SPEED_MS[sp]); ys.push(r.v[k]); (bySubj[r.id] = bySubj[r.id] || []).push([SPEED_MS[sp], r.v[k]]); }
      const pooled = ols(xs, ys);
      // within-subject slope from the 3 speeds (n=3, so computed directly)
      const w = []; for (const s in bySubj) { const p = bySubj[s].filter(q => Number.isFinite(q[1])); if (p.length === 3) { const mx = mean(p.map(q => q[0])), my = mean(p.map(q => q[1])); w.push(p.reduce((a, q) => a + (q[0] - mx) * (q[1] - my), 0) / p.reduce((a, q) => a + (q[0] - mx) ** 2, 0)); } }
      regs.vsSpeed[c][k] = { pooled, within: w.length > 2 ? { slope: mean(w), se: sd(w) / Math.sqrt(w.length), n: w.length } : null };
    }
  }
  const roundReg = g => g ? { slope: r4(g.slope), se: r4(g.se), r: r4(g.r), n: g.n } : null;
  summary.regressions = { vsStepRate: {}, vsSpeedPooled: {}, vsSpeedWithinSubject: {} };
  for (const c in regs.vsStepRate) {
    summary.regressions.vsStepRate[c] = {}; summary.regressions.vsSpeedPooled[c] = {}; summary.regressions.vsSpeedWithinSubject[c] = {};
    for (const k in regs.vsStepRate[c]) { summary.regressions.vsStepRate[c][k] = {}; for (const sp in regs.vsStepRate[c][k]) summary.regressions.vsStepRate[c][k][sp] = roundReg(regs.vsStepRate[c][k][sp]); }
    for (const k in regs.vsSpeed[c]) {
      summary.regressions.vsSpeedPooled[c][k] = roundReg(regs.vsSpeed[c][k].pooled);
      const w = regs.vsSpeed[c][k].within; summary.regressions.vsSpeedWithinSubject[c][k] = w ? { slope: r4(w.slope), se: r4(w.se), n: w.n } : null;
    }
  }
  fs.mkdirSync(path.dirname(OUT), { recursive: true });
  fs.writeFileSync(OUT, JSON.stringify(summary));

  // ------------------------------------------------ CSV
  {
    const baseCols = ['subject', 'speed', 'height_m', 'mass_kg', 'footstrike', 'groundMarkerY_mm', 'legs', 'strides'];
    const lines = [baseCols.concat(KEYS).join(',')];
    for (const r of subjRows) {
      const legs = ['R', 'L'].filter(s => r.legRecs[s]).map(s => r.legRecs[s]);
      if (!legs.length) continue;
      const a = avgLegs(legs), num = x => Number.isFinite(x) ? Number(x.toPrecision(5)) : '';
      lines.push([r.id, SPEED_LABEL[r.sp], r.height, r.mass, `${r.fsOf.R || '?'}/${r.fsOf.L || '?'}`, num(r.ground), ['R', 'L'].filter(s => r.legRecs[s]).join(''), legs.reduce((s, l) => s + l.nStrides, 0)].concat(KEYS.map(k => num(a[k]))).join(','));
    }
    fs.mkdirSync(path.dirname(CSV), { recursive: true });
    fs.writeFileSync(CSV, lines.join('\n') + '\n');
  }

  // ------------------------------------------------ report
  const rep = [];
  const P = s => rep.push(s);
  const ms = (c, sp, k, d = 2) => { const s = summary.cohorts[c][SPEED_LABEL[sp]].scalars[k]; return s && s.mean !== null ? `${s.mean.toFixed(d)} ± ${s.sd === null ? '-' : s.sd.toFixed(d)}` : '-'; };
  P(`subjects ${ids.length} (${men} men, ${ids.length - men} women); foot-strike values seen: ${[...fsNames].join(', ')}`);
  P(`contacts: ${log.contactsRaw} force contacts above threshold, ${log.contactsKept} kept as clean and leg-assigned (${log.ambiguous} ambiguous leg); alternating-leg consecutive pairs ${log.alternation.ok}/${log.alternation.total}; ` +
    `COP cross-check agrees with the lowest-foot assignment in ${log.copAgree}/${log.copTotal}`);
  P(`strides used: ${SPEEDS.map(sp => SPEED_LABEL[sp] + ': ' + log.strides[sp]).join(', ')}; rejected strides ${JSON.stringify(log.rej)}; leg-trials with < 5 clean strides are excluded`);
  P(`cluster-fit residual (mean over trials of the per-frame RMS, mm): pelvis ${mean(fin(log.rms.pelvis)).toFixed(1)}, thigh ${mean(fin(log.rms.thigh)).toFixed(1)}, shank ${mean(fin(log.rms.shank)).toFixed(1)}`);
  P(`force / marker alignment: mean height (mm) of the lowest foot marker above ground level, rearfoot legs, relative to the force event; offsets in ms: ${TAU.join(' ')}`);
  for (const sp of SPEEDS) for (const ev of ['fs', 'to']) P(`  ${SPEED_LABEL[sp]} m/s ${ev === 'fs' ? 'strike  ' : 'toe-off '} ` + log.prof[sp][ev].map(v => mean(v).toFixed(1).padStart(6)).join(' '));
  P(`ground marker height (mm above the plate, subject medians): ${mean(log.groundMm).toFixed(1)} ± ${sd(log.groundMm).toFixed(1)}`);
  if (log.hjcMinusGtr.length) P(`standing sagittal offset HJC minus trochanter marker (mm): x ${mean(log.hjcMinusGtr.map(p => p[0])).toFixed(1)}, y ${mean(log.hjcMinusGtr.map(p => p[1])).toFixed(1)}`);
  const so = k => { const v = fin(log.staticOffsets.map(s => s[k])); return `${mean(v).toFixed(2)} ± ${sd(v).toFixed(2)}`; };
  P(`standing (static): kneeDeg ${so('stand_kneeDeg')}, thighDeg ${so('stand_thighDeg')}, hip height / height ${so('stand_hipHeightOverHeight')}, leg length proxy ${so('stand_legLength')} m (hjc variant kneeDeg ${so('hjc_stand_kneeDeg')}, thighDeg ${so('hjc_stand_thighDeg')})`);

  const TABLE = [['strideTime', 's', 3], ['stridesPerMin', '/min', 1], ['stepRate', 'steps/min', 1], ['contactMs', 'ms', 0], ['stancePct', '% stride', 1], ['beltSpeed', 'm/s', 2],
    ['fs_ankleAhead', 'm', 3], ['fs_ankleAheadNorm', '/height', 3], ['fs_heelAhead', 'm', 3], ['fs_hipHeight', 'm', 3], ['fs_landDeg', 'deg', 1], ['fs_thighDeg', 'deg', 1], ['fs_kneeDeg', 'deg', 1], ['fs_shankDeg', 'deg', 1], ['fs_footDeg', 'deg', 1],
    ['to_ankleAhead', 'm', 3], ['to_ankleAheadNorm', '/height', 3], ['to_heelAhead', 'm', 3], ['to_hipHeight', 'm', 3], ['to_landDeg', 'deg', 1], ['to_thighDeg', 'deg', 1], ['to_kneeDeg', 'deg', 1], ['to_shankDeg', 'deg', 1], ['to_footDeg', 'deg', 1],
    ['peakThighFwd', 'deg', 1], ['peakThighBack', 'deg', 1], ['peakKneeSwing', 'deg', 1], ['peakKneeStance', 'deg', 1], ['minHipHeight', 'm', 3], ['hipDrop', 'm', 3], ['hipRange', 'm', 3],
    ['hjc_fs_ankleAhead', 'm', 3], ['hjc_fs_ankleAheadNorm', '/height', 3], ['hjc_fs_hipHeight', 'm', 3], ['hjc_fs_landDeg', 'deg', 1], ['hjc_fs_thighDeg', 'deg', 1], ['hjc_fs_kneeDeg', 'deg', 1],
    ['hjc_to_landDeg', 'deg', 1], ['hjc_to_thighDeg', 'deg', 1], ['hjc_peakThighFwd', 'deg', 1], ['hjc_peakThighBack', 'deg', 1], ['hjc_peakKneeSwing', 'deg', 1], ['hjc_peakKneeStance', 'deg', 1], ['hjc_hipDrop', 'm', 3], ['hjc_hipRange', 'm', 3],
    ['proc_kneeFS', 'deg', 1], ['diff_kneeFS', 'deg', 1], ['proc_kneePeakStance', 'deg', 1], ['diff_kneePeakStance', 'deg', 1], ['proc_kneePeakSwing', 'deg', 1], ['diff_kneePeakSwing', 'deg', 1],
    ['stand_kneeDeg', 'deg', 1], ['stand_thighDeg', 'deg', 1], ['stand_hipHeightOverHeight', '-', 3], ['stand_legLength', 'm', 3]];
  for (const c of ['all', 'rearfoot']) {
    P(`\ncohort ${c}   (subjects / legs per speed: ` + SPEEDS.map(sp => `${SPEED_LABEL[sp]}: ${summary.cohorts[c][SPEED_LABEL[sp]].n}/${summary.cohorts[c][SPEED_LABEL[sp]].legs}`).join(', ') + ')');
    P('metric'.padEnd(26) + SPEEDS.map(sp => (SPEED_LABEL[sp] + ' m/s').padStart(20)).join(''));
    for (const [k, u, d] of TABLE) P(`${k} [${u}]`.padEnd(26) + SPEEDS.map(sp => ms(c, sp, k, d).padStart(20)).join(''));
  }
  // left / right
  P('\nleft vs right leg (all cohort, subject-level; mean R, mean L, mean paired R-L ± SE, n):');
  for (const k of ['strideTime', 'contactMs', 'fs_landDeg', 'fs_ankleAheadNorm', 'fs_thighDeg', 'fs_kneeDeg', 'peakThighFwd', 'peakKneeSwing', 'hipDrop', 'stand_kneeDeg']) {
    P('  ' + k.padEnd(18) + SPEEDS.map(sp => {
      const R = new Map(legVals.R[sp].map(o => [o.id, o.rec[k]])), Lm = new Map(legVals.L[sp].map(o => [o.id, o.rec[k]])), d = [], a = [], b = [];
      for (const [id, v] of R) if (Lm.has(id) && Number.isFinite(v) && Number.isFinite(Lm.get(id))) { d.push(v - Lm.get(id)); a.push(v); b.push(Lm.get(id)); }
      const dd = Math.abs(mean(a)) > 20 ? 1 : 3;
      return `${SPEED_LABEL[sp]}: ${mean(a).toFixed(dd)} / ${mean(b).toFixed(dd)} / ${mean(d).toFixed(dd)} ± ${(sd(d) / Math.sqrt(d.length)).toFixed(dd)} (${d.length})`;
    }).join('   '));
  }
  // regressions
  P('\nregressions across subjects, slope ± SE per (step/min)  r  n');
  for (const c of ['all', 'rearfoot']) {
    P(`\n  cohort ${c}: y against step rate`);
    for (const [k, lab] of RATE_Y) P('    ' + lab.padEnd(22) + SPEEDS.map(sp => `${SPEED_LABEL[sp]}: ${fmtReg(regs.vsStepRate[c][k][SPEED_LABEL[sp]])}`).join('  |  '));
  }
  P('\npooled over speeds, slope ± SE per (m/s)  r  n   [pooled treats the 3 speeds of a subject as independent]; within-subject = mean ± SE of the per-subject 3-point slopes');
  for (const c of ['all', 'rearfoot']) {
    P(`\n  cohort ${c}: y against speed`);
    for (const [k, lab] of SPEED_Y) { const g = regs.vsSpeed[c][k]; P('    ' + lab.padEnd(28) + 'pooled ' + fmtReg(g.pooled).padEnd(40) + 'within-subject ' + (g.within ? `${g.within.slope.toPrecision(3)} ± ${g.within.se.toPrecision(2)} (n=${g.within.n})` : '-')); }
  }
  console.log(rep.join('\n'));
  console.log('\nwrote ' + path.relative(process.cwd(), OUT) + ' and ' + CSV);
}

main();
