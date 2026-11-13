// Compares the model with measured running (Fukuchi et al. 2017, summarised in tests/data/fukuchi2017-summary.json).
// A report, not a gate: always exits 0 unless it crashes.
//   node tests/validate.js [--cohort all|rearfoot] [--speed 3.5] [--set key=value ...] [--json out.json]
'use strict';
const fs = require('fs');
const path = require('path');
const { computeMetrics, G, NSTANCE, STANCE_CURVES } = require('./lib/gait-metrics.js');

// load the engine as a function body (same pattern as check.js)
const { simulate } = new Function(fs.readFileSync(path.join(__dirname, '..', 'src', 'physics.js'), 'utf8') + '\nreturn { simulate };')();

const NGRID = 101;

// periodic linear interpolation at fractional index x of arr (length N; arr[N] wraps to arr[0]); NaN neighbours are ignored
function sample(arr, x) {
  const N = arr.length, i = Math.floor(x), f = x - i, a = arr[i % N], b = arr[(i + 1) % N];
  if (f === 0) return a;
  if (Number.isFinite(a) && Number.isFinite(b)) return a + (b - a) * f;
  return Number.isFinite(a) ? a : Number.isFinite(b) ? b : NaN;
}

// model result -> stride record in dataset conventions
function modelRecord(s, M) {
  const { N, T, K, S } = s, D = 180 / Math.PI;
  const i0 = K.stance.findIndex(v => v === 1);
  if (i0 < 0) throw new Error('model has no stance samples');
  const rot = f => { const out = new Float64Array(N); for (let j = 0; j < N; j++) out[j] = f((i0 + j) % N); return out; };
  const st = i => K.stance[i] === 1;
  const fields = {
    grfX: rot(i => S.Fx[i] * G), grfY: rot(i => S.Fy[i] * G),
    hipAng: rot(i => K.th[i] * D), kneeAng: rot(i => K.ph[i] * D),
    ankleAng: rot(i => st(i) ? (-(K.th[i] - K.ph[i]) + K.footAng[i]) * D : NaN),
    hipMom: rot(i => S.tauH[i] / M), kneeMom: rot(i => S.tauK[i] / M), ankleMom: rot(i => -S.tauA[i] / M),
    hipPow: rot(i => S.Ph[i] / M), kneePow: rot(i => S.Pk[i] / M), anklePow: rot(i => S.Pa[i] / M)
  };
  const kTO = Math.round(s.info.tc / T * (NGRID - 1));
  const rec = { T, n: NGRID, kTO };
  for (const k in fields) {
    const out = new Array(NGRID);
    for (let g = 0; g < NGRID; g++) out[g] = sample(fields[k], g / (NGRID - 1) * N);
    rec[k] = out;
  }
  // forces are zero outside stance; the ankle has no moment or power in swing (foot mass is neglected); ankle angle is stance-only
  // a grid point at or just past the last stance sample (kTO is rounded) takes that last stance sample
  const nS = N - i0;
  for (let g = 0; g <= kTO; g++) for (const k of ['ankleAng', 'ankleMom', 'anklePow']) if (!Number.isFinite(rec[k][g])) rec[k][g] = fields[k][nS - 1];
  for (let g = 0; g < NGRID; g++) {
    if (g > kTO) { rec.grfX[g] = 0; rec.grfY[g] = 0; rec.ankleAng[g] = NaN; }
    if (!Number.isFinite(rec.ankleMom[g])) rec.ankleMom[g] = 0;
    if (!Number.isFinite(rec.anklePow[g])) rec.anklePow[g] = 0;
  }
  return rec;
}

// Geometry rows: model quantities from the simulation result s (stance block K.stance === 1) with the measured reference mean/SD.
// fl = landing summary forceLine block (curves at 0,5..100 % of stance), ls = landing summary scalars, ro = rollover summary (scalars + stanceCurves).
function geometryRows(s, H, fl, ls, ro) {
  const K = s.K, D = 180 / Math.PI;
  const i0 = K.stance.findIndex(v => v === 1);
  let nS = 0;
  while (i0 + nS < K.stance.length && K.stance[i0 + nS] === 1) nS++;
  const idx = p => i0 + Math.round(p / 100 * (nS - 1));
  const sc = (o, k) => (o && o[k] && Number.isFinite(o[k].mean)) ? { mean: o[k].mean, sd: o[k].sd } : null;
  const cv = (curves, name, p) => {                  // curve point at p % of stance (curve grid is every 5 %); null if missing
    const c = curves[name]; if (!c) return null;
    const j = p / 5, m = c.mean[j], sd = c.sd[j];
    return Number.isFinite(m) ? { mean: m, sd: Number.isFinite(sd) ? sd : null } : null;
  };
  const rows = [];
  const add = (name, unit, model, ref, from) => rows.push({ name, unit, model, ref, from });
  const P9 = [10, 20, 30, 40, 50, 60, 70, 80, 90];
  // pressure point ahead of the hip: forceLine.curves.copAheadHipHjc (hip joint centre, Harrington)
  for (const p of P9) add('copAheadHip_' + p, 'm', K.copX[idx(p)], cv(fl.curves, 'copAheadHipHjc', p), 'landing.forceLine.curves.copAheadHipHjc');
  // pressure point ahead of the ankle: rollover.stanceCurves.copAheadAnkle_jc (ankle joint centre = midpoint of malleoli)
  for (const p of P9) add('copAheadAnkle_' + p, 'm', K.copX[idx(p)] - K.footX[idx(p)], cv(ro.stanceCurves, 'copAheadAnkle_jc', p), 'rollover.stanceCurves.copAheadAnkle_jc');
  // ankle ahead of hip = (pressure ahead of hip) - (pressure ahead of ankle) from the two measured means;
  // its SD is the SD of pressure-ahead-of-hip (the two curves come from different files, so no covariance is available)
  for (const p of P9) {
    const a = cv(fl.curves, 'copAheadHipHjc', p), b = cv(ro.stanceCurves, 'copAheadAnkle_jc', p);
    add('ankleAheadHip_' + p, 'm', K.footX[idx(p)], a && b ? { mean: a.mean - b.mean, sd: a.sd } : null, 'landing copAheadHipHjc - rollover copAheadAnkle_jc (SD of copAheadHipHjc)');
  }
  // hip height above the plate (hip joint centre) over mean height
  for (const p of [0, 50, 100]) {
    const r = cv(fl.curves, 'hipHeightHjc', p);
    add('hipHeightOverH_' + p, '', K.hipY[idx(p)] / H, r ? { mean: r.mean / H, sd: r.sd === null ? null : r.sd / H } : null, 'landing.forceLine.curves.hipHeightHjc / summary height');
  }
  for (const p of [0, 10, 20]) add('footDeg_' + p, 'deg', K.footAng[idx(p)] * D, cv(ro.stanceCurves, 'footDeg', p), 'rollover.stanceCurves.footDeg');
  // foot flat: first stance % where the model foot angle is within 2 deg of its value at 45 %
  const fRef = K.footAng[idx(45)] * D;
  let flat = NaN;
  for (let j = 0; j < nS; j++) if (Math.abs(K.footAng[i0 + j] * D - fRef) <= 2) { flat = 100 * j / (nS - 1); break; }
  add('footFlatAt', '% stance', flat, sc(ro.scalars, 'footFlatPct'), 'rollover.scalars.footFlatPct');
  let off = NaN;
  for (let j = 0; j < nS; j++) { const p = 100 * j / (nS - 1); if (p > 50 && K.heelZ[i0 + j] > 0.010) { off = p; break; } }
  add('heelOffAt', '% stance', off, sc(ro.scalars, 'heelOffPct'), 'rollover.scalars.heelOffPct');
  add('landDeg', 'deg', s.info.landDeg, sc(ls, 'hjc_fs_landDeg'), 'landing.cohorts.scalars.hjc_fs_landDeg');
  add('footAtStrikeDeg', 'deg', K.footAng[i0] * D, sc(ls, 'fs_footDeg'), 'landing.cohorts.scalars.fs_footDeg');
  let hmin = Infinity;
  for (let j = 0; j < nS; j++) hmin = Math.min(hmin, K.hipY[i0 + j]);
  const hd = sc(ls, 'hjc_hipDrop');   // hip drop, m: hip height at strike minus the lowest in stance (hip joint centre)
  add('hipDropOverH', '', (K.hipY[i0] - hmin) / H, hd ? { mean: hd.mean / H, sd: hd.sd / H } : null, 'landing.cohorts.scalars.hjc_hipDrop / summary height');
  return rows;
}

const UNITS = {
  contactMs: 'ms', stancePct: '%', peakFz: 'BW', peakFzAt: '% stance', peakBrake: 'BW', peakProp: 'BW', brakeImpulse: 'm/s', propImpulse: 'm/s',
  brakeEndsAt: '% stance', vertImpulseCheck: '', hipROM: 'deg', sumPos: 'J/kg',
  ankleStiffLoad: 'Nm/kg/deg', ankleEngageDeg: 'deg', ankleStiffR2: '', ankleStiffUnload: 'Nm/kg/deg', anklePeakDFAt: '% stance'
};
function unitOf(k) {
  if (UNITS[k] !== undefined) return UNITS[k];
  if (/^(knee|ankle)(At|Peak)/.test(k)) return 'deg';
  if (/Mom/.test(k)) return 'N m/kg';
  if (/(Pos|Neg)(Stance|Swing)$/.test(k)) return 'J/kg';
  if (/^fxOverFz/.test(k)) return '-';
  if (/^hipLever/.test(k)) return 'm';
  return '';
}
const fmt = x => {
  if (x === null || x === undefined || !Number.isFinite(x)) return '-';
  const t = Math.abs(x) >= 100 ? x.toFixed(1) : Math.abs(x) >= 10 ? x.toFixed(2) : x.toFixed(3);
  return /^-0\.0*$/.test(t) ? t.slice(1) : t;
};

function parseArgs(argv) {
  const a = { cohort: 'rearfoot', speed: null, set: {}, json: null };
  for (let i = 0; i < argv.length; i++) {
    const t = argv[i];
    if (t === '--cohort') a.cohort = argv[++i];
    else if (t === '--speed') a.speed = argv[++i];
    else if (t === '--json') a.json = argv[++i];
    else if (t === '--set') {
      const [k, v] = String(argv[++i]).split('=');
      if (!k || v === undefined || !Number.isFinite(Number(v))) throw new Error(`--set needs key=number, got "${argv[i]}"`);
      a.set[k] = Number(v);
    } else throw new Error(`unknown option ${t}`);
  }
  if (a.cohort !== 'all' && a.cohort !== 'rearfoot') throw new Error('--cohort must be all or rearfoot');
  return a;
}

function main() {
  const args = parseArgs(process.argv.slice(2));
  const sum = JSON.parse(fs.readFileSync(path.join(__dirname, 'data', 'fukuchi2017-summary.json'), 'utf8'));
  const speeds = Object.keys(sum.cadenceStridesPerMin).filter(sp => !args.speed || Number(sp) === Number(args.speed));
  if (!speeds.length) throw new Error(`no such speed: ${args.speed}`);
  const { heightM: H, massKg: M, men, women } = sum.subjects;
  console.log(`Model vs Fukuchi et al. 2017, cohort "${args.cohort}", ${sum.subjects.n} subjects (${men} men, ${women} women), H ${H} m, M ${M} kg`);
  if (Object.keys(args.set).length) console.log('extra model options: ' + JSON.stringify(args.set));

  // rows that can never match by construction of the model (see the note printed under each table)
  const STRUCTURAL = new Set(['hipExtMomStance', 'hipFlxMomStance', 'hipPosStance', 'hipNegStance', 'hipExtMomAll', 'anklePosSwing', 'ankleNegSwing']);
  const isStructural = k => STRUCTURAL.has(k) || k.startsWith('hipLever_');
  const land = JSON.parse(fs.readFileSync(path.join(__dirname, 'data', 'fukuchi2017-landing.json'), 'utf8'));
  const roll = JSON.parse(fs.readFileSync(path.join(__dirname, 'data', 'fukuchi2017-rollover.json'), 'utf8'));
  const totals = { ok: 0, '~': 0, OUT: 0 }, gtotals = { ok: 0, '~': 0, OUT: 0 }, report = { cohort: args.cohort, options: args.set, subjects: sum.subjects, speeds: {} };
  for (const sp of speeds) {
    const meas = sum.cohorts[args.cohort][sp];
    const o = { sex: 'male', H, M, mode: 'real', speed: Number(sp), cadence: 2 * sum.cadenceStridesPerMin[sp], share: true, elastic: 0, eff: 0.25, gravity: true,
                arcDeg: 60, phmaxDeg: 140, N: 1440, landDeg: null, contactMs: null, contactLen: null, realKnee: null, ...args.set };
    const s = simulate(o), rec = modelRecord(s, M), model = computeMetrics(rec);
    console.log(`\n=== ${sp} m/s: ${meas.n} subjects, cadence ${o.cadence.toFixed(2)} steps/min, T ${s.T.toFixed(4)} s ===`);
    console.log('metric'.padEnd(18) + 'unit'.padEnd(10) + 'model'.padStart(9) + '  ' + 'measured mean ± SD'.padEnd(22) + 'z'.padStart(7) + '  flag');
    const counts = { ok: 0, '~': 0, OUT: 0 }, rows = {};
    for (const k in model.scalars) {
      const m = model.scalars[k], ref = meas.scalars[k];
      let z = NaN, flag = '';
      if (ref && ref.sd && Number.isFinite(m) && k !== 'vertImpulseCheck') {
        z = (m - ref.mean) / ref.sd;
        if (isStructural(k)) flag = 'n/a*';          // shown, but not scored
        else { flag = Math.abs(z) <= 1 ? 'ok' : Math.abs(z) <= 2 ? '~' : 'OUT'; counts[flag]++; }
      }
      rows[k] = { model: m, mean: ref ? ref.mean : null, sd: ref ? ref.sd : null, z: Number.isFinite(z) ? z : null, flag };
      console.log(k.padEnd(18) + unitOf(k).padEnd(10) + fmt(m).padStart(9) + '  ' + (fmt(ref && ref.mean) + ' ± ' + fmt(ref && ref.sd)).padEnd(22) +
                  (Number.isFinite(z) ? z.toFixed(1) : '-').padStart(7) + '  ' + flag);
    }
    console.log('* not scored: the planted leg is treated as weightless (no stance hip moment), the swing ankle is not modelled, and the measured hip lever includes the leg\'s own inertia.');
    console.log('\nStance curves, model vs measured mean (RMSE in native units, and as % of the measured curve range):');
    const curves = {};
    for (const c of STANCE_CURVES) {
      const mc = model.stanceCurves[c], ec = meas.stanceCurves[c].mean;
      let se = 0, cnt = 0, lo = Infinity, hi = -Infinity;
      for (let i = 0; i < NSTANCE; i++) if (Number.isFinite(mc[i]) && Number.isFinite(ec[i])) { se += (mc[i] - ec[i]) ** 2; cnt++; lo = Math.min(lo, ec[i]); hi = Math.max(hi, ec[i]); }
      const rmse = Math.sqrt(se / cnt), pct = 100 * rmse / (hi - lo);
      curves[c] = { rmse, rangePct: pct, model: mc };
      console.log('  ' + c.padEnd(10) + ('RMSE ' + fmt(rmse)).padEnd(16) + pct.toFixed(1).padStart(6) + ' % of range');
    }
    const i = s.info, scalarCounts = { ...counts };

    // ---- geometry through the contact: model vs landing / rollover summaries (same cohort and speed) ----
    const gcounts = { ok: 0, '~': 0, OUT: 0 }, grows = {};
    const geo = geometryRows(s, H, land.forceLine.cohorts[args.cohort][sp], land.cohorts[args.cohort][sp].scalars, roll.cohorts[args.cohort][sp]);
    console.log(`\nGeometry through the contact at ${sp} m/s (pressure point, ankle and hip positions are forward of the hip, m; hip height above the plate, divided by H):`);
    console.log('metric'.padEnd(20) + 'unit'.padEnd(10) + 'model'.padStart(9) + '  ' + 'measured mean ± SD'.padEnd(22) + 'z'.padStart(7) + '  flag');
    for (const g of geo) {
      let z = NaN, flag = '';
      if (g.ref && Number.isFinite(g.ref.mean) && g.ref.sd && Number.isFinite(g.ref.sd) && Number.isFinite(g.model)) {
        z = (g.model - g.ref.mean) / g.ref.sd;
        flag = Math.abs(z) <= 1 ? 'ok' : Math.abs(z) <= 2 ? '~' : 'OUT';
        gcounts[flag]++;
      }
      grows[g.name] = { model: g.model, mean: g.ref ? g.ref.mean : null, sd: g.ref ? g.ref.sd : null, z: Number.isFinite(z) ? z : null, flag, from: g.from };
      console.log(g.name.padEnd(20) + g.unit.padEnd(10) + fmt(g.model).padStart(9) + '  ' + (fmt(g.ref && g.ref.mean) + ' ± ' + fmt(g.ref && g.ref.sd)).padEnd(22) +
                  (Number.isFinite(z) ? z.toFixed(1) : '-').padStart(7) + '  ' + flag);
    }
    const spr = ['ankleStiffLoad', 'ankleEngageDeg', 'anklePeakDFAt'].map(k => `${k} ${fmt(model.scalars[k])} (measured ${fmt(meas.scalars[k] && meas.scalars[k].mean)})`);
    console.log('Ankle spring, model vs measured (scored in the first table): ' + spr.join(', '));
    for (const k in counts) { counts[k] += gcounts[k]; totals[k] += counts[k]; gtotals[k] += gcounts[k]; }
    console.log(`\nCounts at ${sp} m/s: ok ${counts.ok}, ~ ${counts['~']}, OUT ${counts.OUT}`);
    console.log(`  of which scalar table: ok ${scalarCounts.ok}, ~ ${scalarCounts['~']}, OUT ${scalarCounts.OUT}; geometry table: ok ${gcounts.ok}, ~ ${gcounts['~']}, OUT ${gcounts.OUT}`);
    console.log(`Model state: contact ${(i.tc * 1000).toFixed(1)} ms, skew ${i.skew.toFixed(3)}, landDeg ${i.landDeg.toFixed(1)}, unbalanced ${i.unbalanced}`);
    report.speeds[sp] = { n: meas.n, modelState: { contactMs: i.tc * 1000, skew: i.skew, landDeg: i.landDeg, unbalanced: i.unbalanced, peakGRF: i.peakGRF }, counts, scalarCounts, geometryCounts: gcounts, scalars: rows, geometry: grows, stanceCurves: curves };
  }
  console.log(`\nOverall: ok ${totals.ok}, ~ ${totals['~']}, OUT ${totals.OUT}  (of ${totals.ok + totals['~'] + totals.OUT} flagged metrics)`);
  console.log(`  of which geometry tables: ok ${gtotals.ok}, ~ ${gtotals['~']}, OUT ${gtotals.OUT}`);
  report.totals = totals; report.geometryTotals = gtotals;
  if (args.json) { fs.writeFileSync(args.json, JSON.stringify(report, null, 1)); console.log('wrote ' + args.json); }
}

module.exports = { modelRecord, simulate };
if (require.main === module) main();
