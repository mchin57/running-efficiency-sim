// Reduces the Fukuchi 2017 processed gait files (subjects 1-28) to across-subject summary statistics:
// tests/data/fukuchi2017-summary.json. Run: node tools/reduce-fukuchi.js   (needs data/ from tools/fetch-fukuchi.js)
'use strict';
const fs = require('fs');
const path = require('path');
const { computeMetrics, STANCE_CURVES, STRIDE_CURVES } = require('../tests/lib/gait-metrics.js');

const DIR = path.join(__dirname, '..', 'data', 'fukuchi2017');
const OUT = path.join(__dirname, '..', 'tests', 'data', 'fukuchi2017-summary.json');
const SPEEDS = ['25', '35', '45'], SPEED_LABEL = { 25: '2.5', 35: '3.5', 45: '4.5' };
const CADENCE = { '2.5': 80.82, '3.5': 85.68, '4.5': 91.74 };   // strides/min, Fukuchi et al. 2017 Table 4
const NSUBJ = 28, REAR = 'Rearfoot';

function readTable(file) {
  const lines = fs.readFileSync(file, 'utf8').split(/\r?\n/).filter(l => l.trim() !== '');
  const head = lines[0].split('\t').map(s => s.trim());
  const col = {};
  head.forEach((h, i) => { col[h] = i; });
  const rows = lines.slice(1).map(l => l.split('\t').map(s => parseFloat(s)));
  return { col, rows };
}
function column(tab, name) {
  const i = tab.col[name];
  return i === undefined ? null : tab.rows.map(r => r[i]);
}

// stride record for one leg and speed, or null with a reason
function legRecord(tab, side, sp) {
  const T = 60 / CADENCE[SPEED_LABEL[sp]];
  const names = { hipAng: 'hipAngZ', kneeAng: 'kneeAngZ', ankleAng: 'ankleAngZ', hipMom: 'hipMomZ', kneeMom: 'kneeMomZ', ankleMom: 'ankleMomZ',
                  grfX: 'grfX', grfY: 'grfY', hipPow: 'hipPow', kneePow: 'kneePow', anklePow: 'anklePow' };
  const rec = { T, n: 101 };
  for (const k in names) {
    const c = column(tab, side + names[k] + sp);
    if (!c) return { why: 'missing' };
    if (c.every(Number.isNaN)) return { why: 'allNaN' };
    if (c.some(x => Math.abs(x) > 1000)) return { why: 'garbage' };
    rec[k] = c;
  }
  if (rec.hipAng.length !== 101) return { why: 'rows' };
  const jointOk = ['hipAng', 'kneeAng', 'ankleAng', 'hipMom', 'kneeMom', 'ankleMom', 'hipPow', 'kneePow', 'anklePow'].every(k => rec[k].every(Number.isFinite));
  if (!jointOk) return { why: 'nan' };
  let kTO = -1;
  while (kTO + 1 < 101 && Number.isFinite(rec.grfY[kTO + 1])) kTO++;
  if (kTO < 5 || kTO > 80) return { why: 'stance' };
  if (!rec.grfX.slice(0, kTO + 1).every(Number.isFinite)) return { why: 'grfXnan' };
  rec.kTO = kTO;
  rec.grfX = rec.grfX.map((x, i) => i <= kTO ? x : 0);
  rec.grfY = rec.grfY.map((x, i) => i <= kTO ? x : 0);
  return { rec };
}

const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
const sd = a => { if (a.length < 2) return NaN; const m = mean(a); return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / (a.length - 1)); };
const r4 = x => Number.isFinite(x) ? Number(x.toPrecision(4)) : null;

function main() {
  const infoLines = fs.readFileSync(path.join(DIR, 'RBDSinfo.txt'), 'utf8').split(/\r?\n/).filter(Boolean);
  const head = infoLines[0].split('\t'), ix = n => head.indexOf(n);
  const subj = {};
  for (const l of infoLines.slice(1)) {
    const f = l.split('\t'), id = parseInt(f[ix('Subject')], 10);
    if (id >= 1 && id <= NSUBJ && !subj[id]) subj[id] = f;
  }
  const ids = Object.keys(subj).map(Number).sort((a, b) => a - b);
  if (ids.length !== NSUBJ) console.warn(`warning: ${ids.length} subjects in RBDSinfo, expected ${NSUBJ}`);
  const heights = ids.map(i => parseFloat(subj[i][ix('Height')]) / 100), masses = ids.map(i => parseFloat(subj[i][ix('Mass')]));
  const men = ids.filter(i => subj[i][ix('Gender')].trim() === 'M').length;

  const fsValues = new Set(), skipped = {};
  // per cohort, per speed: list of per-subject results
  const cohorts = { all: {}, rearfoot: {} };
  for (const c in cohorts) for (const sp of SPEEDS) cohorts[c][sp] = [];
  const legCount = { all: {}, rearfoot: {} };
  for (const c in legCount) for (const sp of SPEEDS) legCount[c][sp] = 0;

  for (const id of ids) {
    const tab = readTable(path.join(DIR, `RBDS${String(id).padStart(3, '0')}processed.txt`));
    for (const sp of SPEEDS) {
      const per = { all: [], rearfoot: [] };
      for (const side of ['R', 'L']) {
        const fs0 = (subj[id][ix(side + 'FSI' + sp)] || '').trim();
        if (fs0) fsValues.add(fs0);
        const res = legRecord(tab, side, sp);
        if (!res.rec) { const k = res.why; skipped[k] = (skipped[k] || 0) + 1; continue; }
        const m = computeMetrics(res.rec);
        per.all.push(m);
        if (fs0 === REAR) per.rearfoot.push(m);
      }
      for (const c in per) if (per[c].length) { cohorts[c][sp].push(per[c]); legCount[c][sp] += per[c].length; }
    }
  }

  // average the legs of a subject, then mean / SD / n across subjects
  const avgLegs = legs => {
    const out = { scalars: {}, stanceCurves: {}, strideCurves: {} };
    const avg = arrs => arrs[0].map((_, j) => { const v = arrs.map(a => a[j]).filter(Number.isFinite); return v.length ? mean(v) : NaN; });
    for (const k in legs[0].scalars) { const v = legs.map(l => l.scalars[k]).filter(Number.isFinite); out.scalars[k] = v.length ? mean(v) : NaN; }
    for (const k of STANCE_CURVES) out.stanceCurves[k] = avg(legs.map(l => l.stanceCurves[k]));
    for (const k of STRIDE_CURVES) out.strideCurves[k] = avg(legs.map(l => l.strideCurves[k]));
    return out;
  };
  const stat = vals => { const v = vals.filter(Number.isFinite); return { mean: r4(mean(v)), sd: r4(sd(v)), n: v.length }; };
  const curveStat = subjects => {
    const m = [], s = [];
    for (let j = 0; j < subjects[0].length; j++) { const v = subjects.map(c => c[j]).filter(Number.isFinite); m.push(r4(mean(v))); s.push(r4(sd(v))); }
    return { mean: m, sd: s };
  };

  const summary = {
    source: { citation: 'Fukuchi RK, Fukuchi CA, Duarte M (2017). A public dataset of running biomechanics and the effects of running speed on lower extremity kinematics and kinetics. PeerJ 5:e3298',
              doi: '10.7717/peerj.3298', figshare: 'https://doi.org/10.6084/m9.figshare.4543435 (article 4543435)', licence: 'CC BY 4.0', note: 'derived summary statistics' },
    generated: new Date().toISOString(),
    subjects: { n: ids.length, heightM: r4(mean(heights)), massKg: r4(mean(masses)), men, women: ids.length - men },
    cadenceStridesPerMin: CADENCE,
    cohorts: { all: {}, rearfoot: {} }
  };
  for (const c in cohorts) for (const sp of SPEEDS) {
    const subs = cohorts[c][sp].map(avgLegs), key = SPEED_LABEL[sp];
    if (!subs.length) { summary.cohorts[c][key] = { n: 0, scalars: {}, stanceCurves: {}, strideCurves: {} }; continue; }
    const e = { n: subs.length, scalars: {}, stanceCurves: {}, strideCurves: {} };
    for (const k in subs[0].scalars) e.scalars[k] = stat(subs.map(x => x.scalars[k]));
    for (const k of STANCE_CURVES) e.stanceCurves[k] = curveStat(subs.map(x => x.stanceCurves[k]));
    for (const k of STRIDE_CURVES) e.strideCurves[k] = curveStat(subs.map(x => x.strideCurves[k]));
    summary.cohorts[c][key] = e;
  }
  fs.mkdirSync(path.dirname(OUT), { recursive: true });
  fs.writeFileSync(OUT, JSON.stringify(summary));

  console.log(`subjects ${ids.length} (${men} men, ${ids.length - men} women), mean height ${summary.subjects.heightM} m, mass ${summary.subjects.massKg} kg`);
  console.log('foot-strike values seen: ' + [...fsValues].join(', ') + `  (rearfoot cohort = "${REAR}")`);
  console.log('legs skipped: ' + (Object.keys(skipped).length ? JSON.stringify(skipped) : 'none'));
  const show = ['contactMs', 'stancePct', 'peakFz', 'peakFzAt', 'peakBrake', 'brakeImpulse', 'propImpulse', 'vertImpulseCheck', 'kneePeakSwing', 'kneeExtMom', 'anklePFMom', 'hipExtMomAll', 'hipFlxMomAll', 'sumPos'];
  for (const c of ['all', 'rearfoot']) {
    console.log(`\ncohort ${c}   (subjects/legs per speed: ` + SPEEDS.map(sp => `${SPEED_LABEL[sp]}: ${summary.cohorts[c][SPEED_LABEL[sp]].n}/${legCount[c][sp]}`).join(', ') + ')');
    console.log('metric'.padEnd(18) + SPEEDS.map(sp => (SPEED_LABEL[sp] + ' m/s').padStart(22)).join(''));
    for (const k of show) console.log(k.padEnd(18) + SPEEDS.map(sp => { const s = summary.cohorts[c][SPEED_LABEL[sp]].scalars[k]; return (s ? `${s.mean} ± ${s.sd}` : '-').padStart(22); }).join(''));
  }
  console.log('\nwrote ' + path.relative(process.cwd(), OUT));
}

main();
