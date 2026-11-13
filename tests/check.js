// Sanity checks for the physics engine. Run with: node tests/check.js
const fs = require('fs');
const path = require('path');
// load the engine as a function body (a direct eval would work too, but V8 can't optimise code inside one: ~10x slower)
const { simulate, realisticPlan, legModel, typicalCadence, G: g } = new Function(
  fs.readFileSync(path.join(__dirname, '..', 'src', 'physics.js'), 'utf8') + '\nreturn { simulate, realisticPlan, legModel, typicalCadence, G };')();

const man = { sex: 'male', H: 1.75, M: 90.3, cadence: 180, arcDeg: 60, phmaxDeg: 140, gravity: true, share: true, elastic: 0, eff: 0.25 };
let failures = 0;
function near(label, got, want, tol, digits = 1) {
  const ok = Math.abs(got - want) <= tol;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${label}: ${got.toFixed(digits)} (expected ${want} ± ${tol})`);
}
function below(label, got, limit, digits = 4) {
  const ok = got <= limit;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${label}: ${got.toFixed(digits)} (at most ${limit})`);
}

// Simple pendulum modes (hip fixed in the air)
near('straight leg swing, W', simulate({ ...man, mode: 'straight' }).mech, 187.5, 0.5);
near('bend then straighten, W', simulate({ ...man, mode: 'fold' }).mech, 1025.9, 0.5);
near('fold, each joint pays its own way, W', simulate({ ...man, mode: 'fold', share: false }).mech, 1268.0, 0.5);
{ // the fold has no jump in knee torque where it ends: the step there is no bigger than anywhere else in the stride
  const s = simulate({ ...man, mode: 'fold' }), N = s.N, h = N / 2, step = i => Math.abs(s.S.tauK[(i + 1) % N] - s.S.tauK[i]);
  let elsewhere = 0; for (let i = 0; i < N; i++) if (i !== h - 1 && i !== N - 1) elsewhere = Math.max(elsewhere, step(i));
  below('fold: knee torque step at the end of the fold vs largest step elsewhere, N m', step(h - 1) - elsewhere, 0.01, 2);
}

// Realistic stride at 3 m/s
const r = simulate({ ...man, mode: 'real', speed: 3 });
near('realistic swing, W', r.mech, 160.9, 2);
near('realistic stance, W', r.stance.mech, 284.0, 3);
near('contact time, ms', r.info.tc * 1000, 293.4, 2);
near('peak ground force, x body weight', r.info.peakGRF, 1.97, 0.02, 2);
// joint power equals the ground force's power on the body (the planted foot point never moves)
below('stance energy balance error, %', r.stance.balanceErr * 100, 0.01);
// steady running: the body ends each step with the energy it started with, so stance work nets to zero
below('net stance work per step, J', Math.abs(r.stance.net), 0.05);
// braking and propulsion cancel through the force timing, not by tilting the force off the hip
below('residual tilt of the ground force', Math.abs(r.stance.tilt), 1e-4, 6);
below('stance hip work per step, J', r.stance.hipPos + r.stance.hipNeg, 0.01);
{ // one body motion: the drawn hip's vertical acceleration matches the force the swing leg feels
  const K = r.K, N = r.N, dt = r.dt; let worst = 0, peak = 0;
  for (let i = 0; i < N; i++) {
    const n = (i + 1) % N, p = (i - 1 + N) % N;
    const acc = (K.hipY[n] - 2 * K.hipY[i] + K.hipY[p]) / (dt * dt);
    worst = Math.max(worst, Math.abs(acc - (K.gEff[i] - g))); peak = Math.max(peak, K.gEff[i]);
  }
  below('hip vertical acceleration vs ground force, share of peak', worst / peak, 0.005);
}
{ // the planted leg reaches the foot exactly, and the foot is rigid
  const K = r.K, L = r.leg; let gap = 0, foot = 0;
  for (let i = 0; i < r.N; i++) if (K.stance[i]) {
    const ps = K.th[i] - K.ph[i];
    gap = Math.max(gap, Math.hypot(L.L1 * Math.sin(K.th[i]) + L.Ls * Math.sin(ps) - K.footX[i], K.hipY[i] - L.L1 * Math.cos(K.th[i]) - L.Ls * Math.cos(ps) - K.ankleZ[i]));
    foot = Math.max(foot, Math.abs(Math.hypot(K.ballX[i] - K.footX[i], K.ballZ[i] - K.ankleZ[i]) - Math.hypot(0.45 * L.footLen, L.ha)));
  }
  below('leg-to-ankle gap in stance, m', gap, 1e-9, 12);
  below('change in ankle-to-ball distance, m', foot, 1e-9, 12);
}

// Swing grows with speed as in measured running (Clark, Meng & Stearne 2020): peak thigh flexion 0.102 v + 0.273 rad,
// and at the measured step rate (0.31 v + 1.58 steps/s) swing cost rises with every step up in speed
{
  let prev = 0, rising = true;
  for (const v of [3, 4, 5, 6]) {
    const s = simulate({ ...man, mode: 'real', speed: v, cadence: 60 * (0.31 * v + 1.58) });
    let flex = -1; for (let i = 0; i < s.N; i++) if (!s.K.stance[i]) flex = Math.max(flex, s.K.th[i]);
    near(`peak swing thigh flexion at ${v} m/s, degrees`, flex * 180 / Math.PI, +((0.102 * v + 0.273) * 180 / Math.PI).toFixed(1), 0.5);
    if (s.mech <= prev) rising = false; prev = s.mech;
  }
  if (!rising) failures++;
  console.log(`${rising ? 'ok  ' : 'FAIL'} swing cost rises with speed at the measured step rate (6 m/s: ${prev.toFixed(0)} W)`);
}

// The swing leg starts with exactly the speed the stance leg ended with (no free energy at toe-off)
for (const v of [3, 6]) {
  const o = { ...man, mode: 'real', speed: v, cadence: 180 }, s = simulate(o);
  const P = realisticPlan(o, legModel('male', 1.75, 90.3), 120 / 180), hh = 1e-5;
  const stanceEnd = (P.stanceAt(P.tc).th - P.stanceAt(P.tc - hh).th) / hh;
  near(`thigh speed: swing start minus stance end at ${v} m/s, deg/s`, (s.K.thd[0] - stanceEnd) * 180 / Math.PI, 0, 2);
}
{ // slow running (2 m/s at its typical cadence): the heel waits for the knee's deepest bend before lifting, so the centre of pressure moves smoothly
  const s = simulate({ ...man, mode: 'real', speed: 2, cadence: typicalCadence(2) }), K = s.K; let jump = 0;
  for (let i = 1; i < s.N; i++) if (K.stance[i] && K.stance[i - 1]) jump = Math.max(jump, Math.abs(K.copX[i] - K.copX[i - 1] + 2 * s.dt));
  below('2 m/s at its typical cadence: largest centre-of-pressure step along the foot between samples, cm', jump * 100, 1, 2);
  below('2 m/s at its typical cadence: stance hip work per step, J', s.stance.hipPos + s.stance.hipNeg, 0.01);
}
{ // touchdown knee: the typical stride lands with the typical bend; a leg that would lock straight lands with more
  const d = simulate({ ...man, mode: 'real', speed: 3 }).info;
  near('default stride: touchdown knee bend, degrees', d.kneeTDDeg, 10, 1e-9, 3);
  console.log(`${d.kneeTDRaised === false ? 'ok  ' : 'FAIL'} default stride: touchdown knee bend is not raised (${d.kneeTDRaised})`);
  if (d.kneeTDRaised !== false) failures++;
  const H = 1.75, s = simulate({ ...man, H, M: 25 * H * H, mode: 'real', speed: 4.5, cadence: 220, landDeg: 0 });
  const raised = s.info.kneeTDRaised === true;
  if (!raised) failures++;
  console.log(`${raised ? 'ok  ' : 'FAIL'} straight-landing, short stride (4.5 m/s, 220 steps/min, landing 0): touchdown knee bend is raised (${s.info.kneeTDDeg.toFixed(1)} deg)`);
  below('  ... and the leg still reaches the foot (reachShort), m', s.info.reachShort, 0, 6);
  below('  ... and the energy balance error', s.stance.balanceErr, 1e-4, 8);
  const f = simulate({ ...man, H, M: 25 * H * H, mode: 'real', speed: 6, cadence: 140, landDeg: 30, contactMs: 250 });
  let mk = 1e9; for (const ph of f.K.ph) mk = Math.min(mk, ph);
  const okk = mk * 180 / Math.PI >= -6;
  if (!okk) failures++;
  console.log(`${okk ? 'ok  ' : 'FAIL'} far-forward landing at 6 m/s, 140 steps/min: smallest knee angle in the stride, degrees: ${(mk * 180 / Math.PI).toFixed(1)} (at least -6)`);
}

// Toe roll and contact length: with the foot rolling off the toe tip, contact gets close to measured lengths
// (contact time 0.78 v^-0.89 s; Clark, Meng & Stearne 2020). It matches at 3 m/s and at 6 m/s.
for (const [v, tol] of [[3, 15], [6, 3]]) {
  const s = simulate({ ...man, mode: 'real', speed: v }), want = 780 * v ** -0.89;
  below(`contact time at ${v} m/s, % short of measured ${want.toFixed(0)} ms`, 100 * (1 - s.info.tc * 1000 / want), tol, 1);
}
{ // contact length and contact time are the same setting: either gives the same stride
  const a = simulate({ ...man, mode: 'real', speed: 4, contactMs: 200 }), b = simulate({ ...man, mode: 'real', speed: 4, contactLen: 0.8 });
  near('contact set as 0.8 m at 4 m/s vs as 200 ms: stance power difference, W', b.stance.mech - a.stance.mech, 0, 0.01, 3);
  const s = simulate({ ...man, mode: 'real', speed: 6 });
  below('toe joint is part of the energy balance at 6 m/s (toe roll in use), %', s.stance.balanceErr * 100, 0.01);
  near('typical cadence at 3 m/s, steps/min (recreational runners, Fukuchi et al. 2017)', typicalCadence(3), 166.7, 0.1);
}
{ // the toes carry no moment until the pressure is on them (in front of the ball)
  const s = simulate({ ...man, mode: 'real', speed: 3 }), K = s.K; let m = 0;
  for (let i = 0; i < s.N; i++) if (K.stance[i] && K.copX[i] <= K.ballX[i]) m = Math.max(m, Math.abs(s.S.tauM[i]));
  below('toe-joint moment while the pressure is behind the ball, N m', m, 1e-9, 3);
}
{ // the heel-lift start moves smoothly with contact time (it used to snap to a grid, breaking the balance search)
  const o = { ...man, sex: 'male', H: 1.5, M: 56.25, mode: 'real', speed: 2, cadence: 220, landDeg: 0 }, L = legModel('male', 1.5, 56.25);
  const starts = [114, 115, 116, 117].map(ms => realisticPlan({ ...o, contactMs: ms }, L, 120 / 220).liftStart);
  const distinct = new Set(starts.map(x => x.toFixed(9))).size;
  if (distinct < 4) failures++;
  console.log(`${distinct < 4 ? 'FAIL' : 'ok  '} heel-lift start changes smoothly with contact time (${distinct} distinct of 4)`);
  near('that setting at 115.5 ms balances (skew solves)', realisticPlan({ ...o, contactMs: 115.5 }, L, 120 / 220).unbalanced ? 1 : 0, 0, 0, 0);
  const d = simulate({ ...man, mode: 'real', speed: 3 }).info;
  // (the longer contact is now limited by the flight time running out ('air'), the shorter one by balance)
  if (d.minWhy !== 'balance' || d.maxWhy !== 'air') failures++;
  console.log(`${d.minWhy === 'balance' && d.maxWhy === 'air' ? 'ok  ' : 'FAIL'} default contact range: shortest set by balance, longest by flight time (${d.minWhy}, ${d.maxWhy})`);
}
{ // at the cadence measured for each speed, the typical stride balances (no warning) from 2 to 6 m/s
  let flagged = [];
  for (const v of [2, 3, 4, 5, 6]) if (simulate({ ...man, mode: 'real', speed: v, cadence: typicalCadence(v), N: 480 }).info.unbalanced) flagged.push(v);
  if (flagged.length) failures++;
  console.log(`${flagged.length ? 'FAIL' : 'ok  '} typical strides at measured cadences balance (${flagged.length ? 'flagged at ' + flagged.join(', ') + ' m/s' : 'all'})`);
}

// Every combination of settings stays finite, balanced and sensible
let bad = 0, worstBal = 0, worstKnee = 0, worstSwing = 0, worstTorque = 0, minKneeTO = 99, worstHip = 0, unbalanced = 0, n = 0;
for (const sex of ['male', 'female']) for (const H of [1.5, 1.75, 2.0]) for (const v of [2, 3, 4.5, 6]) for (const cadence of [140, 180, 220])
  for (const landDeg of [null, 0, 15, 30]) for (const contactMs of [null, 60, 250, 400]) for (const realKnee of [null, 40, 140]) {
    if (H !== 1.75 && (realKnee !== null || contactMs !== null)) continue;   // body size: a smaller sweep
    const s = simulate({ ...man, sex, H, M: 25 * H * H, mode: 'real', speed: v, cadence, landDeg, contactMs, realKnee, N: 480 }); n++;
    if (!isFinite(s.mech) || !isFinite(s.stance.mech)) bad++;
    worstBal = Math.max(worstBal, s.stance.balanceErr);
    for (const ph of s.K.ph) worstKnee = Math.min(worstKnee, ph);
    worstSwing = Math.max(worstSwing, s.mech / (25 * H * H));
    minKneeTO = Math.min(minKneeTO, s.info.kneeTODeg);
    if (s.info.unbalanced) { unbalanced++; continue; }
    worstHip = Math.max(worstHip, s.stance.hipPos + s.stance.hipNeg);
    for (let i = 0; i < s.N; i++) if (s.K.count[i]) worstTorque = Math.max(worstTorque, Math.abs(s.S.tauH[i]) / (25 * H * H));
  }
if (bad) failures++;
console.log(`${bad ? 'FAIL' : 'ok  '} all ${n} realistic-stride setting combinations finite (${bad} bad)`);
below('worst stance energy balance error over all settings, %', worstBal * 100, 0.01);
below('worst knee hyperextension over all settings, degrees', -worstKnee * 180 / Math.PI, 6, 1);
// the leg always reaches the ground at push-off with some knee bend (the landing moves forward if it has to)
below('straightest knee at toe-off over all settings, degrees below 10', 10 - minKneeTO, 0.05, 2);
// extreme settings (e.g. 140 degree knee bend at 220 steps/min) are genuinely costly; this catches blow-ups, not them
below('worst swing power over all settings, W/kg', worstSwing, 15, 1);
// where braking and propulsion can cancel, the force passes through the hip; elsewhere the page flags the setting
below('worst swing hip torque where the contact balances, N m/kg', worstTorque, 10, 1);   // e.g. a 2 m runner at 6 m/s, 220 steps/min
below('worst stance hip work where the contact balances, J', worstHip, 0.5, 2);
console.log(`     (${unbalanced} of ${n} combinations can't balance and are flagged on the page)`);

process.exit(failures ? 1 : 0);
