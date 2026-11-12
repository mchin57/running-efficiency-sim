// Sanity checks for the physics engine. Run with: node tests/check.js
const fs = require('fs');
const path = require('path');
eval(fs.readFileSync(path.join(__dirname, '..', 'src', 'physics.js'), 'utf8'));
const g = 9.81;   // the engine's own G is block-scoped inside the eval

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
near('realistic swing, W', r.mech, 144.6, 2);
near('realistic stance, W', r.stance.mech, 289.0, 3);
near('contact time, ms', r.info.tc * 1000, 247.3, 2);
near('peak ground force, x body weight', r.info.peakGRF, 2.40, 0.02, 2);
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

// Every combination of settings stays finite, balanced and sensible
let bad = 0, worstBal = 0, worstKnee = 0;
for (const v of [2, 3, 4.5, 6]) for (const cadence of [140, 180, 220]) for (const landDeg of [null, 0, 15, 30])
  for (const contactMs of [null, 60, 250, 400]) for (const realKnee of [null, 40, 140]) {
    const s = simulate({ ...man, mode: 'real', speed: v, cadence, landDeg, contactMs, realKnee, N: 480 });
    if (!isFinite(s.mech) || !isFinite(s.stance.mech)) bad++;
    worstBal = Math.max(worstBal, s.stance.balanceErr);
    for (const ph of s.K.ph) worstKnee = Math.min(worstKnee, ph);
  }
if (bad) failures++;
console.log(`${bad ? 'FAIL' : 'ok  '} all realistic-stride setting combinations finite (${bad} bad)`);
below('worst stance energy balance error over all settings, %', worstBal * 100, 0.01);
below('worst knee hyperextension over all settings, degrees', -worstKnee * 180 / Math.PI, 6, 1);

process.exit(failures ? 1 : 0);
