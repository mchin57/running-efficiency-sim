// Sanity checks for the physics engine. Run with: node tests/check.js
const fs = require('fs');
const path = require('path');
eval(fs.readFileSync(path.join(__dirname, '..', 'src', 'physics.js'), 'utf8'));

const man = { sex: 'male', H: 1.75, M: 90.3, cadence: 180, arcDeg: 60, phmaxDeg: 140, gravity: true, share: true, elastic: 0, eff: 0.25 };
let failures = 0;
function near(label, got, want, tol) {
  const ok = Math.abs(got - want) <= tol;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${label}: ${got.toFixed(1)} (expected ${want} ± ${tol})`);
}

// Simple pendulum modes (hip fixed in the air)
near('straight leg swing, W', simulate({ ...man, mode: 'straight' }).mech, 187.4, 0.5);
near('bend then straighten, W', simulate({ ...man, mode: 'fold' }).mech, 867.9, 0.5);
near('fold, each joint pays its own way, W', simulate({ ...man, mode: 'fold', share: false }).mech, 1061.5, 0.5);

// Realistic stride at 3 m/s
const r = simulate({ ...man, mode: 'real', speed: 3 });
near('realistic swing, W', r.mech, 146, 2);
near('realistic stance, W', r.stance.mech, 325, 3);
near('contact time, ms', r.info.tc * 1000, 247, 2);
near('peak ground force, x body weight', r.info.peakGRF, 2.12, 0.02);
near('stance energy balance error, %', r.stance.balanceErr * 100, 4, 2);

// Every combination of settings stays finite
let bad = 0;
for (const v of [2, 3, 4.5, 6]) for (const cadence of [140, 180, 220]) for (const landDeg of [null, 0, 15, 30])
  for (const contactMs of [null, 60, 250, 400]) for (const realKnee of [null, 40, 140]) {
    const s = simulate({ ...man, mode: 'real', speed: v, cadence, landDeg, contactMs, realKnee, N: 480 });
    if (!isFinite(s.mech) || !isFinite(s.stance.mech)) bad++;
  }
if (bad) failures++;
console.log(`${bad ? 'FAIL' : 'ok  '} all realistic-stride setting combinations finite (${bad} bad)`);

process.exit(failures ? 1 : 0);
