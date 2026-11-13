// ---- Leg-swing physics: two-segment leg hung from the hip, prescribed motion, inverse dynamics ----
const G = 9.81;
const DRAG_TOL = 0.01;   // m: how far the simplified foot may dip into the ground before it counts as dragging
// de Leva (1996) segment parameters: [mass fraction of body, COM from proximal joint (fraction of length), sagittal radius of gyration (fraction of length)]
const SEG = {
  male:   { th: [0.1416, 0.4095, 0.329], sh: [0.0433, 0.4459, 0.255], ft: 0.0137 },
  female: { th: [0.1478, 0.3612, 0.369], sh: [0.0481, 0.4416, 0.271], ft: 0.0129 }
};

function legModel(sex, H, M) {
  const p = SEG[sex];
  const L1 = 0.245 * H, Ls = 0.246 * H, ha = 0.039 * H, footLen = 0.152 * H; // Winter length ratios
  const m1 = p.th[0] * M, c1 = p.th[1] * L1, I1 = m1 * (p.th[2] * L1) ** 2;
  const ms = p.sh[0] * M, mf = p.ft * M, dfoot = Ls + ha / 2;
  const m2 = ms + mf, c2 = (ms * p.sh[1] * Ls + mf * dfoot) / m2;
  const I2 = ms * (p.sh[2] * Ls) ** 2 + ms * (p.sh[1] * Ls - c2) ** 2 + mf * (dfoot - c2) ** 2;
  return { L1, Ls, ha, footLen, m1, c1, I1, m2, c2, I2, legMass: m1 + m2, legLen: L1 + Ls + ha };
}

// ---------- simple modes: hip pinned, theta = -A cos(wt), forward swing for 0<t<T/2 ----------
function simpleKinematics(o, T, N) {
  const w = 2 * Math.PI / T, A = (o.arcDeg / 2) * Math.PI / 180, phmax = o.phmaxDeg * Math.PI / 180;
  const K = blankK(N);
  for (let i = 0; i < N; i++) {
    const t = i * T / N;
    K.th[i] = -A * Math.cos(w * t); K.thd[i] = A * w * Math.sin(w * t); K.thdd[i] = A * w * w * Math.cos(w * t);
    if (o.mode === 'bent') K.ph[i] = phmax;
    else if (o.mode === 'fold' && t < T / 2) {
      // knee bend phmax * 64 s^3 (1-s)^3: bend, speed and acceleration are all zero at both ends of the forward swing,
      // so the knee torque has no jumps where the fold starts and ends
      const Th = T / 2, s = t / Th, g = s * (1 - s), g1 = 1 - 2 * s;
      K.ph[i] = phmax * 64 * g ** 3;
      K.phd[i] = phmax * 192 * g * g * g1 / Th;
      K.phdd[i] = phmax * 384 * g * (g1 * g1 - g) / (Th * Th);
    }
    K.gEff[i] = o.gravity ? G : 0; K.count[i] = 1;
  }
  return K;
}
function blankK(N) {
  const z = () => new Float64Array(N);
  return { th: z(), thd: z(), thdd: z(), ph: z(), phd: z(), phdd: z(), gEff: z(), count: z(), stance: z(),
           hipY: z(), hipVy: z(), grf: z(), footPitch: z(), copX: z(),
           footX: z(), ankleZ: z(), footAng: z(), footAngRate: z(), toeAngRate: z(), heelX: z(), heelZ: z(), ballX: z(), ballZ: z(), toeX: z(), toeZ: z() };
}

// ---------- realistic stride: real ground, planted foot that rolls heel to toe, bouncing body ----------
// The body (lumped at the hip) is driven by the ground force. Per unit body mass it is
//   F/m = Fpk [sin(pi s) + kappa sin(2 pi s)]   during contact (s = fraction of the contact), 0 in the air.
// The kappa term moves the force earlier or later in the contact without changing its impulse. kappa is chosen so that,
// with the force pointing from the centre of pressure at the hip, braking and propulsion cancel (steady speed).
// The hip height follows from that force in closed form, so the swing leg, the stance leg and the drawing all share one
// body motion. During contact the foot is planted: it lands on the heel with the toes up and lowers flat, stays flat,
// then the heel lifts as the rearfoot turns about the ball (toe joint) while the toes stay flat and the pressure moves out
// to the toe tip, and finally, if the leg needs more reach, the toes roll up about their tip. The knee and thigh angles
// are whatever geometry makes the leg reach from the hip down to the ankle. The foot lowers at the rate that makes the
// knee bend at a typical speed just after touchdown; later the heel lifts just enough to stop the knee straightening past
// a floor that ends near 22 degrees at toe-off.
// Typical values for a speed v (m/s) come from measured running. Clark, Meng & Stearne 2020 (154 trials at 3.1-10 m/s):
// contact time 0.78 v^-0.89 s (scaled with height from a 1.75 m runner), peak swing thigh flexion 0.102 v + 0.273 rad.
// The typical landing point is measured directly from the hip and ankle markers of 28 runners at 2.5-4.5 m/s (Fukuchi et al.
// 2017 public dataset; tools/landing-fukuchi.js): 13.9 degrees at 3.5 m/s, 3.3 degrees per m/s, continued above 4.5 m/s at the
// slope implied by Clark's thigh fit (1.73 degrees per m/s). Clark's thigh-at-touchdown fit is no longer used because it runs
// 2-7 degrees above these runners. Touchdown knee bend (12 degrees), foot-flat time, heel-lift time and the pressure travel
// follow the heel strikers among those runners (see tests/validate.js). Cadence: a straight-line fit, 133.9 + 10.9 v steps/min, to 28 recreational
// runners at 2.5, 3.5 and 4.5 m/s (Fukuchi et al. 2017, Table 4: 80.8, 85.7 and 91.7 strides/min).
const REAL = { kneeTD: 12, kneeTO: 22, kneeMinTO: 10, kneeRateTD: 300, kneeRateTOMax: 600,
               strikeDeg: 12, flatAt: 0.12, rollStart: 0.55, liftKnee: 30, copMin: 0.10, copEase: 1.5, kneeFloorTD: 3, kneeTDMax: 30, kneeRateTDSwingMax: 450, holdAt: 0.85, rollMax: 65, skewMax: 0.5,
               toeLen: 0.30, tipAt: 0.6, toeRollMax: 45,           // toes are 30% of the foot; pressure reaches the tip 60% into the lift
               landing: { at: 3.5, deg: 13.9, slope: 3.3, kneeSpeed: 4.5, slopeFast: 1.73, min: 5 }, contactTime: [0.78, -0.89], refHeight: 1.75, cadence: [10.92, 133.9],
               thighFlex: [0.102, 0.273], retractMin: 4,   // peak swing thigh flexion = 0.102 v + 0.273 rad (v in m/s)
               termKnee: 14, clearBoostMax: 15 };          // knee bend at 88% of swing; most extra thigh swing for foot clearance
const sstep = s => { s = Math.min(1, Math.max(0, s)); return s * s * s * (s * (6 * s - 15) + 10); };
// typical cadence (steps per minute, both feet) at speed v
function typicalCadence(v) { return REAL.cadence[0] * v + REAL.cadence[1]; }
// typical landing point (degrees from vertical, hip to ankle at touchdown) at speed v
function typicalLanding(v) { const L = REAL.landing, base = L.deg + L.slope * (Math.min(v, L.kneeSpeed) - L.at); return Math.max(L.min, base + L.slopeFast * Math.max(0, v - L.kneeSpeed)); }

// Everything about the stride that does not depend on the late-swing knee bend: contact timing, ground force,
// body motion and the stance-leg geometry.
function planFor(o, leg, T, kneeTDdeg = REAL.kneeTD) {
  const v = o.speed, Ts = T / 2, D = Math.PI / 180;
  const { L1, Ls, ha, footLen } = leg;
  // ankle-to-ball, heel-to-ankle, heel-to-ball and ball-to-toe-tip distances along the sole
  const ball = 0.45 * footLen, hx = 0.25 * footLen, sole = hx + ball, Lt = REAL.toeLen * footLen;
  const reach = ph => Math.sqrt(L1 * L1 + Ls * Ls + 2 * L1 * Ls * Math.cos(ph));      // hip-to-ankle distance at knee bend ph
  const kneeFor = d => Math.acos(Math.max(-1, Math.min(1, (d * d - L1 * L1 - Ls * Ls) / (2 * L1 * Ls))));
  // ankle position relative to the heel when the toes are up by b, and relative to the ball when the heel is up by al
  // (rigid rearfoot turning about the point on the ground), with their rates of change per radian
  const heelOff = b => [hx * Math.cos(b) - ha * Math.sin(b), hx * Math.sin(b) + ha * Math.cos(b)];
  const heelOffD = b => [-hx * Math.sin(b) - ha * Math.cos(b), hx * Math.cos(b) - ha * Math.sin(b)];
  const ballOff = al => [-ball * Math.cos(al) + ha * Math.sin(al), ball * Math.sin(al) + ha * Math.cos(al)];
  // ball (toe joint) relative to the toe tip when the toes are rolled up by g
  const toeOff = g => [-Lt * Math.cos(g), Lt * Math.sin(g)];
  // landing point: angle from vertical of the line from the hip to the ground under the ankle at touchdown. Typical: the
  // measured landing point (typicalLanding) with the knee at the typical touchdown bend.
  const b0 = REAL.strikeDeg * D, strike = heelOff(b0);
  // The typical landing point always uses the typical touchdown knee (REAL.kneeTD), even when this plan lands with more bend.
  const dTyp = reach(REAL.kneeTD * D);
  const hipAtTypical = x => Math.sqrt(Math.max(dTyp * dTyp - x * x, 1e-6)) + strike[1];
  const landTypical = typicalLanding(v);
  // the ankle's horizontal distance ahead of the hip that gives that angle with the typical touchdown posture
  let aTyp = 0.2;
  for (let it = 0; it < 60; it++) aTyp = Math.min(0.9 * dTyp, Math.tan(landTypical * D) * hipAtTypical(aTyp));
  const landWanted = o.landDeg != null ? o.landDeg : landTypical;
  const dTD = reach(kneeTDdeg * D);                          // the real touchdown posture of this plan
  const hipAtTD = x => Math.sqrt(Math.max(dTD * dTD - x * x, 1e-6)) + strike[1];
  let a = Math.min(aTyp, 0.9 * dTD);
  for (let it = 0; it < 60; it++) a = Math.min(0.9 * dTD, Math.tan(landWanted * D) * hipAtTD(a));
  let yTD = hipAtTD(a);                                     // hip height at touchdown: heel down, toes up, knee at 10 degrees
  // ground points measured from the point on the ground under the hip, tau after touchdown
  const heelAt = tau => a - strike[0] - v * tau, ballAt = tau => heelAt(tau) + sole, tipAt = tau => ballAt(tau) + Lt;
  const tFlat = REAL.flatAt;
  // vertical body motion in contact, in closed form. The touchdown velocity is the one that makes the motion repeat every step.
  const motion = (tc, kap) => {
    const tf = Ts - tc, k = Math.PI / tc, Fpk = Math.PI * G * Ts / (2 * tc);   // Fpk: supports body weight over each step
    const vTD = -G * tf / 2 - G * tc * kap / 4, vTO = vTD + G * tf;
    return { vTD, vTO, Fpk,
      y: tau => yTD + vTD * tau + Fpk * ((tau - Math.sin(k * tau) / k) / k + kap * (tau - Math.sin(2 * k * tau) / (2 * k)) / (2 * k)) - G * tau * tau / 2,
      vy: tau => vTD + Fpk * ((1 - Math.cos(k * tau)) / k + kap * (1 - Math.cos(2 * k * tau)) / (2 * k)) - G * tau,
      f: tau => Fpk * (Math.sin(k * tau) + kap * Math.sin(2 * k * tau)) };   // ground force per unit mass
  };
  // knee bend with the foot flat
  const kneeFlat = (mo, tau) => kneeFor(Math.hypot(heelAt(tau) + hx, mo.y(tau) - ha));
  // The heel starts to lift at 55% of contact, or earlier if the flat-foot knee, once it has passed its deepest bend and is
  // straightening again, would straighten past 30 degrees first. The centre of pressure gets at least 10% of the contact
  // to travel from heel to ball before the heel can lift.
  // (Every step is continuous in the inputs, so the force skew that depends on it is too: the deepest bend is located
  // exactly, not snapped to the search grid.)
  const liftStart = (tc, mo) => {
    const t0 = tFlat * tc, t1 = REAL.rollStart * tc, lim = REAL.liftKnee * D, n = 40, at = j => t0 + (t1 - t0) * j / n;
    const kn = t => kneeFlat(mo, t);
    let jPk = 0; for (let j = 1; j <= n; j++) if (kn(at(j)) > kn(at(jPk))) jPk = j;
    // refine the deepest bend between the neighbouring grid points (golden-section search)
    let a = at(Math.max(0, jPk - 1)), b = at(Math.min(n, jPk + 1));
    const gr = (Math.sqrt(5) - 1) / 2;
    for (let it = 0; it < 60; it++) { const c = b - gr * (b - a), d = a + gr * (b - a); if (kn(c) > kn(d)) b = d; else a = c; }
    const tPk = Math.min(t1, Math.max(t0, (a + b) / 2));
    let tR = t1;
    if (kn(tPk) <= lim) tR = tPk;               // already straighter than 30 degrees at its deepest: lift as it starts to straighten
    else for (let j = jPk; j <= n; j++) if (at(j) > tPk && kn(at(j)) <= lim) {
      let lo = Math.max(tPk, at(j - 1)), hi = at(j);
      for (let it = 0; it < 50; it++) { const m = (lo + hi) / 2; if (kn(m) <= lim) hi = m; else lo = m; }
      tR = hi; break;
    }
    return Math.min(t1, Math.max(tR, t0 + REAL.copMin * tc));
  };
  // centre of pressure: on the heel while the foot lowers, moving to the ball while it is flat (quickly at first, then more
  // slowly), then out along the flat
  // toes to the tip while the heel lifts, and on the tip from then on (where the toes roll, if they do)
  const tipTime = (tc, tR) => tR + REAL.tipAt * (tc - tR);
  const copAt = (tau, tc, tR) => {
    const t0 = tFlat * tc, tT = tipTime(tc, tR);
    return tau < t0 ? heelAt(tau) : tau < tR ? heelAt(tau) + sole * (1 - Math.pow(1 - (tau - t0) / (tR - t0), REAL.copEase))
      : tau < tT ? ballAt(tau) + Lt * sstep((tau - tR) / (tT - tR)) : tipAt(tau);
  };
  // Force skew, body motion and heel-lift start for contact time tc. They depend on each other: for a body motion (set by
  // kappa) and heel-lift start, the horizontal impulse, -integral F * cop / y, is linear in kappa and is zero at `need`.
  // The skew used is the one that is its own need, found by bracketing and root-finding (simple repeated substitution
  // can run away: more skew, earlier heel lift, pressure forward sooner, more skew needed). If no skew within the limit
  // (+-0.5; beyond it the force would go negative) balances, the contact can't be balanced: `need` is then beyond the limit.
  const contact = tc => {
    const k = Math.PI / tc, NS = 240, cap = REAL.skewMax;
    const ev = kap => {
      const mo = motion(tc, kap), tR = liftStart(tc, mo); let i1 = 0, i2 = 0;
      for (let j = 0; j < NS; j++) { const tau = (j + 0.5) * tc / NS, w = copAt(tau, tc, tR) / mo.y(tau); i1 += Math.sin(k * tau) * w; i2 += Math.sin(2 * k * tau) * w; }
      const need = i2 !== 0 ? -i1 / i2 : Infinity;
      return { kap, need, mo, tR, g: need - kap };
    };
    const NK = 8, pts = Array.from({ length: NK + 1 }, (_, j) => ev((2 * j / NK - 1) * cap)), mid = pts[NK / 2];
    // the bracket (sign change of need - kappa) nearest zero skew
    let br = null;
    const order = Array.from({ length: NK }, (_, j) => j).sort((p, q) => Math.abs(p + 0.5 - NK / 2) - Math.abs(q + 0.5 - NK / 2));
    for (const i of order) if (isFinite(pts[i].g) && isFinite(pts[i + 1].g) && pts[i].g * pts[i + 1].g <= 0) { br = [pts[i], pts[i + 1]]; break; }
    if (!br) {   // no balance: use the limit on the side the need points to
      const end = mid.need >= 0 ? pts[NK] : pts[0];
      return { kap: end.kap, need: isFinite(end.need) && Math.abs(end.need) > cap ? end.need : Math.sign(mid.need || 1) * Infinity, mo: end.mo, tR: end.tR };
    }
    let [lo, hi] = br, best = Math.abs(lo.g) < Math.abs(hi.g) ? lo : hi;
    for (let it = 0; it < 40 && Math.abs(best.g) > 1e-12; it++) {   // regula falsi with the Illinois tweak, kept inside the bracket
      let x = hi.kap - hi.g * (hi.kap - lo.kap) / (hi.g - lo.g);
      if (!(x > Math.min(lo.kap, hi.kap) && x < Math.max(lo.kap, hi.kap))) x = (lo.kap + hi.kap) / 2;
      const p = ev(x);
      if (p.g * hi.g < 0) { lo = hi; hi = p; } else { hi = p; lo = { ...lo, g: lo.g / 2 }; }
      best = p;
    }
    // accept only a real root: the bracket can also close on a jump (the heel-lift start can switch abruptly at extreme
    // settings), where need and skew never actually meet
    if (Math.abs(best.g) > 1e-6) return { kap: best.kap, need: Math.sign(best.need || 1) * Infinity, mo: best.mo, tR: best.tR };
    return { kap: best.kap, need: best.kap, mo: best.mo, tR: best.tR };
  };
  const balances = tc => Math.abs(contact(tc).need) <= REAL.skewMax;
  // can the leg still reach the toe tip at push-off with the heel lifted to its working limit, the toes rolled as far as
  // they go and the knee bent at least 10 degrees?
  const blend = 2 * D, rollMax = REAL.rollMax * D, liftCap = rollMax - blend;   // the lift is eased off over its last 2 degrees
  const toeMax = REAL.toeRollMax * D, offCap = (g => { const q = ballOff(liftCap), t = toeOff(g); return [q[0] + t[0], q[1] + t[1]]; })(toeMax);
  const dMin = reach(REAL.kneeMinTO * D);
  const reaches = tc => Math.hypot(tipAt(tc) + offCap[0], contact(tc).mo.y(tc) - offCap[1]) <= dMin;
  // Contact time. Limits: some time in the air (one foot down at a time); the foot must at least pass under the hip; the
  // leg must reach the ground at push-off; and braking and propulsion must be able to cancel (steady speed).
  const tcFloor = () => Math.min(0.95 * Ts, Math.max(0.06, (a + 0.05) / v));
  // If the leg can't reach at push-off even with the shortest contact, the foot lands further forward until it can.
  let landLimited = false;
  if (!reaches(tcFloor())) {
    const setA = x => { a = x; yTD = hipAtTD(a); };
    let lo = a, hi = 0.9 * dTD; setA(hi);
    if (reaches(tcFloor())) {
      for (let it = 0; it < 30; it++) { const mid = (lo + hi) / 2; setA(mid); if (reaches(tcFloor())) hi = mid; else lo = mid; }
      setA(hi);
    }
    landLimited = true;
  }
  const landDeg = Math.atan(a / hipAtTD(a)) / D;
  let reachMax = 0.95 * Ts;
  if (!reaches(reachMax)) {
    let lo = tcFloor(), hi = reachMax;
    for (let it = 0; it < 30; it++) { const mid = (lo + hi) / 2; if (reaches(mid)) lo = mid; else hi = mid; }
    reachMax = lo;
  }
  const reachMin = Math.min(tcFloor(), reachMax);
  // typical: the measured contact time for this speed, scaled with height. A set contact time or contact length is used as
  // set (within the limits); the landing point stays where it is, so the push-off point moves.
  const lenTypical = REAL.contactTime[0] * v ** (REAL.contactTime[1] + 1) * leg.legLen / (0.53 * REAL.refHeight);
  const tcTyp = lenTypical / v;
  const tcWanted = o.contactMs != null ? o.contactMs / 1000 : o.contactLen != null ? o.contactLen / v : tcTyp;
  // stretches of contact time that balance, found on a grid and refined at their edges
  const NG = 40, grid = Array.from({ length: NG + 1 }, (_, j) => reachMin + (reachMax - reachMin) * j / NG), ok = grid.map(balances);
  const edge = (good, bad) => { for (let it = 0; it < 12; it++) { const m = (good + bad) / 2; if (balances(m)) good = m; else bad = m; } return good; };
  const ranges = [];
  for (let j = 0; j <= NG; j++) if (ok[j] && (j === 0 || !ok[j - 1])) {
    let e = j; while (e < NG && ok[e + 1]) e++;
    ranges.push([j > 0 ? edge(grid[j], grid[j - 1]) : grid[j], e < NG ? edge(grid[e], grid[e + 1]) : grid[e]]);
  }
  // use the stretch nearest the wanted contact time; if none balances, the force has to lean off the hip (flagged)
  let tcMin = reachMin, tcMax = reachMax, unbalanced = ranges.length === 0;
  if (!unbalanced) {
    const dist = r => tcWanted < r[0] ? r[0] - tcWanted : tcWanted > r[1] ? tcWanted - r[1] : 0;
    [tcMin, tcMax] = ranges.reduce((b, r) => dist(r) < dist(b) ? r : b);
  }
  // which limit sets each end of the range, for the page: 'pass' (the foot has to get behind the hip), 'min' (60 ms floor),
  // 'balance' (braking and propulsion stop cancelling), 'reach' (the leg can't reach at push-off), 'air' (both feet down)
  const near = (x, y) => Math.abs(x - y) < 1e-6;
  const minWhy = !unbalanced && tcMin > reachMin + 1e-6 ? 'balance' : near(reachMin, 0.06) ? 'min' : 'pass';
  const maxWhy = !unbalanced && tcMax < reachMax - 1e-6 ? 'balance' : reachMax < 0.95 * Ts - 1e-6 ? 'reach' : 'air';
  const tc = Math.max(tcMin, Math.min(tcMax, tcWanted));
  const userSet = o.contactMs != null || o.contactLen != null;
  const tcLimited = userSet && Math.abs(tc - tcWanted) > 5e-4 ? (tc > tcWanted ? 'min' : 'max') : null;
  const tf = Ts - tc, Tsw = T - tc;
  const { kap, mo, tR } = contact(tc), yTO = mo.y(tc), tT = tipTime(tc, tR);
  // knee bend rate from the rate the hip-to-ankle distance d changes: d^2 = L1^2 + Ls^2 + 2 L1 Ls cos(knee)
  const dRateFor = (kneeRate, ph, d) => -kneeRate * L1 * Ls * Math.sin(ph) / d;
  // rate the foot must turn (rad/s) for a given knee rate, when the ankle sits at off(ang) from a ground point gx
  // and moves at offD(ang) per radian of foot turn: d*d' = X*X' + dz*(y' - Z') with X' = -v + w*offD_x, Z' = w*offD_z
  const turnRate = (gx, ang, off, offD, y, vy, kneeRate) => {
    const q = off(ang), qd = offD(ang), X = gx + q[0], dz = y - q[1], d = Math.hypot(X, dz), ph = kneeFor(d);
    return (d * dRateFor(kneeRate, ph, d) + X * v - dz * vy) / (X * qd[0] - dz * qd[1]);
  };
  // touchdown: how fast the foot slaps flat, for a typical knee-bend speed just after touchdown.
  // Toes-up angle b(u) = (1-u)^3 (b0 + (3 b0 - Wtd) u), u = tau / (flatAt tc); Wtd in [0, 4 b0] keeps it monotonic.
  const lowerRate = -turnRate(heelAt(0), b0, heelOff, heelOffD, yTD, mo.vTD, REAL.kneeRateTD * D);   // rad/s, toes coming down
  const Wtd = Math.max(0, Math.min(4 * b0, (isFinite(lowerRate) ? lowerRate : 0) * tFlat * tc));
  // push-off: once the heel starts to lift, the knee may not straighten past a floor curve that runs smoothly from the
  // flat-foot knee at that moment to 22 degrees at toe-off (or to as bent as the leg allows with the heel and toes at their
  // working limits, if it can't reach at 22). The heel lifts just enough to hold the knee on the floor: the ankle sits on a
  // circle of radius rho around the ball and on a circle of radius d around the hip. The lift is smoothly kept between 0 and
  // 65 degrees (no corners, so the joint speeds stay smooth).
  const rho = Math.hypot(ball, ha), del = Math.atan2(ha, ball);
  const ballPos = (tau, g) => { const t = toeOff(g); return [tipAt(tau) + t[0], t[1]]; };
  const liftFor = (tau, d, g) => {   // heel lift (rad) for hip-to-ankle distance d with the toes rolled g (closest the ankle gets if d is too short)
    const B = ballPos(tau, g), y = mo.y(tau), L = Math.hypot(B[0], y - B[1]), c = Math.max(-1, Math.min(1, (rho * rho + L * L - d * d) / (2 * rho * L)));
    return Math.PI - Math.atan2(y - B[1], -B[0]) - Math.acos(c) - del;
  };
  const ramp = x => { if (x <= 0) return 0; if (x >= blend) return x; const s = x / blend; return blend * s * s * s * (6 - 8 * s + 3 * s * s); };  // smooth max(x, 0)
  const hR = 1e-4, kf = t => kneeFlat(mo, t);
  const kR = [kf(tR), (kf(tR + hR) - kf(tR - hR)) / (2 * hR), (kf(tR + hR) - 2 * kf(tR) + kf(tR - hR)) / (hR * hR)];
  // (if the flat-foot knee is still more bent than 22 degrees late in the contact, the floor ends at that bend instead,
  // so the heel still lifts and the knee still arrives at toe-off without extending)
  let kneeTO = Math.max(REAL.kneeTO * D, kf(REAL.holdAt * tc));
  // toe roll at toe-off: none if the heel lift alone lets the leg reach at that knee bend, otherwise as little as does
  // (up to 45 degrees; beyond that the knee straightens instead)
  let toeTO = 0;
  if (liftFor(tc, reach(kneeTO), 0) > liftCap) {
    if (liftFor(tc, reach(kneeTO), toeMax) > liftCap) { toeTO = toeMax; kneeTO = kneeFor(Math.hypot(tipAt(tc) + offCap[0], yTO - offCap[1])); }
    else { let lo = 0, hi = toeMax; for (let it = 0; it < 40; it++) { const m = (lo + hi) / 2; if (liftFor(tc, reach(kneeTO), m) > liftCap) lo = m; else hi = m; } toeTO = hi; }
  }
  const toeAt = tau => tau <= tT ? 0 : toeTO * sstep((tau - tT) / (tc - tT));   // toes roll up about their tip, smoothly
  // Knee speed at toe-off: the one that leaves the thigh momentarily still (hip at its furthest extension, as in real
  // running) while the heel is still rising. With the thigh still, the ankle moves only by the knee bending, at
  // -kneeRate * Ls * (cos psi, sin psi); it must match the ankle's motion from the ground and the heel lift (the toe roll
  // has stopped by then). When the heel is already at its limit that can't be done, and the knee arrives with zero speed.
  const kneeRateTO = (() => {
    const al = liftFor(tc, reach(kneeTO), toeTO), q = ballOff(al), B = ballPos(tc, toeTO), X = B[0] + q[0], dz = yTO - B[1] - q[1], d = Math.hypot(X, dz);
    const ps = Math.atan2(X, dz) + Math.acos(Math.max(-1, Math.min(1, (L1 * L1 + d * d - Ls * Ls) / (2 * L1 * d)))) - kneeFor(d);
    const qd = [ball * Math.sin(al) + ha * Math.cos(al), ball * Math.cos(al) - ha * Math.sin(al)];
    const a1 = -Ls * Math.cos(ps), b1 = -qd[0], e1 = -v, a2 = -Ls * Math.sin(ps), b2 = -qd[1], e2 = -mo.vTO, det = a1 * b2 - a2 * b1;
    const rate = (e1 * b2 - e2 * b1) / det, lift = (a1 * e2 - a2 * e1) / det;
    return al > blend && lift > 0 && rate > 0 && isFinite(rate) ? Math.min(rate, REAL.kneeRateTOMax * D) : 0;
  })();
  const kneeFloor = quintic(kR, [kneeTO, kneeRateTO, 0], tc - tR);
  // Foot angles (anticlockwise +, so toes up +, heel up -) of the rearfoot (ang) and the toes (toeAng), the ankle position,
  // and the heel, ball and toe-tip positions, tau after touchdown. Every piece is smooth across touchdown, lift start,
  // the start of the toe roll and toe-off, so central differences there are accurate.
  function footAt(tau) {
    if (tau < tFlat * tc) {
      const T1 = tFlat * tc, u = tau / T1, B = 3 * b0 - Wtd, b = (1 - u) ** 3 * (b0 + B * u);
      const q = heelOff(b), H = heelAt(tau), c = Math.cos(b), s = Math.sin(b);
      return { X: H + q[0], Z: q[1], ang: b, toeAng: b, heel: [H, 0], ball: [H + sole * c, sole * s], toe: [H + (sole + Lt) * c, (sole + Lt) * s] };
    }
    if (tau < tR) { const H = heelAt(tau); return { X: H + hx, Z: ha, ang: 0, toeAng: 0, heel: [H, 0], ball: [H + sole, 0], toe: [H + sole + Lt, 0] }; }
    const g = toeAt(tau), B = ballPos(tau, g);
    const al = rollMax - ramp(rollMax - ramp(liftFor(tau, reach(evalQ(kneeFloor, tau - tR)[0]), g)));
    const q = ballOff(al);
    return { X: B[0] + q[0], Z: B[1] + q[1], ang: -al, toeAng: -g, heel: [B[0] - sole * Math.cos(al), B[1] + sole * Math.sin(al)], ball: B, toe: [tipAt(tau), 0] };
  }
  const dMax = L1 + Ls;
  function stanceAt(tau) {   // planted-leg geometry tau after touchdown
    const ft = footAt(tau), y = mo.y(tau), dz = y - ft.Z, d = Math.hypot(ft.X, dz);
    const gam = Math.acos(Math.max(-1, Math.min(1, (L1 * L1 + d * d - Ls * Ls) / (2 * L1 * d))));
    return { th: Math.atan2(ft.X, dz) + gam, ph: kneeFor(d), y, short: d - dMax, ...ft };
  }
  const roll = -footAt(tc).ang;
  // least knee bend in the first 30% of contact (a leg that locks straight there can't reach the foot)
  let earlyKneeMin = Infinity;
  for (let j = 1; j <= 40; j++) earlyKneeMin = Math.min(earlyKneeMin, stanceAt(0.3 * tc * j / 40).ph);
  // hip height, vertical velocity and ground force per unit mass at stride time t (t = 0 is this leg's toe-off), either foot down
  function hip(t) {
    const tl = (((t - Tsw) % Ts) + Ts) % Ts;
    if (tl < tc) return [mo.y(tl), mo.vy(tl), mo.f(tl)];
    const u = tl - tc; return [yTO + mo.vTO * u - G * u * u / 2, mo.vTO - G * u, 0];
  }
  let fMax = 0, fAt = 0;
  for (let j = 0; j <= 400; j++) { const s = j / 400, f = Math.sin(Math.PI * s) + kap * Math.sin(2 * Math.PI * s); if (f > fMax) { fMax = f; fAt = s; } }
  return { v, Ts, tc, tf, Tsw, a, ball, sole, landDeg, landTypical, landLimited, tcTypical: tcTyp, lenTypical, tcMin, tcMax, minWhy, maxWhy, tcLimited, unbalanced,
           kap, roll, toeRoll: toeTO, strike: b0, yTD, yTO, vTD: mo.vTD, vTO: mo.vTO,
           liftStart: tR / tc, peakGRF: mo.Fpk * fMax / G, forcePeakAt: fAt, copAt: tau => copAt(tau, tc, tR), stanceAt, hip, kneeTD: kneeTDdeg, earlyKneeMin };
}
// A leg that would lock straight just after landing lands with more knee bend instead: the smallest touchdown knee bend
// (from the typical one up to kneeTDMax) that keeps the early-contact knee above its floor.
function realisticPlan(o, leg, T) {
  const D = Math.PI / 180, ok = P => P.earlyKneeMin >= REAL.kneeFloorTD * D;
  const P0 = planFor(o, leg, T);
  if (ok(P0)) { P0.kneeTDRaised = false; return P0; }
  let hi = REAL.kneeTDMax, Phi = planFor(o, leg, T, hi);
  if (ok(Phi)) {
    let lo = REAL.kneeTD;
    for (let it = 0; it < 8; it++) {
      const mid = (lo + hi) / 2, Pm = planFor(o, leg, T, mid);
      if (ok(Pm)) { hi = mid; Phi = Pm; } else lo = mid;
    }
  }
  Phi.kneeTDRaised = true;
  return Phi;
}

// Body motion and planted-leg samples for one plan at N samples per stride. They don't depend on the swing, so they are
// worked out once and shared by every swing tried (the clearance search) and every knee bend (the sweep).
function stanceSamples(P, T, N) {
  P.samples = P.samples || {};
  if (P.samples[N]) return P.samples[N];
  const { tc, Tsw } = P, h = 5e-5;
  // [value, rate, rate of rate] of the thigh and knee angles in contact, by central differences of the planted-leg geometry
  const derivS = tau => {
    const m = P.stanceAt(tau - h), c = P.stanceAt(tau), p = P.stanceAt(tau + h);
    const d = key => [c[key], (p[key] - m[key]) / (2 * h), (p[key] - 2 * c[key] + m[key]) / (h * h)];
    return { th: d('th'), ph: d('ph'), ang: d('ang'), toeAng: d('toeAng'), g: c };
  };
  // the same at the edges of contact, using only points inside it (dir = -1 looks back from toe-off, +1 forward from
  // touchdown). The geometry beyond the contact is not physical: past toe-off the heel lift runs into its cap, and a
  // difference that straddled it would hand the swing leg speed the stance leg never had.
  const derivEdge = (tau, dir) => {
    const f = [0, 1, 2, 3].map(j => P.stanceAt(tau + dir * j * h)), s = dir / h;
    const d = key => { const [a, b, c, e] = f.map(g => g[key]);
      return [a, s * (-1.5 * a + 2 * b - 0.5 * c), (2 * a - 5 * b + 4 * c - e) / (h * h)]; };
    return { th: d('th'), ph: d('ph'), ang: d('ang'), toeAng: d('toeAng'), g: f[0] };
  };
  const S = blankK(N);
  let kneeLoad = 0, short = -1;
  for (let i = 0; i < N; i++) {
    const t = i * T / N, hp = P.hip(t);
    S.hipY[i] = hp[0]; S.hipVy[i] = hp[1]; S.gEff[i] = hp[2]; S.grf[i] = hp[2] / G;
    if (t < Tsw) continue;
    const tau = t - Tsw, s = derivS(tau), g = s.g;
    S.th[i] = s.th[0]; S.thd[i] = s.th[1]; S.thdd[i] = s.th[2]; S.ph[i] = s.ph[0]; S.phd[i] = s.ph[1]; S.phdd[i] = s.ph[2];
    S.stance[i] = 1; S.footX[i] = g.X; S.ankleZ[i] = g.Z; S.footAng[i] = g.ang; S.footAngRate[i] = s.ang[1]; S.toeAngRate[i] = s.toeAng[1];
    S.heelX[i] = g.heel[0]; S.heelZ[i] = g.heel[1]; S.ballX[i] = g.ball[0]; S.ballZ[i] = g.ball[1]; S.toeX[i] = g.toe[0]; S.toeZ[i] = g.toe[1];
    S.copX[i] = P.copAt(tau);
    kneeLoad = Math.max(kneeLoad, g.ph); short = Math.max(short, g.short);
  }
  let yMax = -1e9, yMin = 1e9; for (let i = 0; i < N; i++) { yMax = Math.max(yMax, S.hipY[i]); yMin = Math.min(yMin, S.hipY[i]); }
  return (P.samples[N] = { S, TO: derivEdge(tc, -1), TD: derivEdge(0, 1), kneeLoad, short, yMax, yMin });
}

function realisticCore(o, leg, T, N, P, clearBoost) {
  const { v, tc, tf, Tsw } = P, { L1, Ls, ha, footLen } = leg, D = Math.PI / 180;
  const SS = stanceSamples(P, T, N), { TO, TD, kneeLoad, short, yMax, yMin } = SS;   // TO, TD: states at toe-off and touchdown
  // Swing curves match the angle and its rate at toe-off and touchdown. They need not match the angular acceleration:
  // the stance leg is treated as massless, so its torques jump at those instants anyway, and copying the large
  // accelerations of the planted-leg geometry would only make the swing leg lurch.
  // thigh in swing: from toe-off forward to its peak flexion, then back (retracting) to touchdown. Peak flexion grows
  // with speed as measured in running (from vertical, 154 trials at 3.1-10 m/s; Clark, Meng & Stearne 2020) and is
  // always at least 4 degrees past the touchdown angle; clearBoost adds extra swing when the foot would scuff on landing
  const thMax = Math.max(REAL.thighFlex[0] * v + REAL.thighFlex[1], TD.th[0] + REAL.retractMin * D) + clearBoost * D, sTh = 0.72;
  const thighSwing = swingCurve([0, sTh * Tsw, Tsw], [TO.th[0], thMax, TD.th[0]], TO.th[1], TD.th[1], 1);
  // knee in swing: from toe-off, through the peak bend and a terminal-swing bend, into touchdown
  const phTypical = Math.min(130, 60 + 12 * v);
  const phMax = o.realKnee != null ? o.realKnee : phTypical;
  // The slope it must end with is capped: beyond kneeRateTDSwingMax the knee would have to bend backwards late in swing to
  // arrive that fast. The remainder is taken up as a jump in knee speed at touchdown (an impact), which costs nothing here
  // because the planted leg is treated as weightless.
  const kneeSwing = swingCurve([0, 0.42 * Tsw, 0.88 * Tsw, Tsw], [TO.ph[0], phMax * D, REAL.termKnee * D, TD.ph[0]], TO.ph[1], Math.min(TD.ph[1], REAL.kneeRateTDSwingMax * D), 1);
  const K = blankK(N);
  for (const key in SS.S) K[key].set(SS.S[key]);   // body motion and the planted leg
  for (let i = 0; i < N; i++) {
    const t = i * T / N;
    if (t >= Tsw) continue;
    const th = thighSwing(t), ph = kneeSwing(t);
    K.th[i] = th[0]; K.thd[i] = th[1]; K.thdd[i] = th[2]; K.ph[i] = ph[0]; K.phd[i] = ph[1]; K.phdd[i] = ph[2]; K.count[i] = 1;
    // foot pitch (toe-down +): picks up from the heel-up angle at toe-off, square to the shin through mid-swing,
    // then toes up for the heel-first landing
    const w0 = sstep(t / (0.2 * Tsw)), wgt = sstep((t / Tsw - 0.45) / 0.45);
    K.footPitch[i] = (1 - w0) * P.roll + w0 * ((1 - wgt) * (-(th[0] - ph[0])) - wgt * P.strike);
  }
  // foot clearance: lowest point of the swinging foot (heel or ball; toes can bend) from early swing until just before
  // touchdown. lateClear covers only the approach to landing, which the wrapper fixes by swinging the thigh further.
  let minClear = 1e9, lateClear = 1e9;
  for (let i = 0; i < N; i++) {
    const t = i * T / N; if (t < 0.12 * Tsw || t > 0.93 * Tsw) continue;
    const ps = K.th[i] - K.ph[i];
    const p = K.footPitch[i], ankleH = K.hipY[i] - L1 * Math.cos(K.th[i]) - Ls * Math.cos(ps);
    const ballH = ankleH - (0.45 * footLen * Math.sin(p) + ha * Math.cos(p));
    const heelH = ankleH - (-0.25 * footLen * Math.sin(p) + ha * Math.cos(p));
    const c = Math.min(ballH, heelH);
    minClear = Math.min(minClear, c); if (t >= 0.65 * Tsw) lateClear = Math.min(lateClear, c);
  }
  K.info = { tc, tf, Tsw, stepLen: v * P.Ts, peakGRF: P.peakGRF, forcePeakAt: P.forcePeakAt, skew: P.kap, rollDeg: P.roll / D, strikeDeg: P.strike / D,
             kneeLoadDeg: kneeLoad / D, kneeTODeg: TO.ph[0] / D, kneeRateTD: TD.ph[1] / D, kneeRateTO: TO.ph[1] / D, reachShort: short, kneeTDDeg: P.kneeTD, kneeTDRaised: P.kneeTDRaised,
             bob: yMax - yMin, hipMax: yMax, phMaxDeg: phMax, phTypical, minClear, lateClear, clearBoost, thighPeakDeg: thMax / D, a: P.a, ball: P.ball,
             landDeg: P.landDeg, landTypical: P.landTypical, landLimited: P.landLimited, contactLen: v * tc, tcTypical: P.tcTypical, lenTypical: P.lenTypical, toeRollDeg: P.toeRoll / D, tcMin: P.tcMin, tcMax: P.tcMax, minWhy: P.minWhy, maxWhy: P.maxWhy,
             tcLimited: P.tcLimited, unbalanced: P.unbalanced, pushBehind: v * tc - P.a };
  return K;
}
// Smooth curve through knots (xs, ys) with set slopes at both ends: a cubic spline, except that the knot `flat` gets zero
// slope so its value is an exact turning point (quintic pieces carry the spline's curvature there). Returns t -> [y, y', y''].
function swingCurve(xs, ys, d0, dn, flat) {
  const sp = clampedSpline(xs, ys, d0, dn);
  const st = xs.map((x, i) => [ys[i], i === flat ? 0 : sp.d[i], sp.dd[i]]);
  const segs = xs.slice(1).map((x, i) => quintic(st[i], st[i + 1], x - xs[i]));
  return t => { let j = 0; while (j < segs.length - 1 && t >= xs[j + 1]) j++; return evalQ(segs[j], t - xs[j]); };
}
// cubic spline through (xs, ys) with set end slopes d0 and dn; returns the slope and second derivative at each knot
function clampedSpline(xs, ys, d0, dn) {
  const n = xs.length, hh = [], A = [], b = [];
  for (let i = 0; i < n - 1; i++) hh.push(xs[i + 1] - xs[i]);
  for (let i = 0; i < n; i++) {
    const row = new Array(n).fill(0);
    if (i === 0) { row[0] = 2 * hh[0]; row[1] = hh[0]; b.push(6 * ((ys[1] - ys[0]) / hh[0] - d0)); }
    else if (i === n - 1) { row[n - 2] = hh[n - 2]; row[n - 1] = 2 * hh[n - 2]; b.push(6 * (dn - (ys[n - 1] - ys[n - 2]) / hh[n - 2])); }
    else { row[i - 1] = hh[i - 1]; row[i] = 2 * (hh[i - 1] + hh[i]); row[i + 1] = hh[i]; b.push(6 * ((ys[i + 1] - ys[i]) / hh[i] - (ys[i] - ys[i - 1]) / hh[i - 1])); }
    A.push(row);
  }
  const dd = solve(A, b);
  const d = xs.map((x, i) => i < n - 1 ? (ys[i + 1] - ys[i]) / hh[i] - hh[i] * (2 * dd[i] + dd[i + 1]) / 6 : dn);
  return { d, dd };
}
function solve(A, b) {
  const n = b.length, M = A.map((row, i) => [...row, b[i]]);
  for (let c = 0; c < n; c++) {
    let piv = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    for (let r = 0; r < n; r++) if (r !== c) { const f = M[r][c] / M[c][c]; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
  }
  return M.map((row, i) => row[n] / row[i]);
}
// If the foot would go more than DRAG_TOL into the ground on its way in to landing (the same allowance the drag warning
// uses, since the foot shape is simplified), the thigh swings a little further forward before pulling back to land, which
// lifts the foot. Find (by bisection) the smallest extra forward swing, up to 15 degrees, that keeps it within that.
// (Bending the knee more late in swing is not a reliable fix: with the thigh forward, a bend near the thigh's own angle
// leaves the shin vertical and the foot at its lowest, and unbending it again before touchdown takes huge torques.)
let planCache = null;   // the plan doesn't depend on the swing knee bend, so the knee-bend sweep reuses it
function realisticKinematics(o, leg, T, N) {
  const key = [o.sex, o.H, o.M, o.speed, o.cadence, o.landDeg, o.contactMs, o.contactLen].join('|');
  if (!planCache || planCache.key !== key) planCache = { key, P: realisticPlan(o, leg, T) };
  const P = planCache.P;
  const core = boost => realisticCore(o, leg, T, N, P, boost);
  const K = core(0);
  if (K.info.lateClear >= -DRAG_TOL) return K;
  const top = core(REAL.clearBoostMax);
  if (top.info.lateClear < -DRAG_TOL) return top;
  let lo = 0, hi = REAL.clearBoostMax, best = top;
  for (let it = 0; it < 6; it++) {
    const mid = (lo + hi) / 2, Km = core(mid);
    if (Km.info.lateClear >= -DRAG_TOL) { hi = mid; best = Km; } else lo = mid;
  }
  return best;
}
// quintic Hermite: matches position, velocity, acceleration at both ends over duration D
function quintic(p0, p1, D) {
  const [x0, v0, a0] = p0, [x1, v1, a1] = p1;
  const c0 = x0, c1 = v0 * D, c2 = a0 * D * D / 2;
  const X = x1 - c0 - c1 - c2, V = v1 * D - c1 - 2 * c2, A = a1 * D * D - 2 * c2;
  const c3 = 10 * X - 4 * V + A / 2, c4 = -15 * X + 7 * V - A, c5 = 6 * X - 3 * V + A / 2;
  return { c: [c0, c1, c2, c3, c4, c5], D };
}
function evalQ(q, t) {
  const s = t / q.D, c = q.c;
  const p = c[0] + s * (c[1] + s * (c[2] + s * (c[3] + s * (c[4] + s * c[5]))));
  const d = c[1] + s * (2 * c[2] + s * (3 * c[3] + s * (4 * c[4] + s * 5 * c[5])));
  const dd = 2 * c[2] + s * (6 * c[3] + s * (12 * c[4] + s * 20 * c[5]));
  return [p, d / q.D, dd / (q.D * q.D)];
}

function simulate(o) {
  // o: sex, H (m), M (kg), mode, phmaxDeg, cadence, arcDeg, gravity, share, elastic(0-1), eff(0-1), speed (m/s, realistic mode)
  const leg = legModel(o.sex, o.H, o.M);
  const { L1, m1, c1, I1, m2, c2, I2 } = leg;
  const f = o.cadence / 120;            // strides per second, per leg
  const T = 1 / f, N = o.N || 1440, dt = T / N;
  const K = o.mode === 'real' ? realisticKinematics(o, leg, T, N) : simpleKinematics(o, T, N);
  const M11 = I1 + m1 * c1 * c1 + m2 * L1 * L1, M22 = I2 + m2 * c2 * c2;
  const S = { th: K.th, ph: K.ph, tauH: new Float64Array(N), tauK: new Float64Array(N), Ph: new Float64Array(N), Pk: new Float64Array(N), K };
  // hip and knee torques that move the leg's own mass along its path, with the hip's vertical acceleration felt as extra gravity
  const legTorques = i => {
    const th = K.th[i], thd = K.thd[i], thdd = K.thdd[i], g = K.gEff[i];
    const ps = th - K.ph[i], psd = thd - K.phd[i], psdd = thdd - K.phdd[i];
    const d = K.ph[i];
    const M12 = m2 * L1 * c2 * Math.cos(d), sd = Math.sin(d);
    const Q1 = M11 * thdd + M12 * psdd + m2 * L1 * c2 * sd * psd * psd + (m1 * c1 + m2 * L1) * g * Math.sin(th);
    const Q2 = M12 * thdd + M22 * psdd - m2 * L1 * c2 * sd * thd * thd + m2 * c2 * g * Math.sin(ps);
    return [Q1 + Q2, Q2];   // [hip, knee]
  };
  let hp = 0, hn = 0, kp = 0, kn = 0, sp = 0, sn = 0, peakKnee = 0, peakHip = 0;
  for (let i = 0; i < N; i++) {
    if (!K.count[i]) { S.tauH[i] = NaN; S.tauK[i] = NaN; S.Ph[i] = NaN; S.Pk[i] = NaN; continue; }
    const thd = K.thd[i], psd = thd - K.phd[i];
    const [tauH, tauK] = legTorques(i);
    const Ph = tauH * thd, Pk = tauK * (psd - thd);
    S.tauH[i] = tauH; S.tauK[i] = tauK; S.Ph[i] = Ph; S.Pk[i] = Pk;
    if (Ph > 0) hp += Ph * dt; else hn -= Ph * dt;
    if (Pk > 0) kp += Pk * dt; else kn -= Pk * dt;
    const Ps = Ph + Pk; if (Ps > 0) sp += Ps * dt; else sn -= Ps * dt;
    peakKnee = Math.max(peakKnee, Math.abs(K.phd[i])); peakHip = Math.max(peakHip, Math.abs(thd));
  }
  // ---------- stance (realistic mode): ground force at the centre of pressure, leg treated as massless ----------
  let stance = null;
  if (o.mode === 'real') {
    S.Pa = new Float64Array(N).fill(NaN); S.tauA = new Float64Array(N).fill(NaN); S.Pm = new Float64Array(N).fill(NaN); S.tauM = new Float64Array(N).fill(NaN);
    S.Fx = new Float64Array(N); S.Fy = new Float64Array(N);
    // horizontal force: along the line from the centre of pressure toward the hip. The force timing (skew) makes braking
    // and propulsion cancel. Only when no contact time can balance (flagged on the page) does c tilt the line by what's left.
    // (Not otherwise: the samples would leave a tiny imbalance, and at very short contacts a tiny tilt makes real hip work.)
    let num = 0, den = 0;
    for (let i = 0; i < N; i++) if (K.stance[i]) { const Fy = K.grf[i]; num += Fy * (-K.copX[i] / K.hipY[i]); den += Fy; }
    const c = K.info.unbalanced && den > 0 ? num / den : 0;
    let hp2 = 0, hn2 = 0, kp2 = 0, kn2 = 0, ap = 0, an = 0, mp = 0, mn = 0, sp2 = 0, sn2 = 0, ext = 0, sumCheck = 0;
    for (let i = 0; i < N; i++) {
      if (!K.stance[i]) continue;
      const Fy = o.M * G * K.grf[i], Fx = Fy * (-K.copX[i] / K.hipY[i] - c);
      S.Fx[i] = Fx / (o.M * G); S.Fy[i] = K.grf[i];
      const th = K.th[i], hipY = K.hipY[i];
      const kx = L1 * Math.sin(th), ky = hipY - L1 * Math.cos(th);
      const axp = K.footX[i], ayp = K.ankleZ[i], cx = K.copX[i];
      const cross = (rx, ry) => -(rx * Fy - ry * Fx);        // moment on the distal part from the proximal one: -(r_cop - r_joint) x F
      const Mh = cross(cx - 0, 0 - hipY), Mk = cross(cx - kx, 0 - ky), Ma = cross(cx - axp, 0 - ayp);
      // toe joint (ball of the foot): moment on the toes from the rearfoot. The toes only carry load once the pressure is
      // on them; while it is on the heel or the rearfoot they are unloaded.
      const Mm = cx > K.ballX[i] + 1e-9 ? cross(cx - K.ballX[i], 0 - K.ballZ[i]) : 0;
      const wT = K.thd[i], wS = K.thd[i] - K.phd[i], wF = K.footAngRate[i], wToe = K.toeAngRate[i];
      const Ph = Mh * wT, Pk = Mk * (wS - wT), Pa = Ma * (wF - wS), Pm = Mm * (wToe - wF);
      S.tauH[i] = Mh; S.tauK[i] = Mk; S.tauA[i] = Ma; S.tauM[i] = Mm; S.Ph[i] = Ph; S.Pk[i] = Pk; S.Pa[i] = Pa; S.Pm[i] = Pm;
      if (Ph > 0) hp2 += Ph * dt; else hn2 -= Ph * dt;
      if (Pk > 0) kp2 += Pk * dt; else kn2 -= Pk * dt;
      if (Pa > 0) ap += Pa * dt; else an -= Pa * dt;
      if (Pm > 0) mp += Pm * dt; else mn -= Pm * dt;
      const Ps = Ph + Pk + Pa + Pm; if (Ps > 0) sp2 += Ps * dt; else sn2 -= Ps * dt;
      // the foot point under the centre of pressure never moves, so joint power must equal the ground force's power on the hip
      const Pext = Fx * o.speed + Fy * K.hipVy[i];
      sumCheck += Math.abs(Ps - Pext) * dt; ext += Math.abs(Pext) * dt;
    }
    const r0 = o.elastic;
    const perStep = o.share ? Math.max(sp2 - r0 * sn2, 0)
                            : Math.max(hp2 - r0 * hn2, 0) + Math.max(kp2 - r0 * kn2, 0) + Math.max(ap - r0 * an, 0) + Math.max(mp - r0 * mn, 0);
    const mechS = perStep * f * 2;
    stance = { hipPos: hp2, hipNeg: hn2, kneePos: kp2, kneeNeg: kn2, anklePos: ap, ankleNeg: an, toePos: mp, toeNeg: mn, pos: sp2, neg: sn2, net: sp2 - sn2, perStep,
               mech: mechS, metab: mechS / o.eff, balanceErr: ext > 0 ? sumCheck / ext : 0, tilt: c };
  }
  const r = o.elastic;
  const perStride = o.share ? Math.max(sp - r * sn, 0)
                            : Math.max(hp - r * hn, 0) + Math.max(kp - r * kn, 0);
  const mech = perStride * f * 2;        // both legs, W
  const metab = mech / o.eff;
  return { stance, leg, S, K, info: K.info || null, T, f, N, dt, hipPos: hp, kneePos: kp, hipNeg: hn, kneeNeg: kn, perStride, mech, metab,
           kcalh: metab * 3600 / 4184, peakKneeDeg: peakKnee * 180 / Math.PI, peakHipDeg: peakHip * 180 / Math.PI };
}
