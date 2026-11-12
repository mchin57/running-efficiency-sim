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
           footX: z(), ankleZ: z(), footAng: z(), footAngRate: z(), heelX: z(), heelZ: z(), ballX: z(), ballZ: z(), toesDown: z() };
}

// ---------- realistic stride: real ground, planted foot that rolls heel to ball, bouncing body ----------
// The body (lumped at the hip) is driven by the ground force. Per unit body mass it is
//   F/m = Fpk [sin(pi s) + kappa sin(2 pi s)]   during contact (s = fraction of the contact), 0 in the air.
// The kappa term moves the force earlier or later in the contact without changing its impulse. kappa is chosen so that,
// with the force pointing from the centre of pressure at the hip, braking and propulsion cancel (steady speed).
// The hip height follows from that force in closed form, so the swing leg, the stance leg and the drawing all share one
// body motion. During contact the foot is planted and rigid: it lands on the heel with the toes up and lowers flat, stays
// flat, then rolls over the ball as the heel lifts. The knee and thigh angles are whatever geometry makes the leg reach
// from the hip down to the ankle. The foot lowers at the rate that makes the knee bend at a typical speed just after
// touchdown; later the heel lifts just enough to stop the knee straightening past a floor that ends near 22 degrees at toe-off.
// Typical landing: ankle touches down 0.256 leg lengths ahead of the hip and leaves 0.544 leg lengths behind it (0.8 leg lengths of contact).
const REAL = { landAhead: 0.256, pushBehind: 0.544, kneeTD: 20, kneeTO: 22, kneeMinTO: 10, kneeRateTD: 300, kneeRateTOMax: 600,
               strikeDeg: 12, flatAt: 0.15, rollStart: 0.45, liftKnee: 30, holdAt: 0.85, rollMax: 65, skewMax: 0.45 };
const sstep = s => { s = Math.min(1, Math.max(0, s)); return s * s * s * (s * (6 * s - 15) + 10); };

// Everything about the stride that does not depend on the late-swing knee bend: contact timing, ground force,
// body motion and the stance-leg geometry.
function realisticPlan(o, leg, T) {
  const v = o.speed, Ts = T / 2, D = Math.PI / 180;
  const { L1, Ls, ha, footLen } = leg;
  const ball = 0.45 * footLen, hx = 0.25 * footLen, sole = hx + ball;   // ankle-to-ball, heel-to-ankle, heel-to-ball (along the sole)
  const reach = ph => Math.sqrt(L1 * L1 + Ls * Ls + 2 * L1 * Ls * Math.cos(ph));      // hip-to-ankle distance at knee bend ph
  const kneeFor = d => Math.acos(Math.max(-1, Math.min(1, (d * d - L1 * L1 - Ls * Ls) / (2 * L1 * Ls))));
  // ankle position relative to the heel when the toes are up by b, and relative to the ball when the heel is up by al
  // (rigid foot turning about the point on the ground), with their rates of change per radian
  const heelOff = b => [hx * Math.cos(b) - ha * Math.sin(b), hx * Math.sin(b) + ha * Math.cos(b)];
  const heelOffD = b => [-hx * Math.sin(b) - ha * Math.cos(b), hx * Math.cos(b) - ha * Math.sin(b)];
  const ballOff = al => [-ball * Math.cos(al) + ha * Math.sin(al), ball * Math.sin(al) + ha * Math.cos(al)];
  // landing point: angle from vertical of the line from the hip to the ground under the ankle at touchdown
  const b0 = REAL.strikeDeg * D, strike = heelOff(b0);
  const dTD = reach(REAL.kneeTD * D);
  const hipAtTD = x => Math.sqrt(Math.max(dTD * dTD - x * x, 1e-6)) + strike[1];
  const aTyp = REAL.landAhead * leg.legLen;
  const landTypical = Math.atan(aTyp / hipAtTD(aTyp)) / D;
  const landDeg = o.landDeg != null ? o.landDeg : landTypical;
  let a = aTyp;
  for (let it = 0; it < 60; it++) a = Math.min(0.9 * dTD, Math.tan(landDeg * D) * hipAtTD(a));
  const yTD = hipAtTD(a);                                   // hip height at touchdown: heel down, toes up, knee at 20 degrees
  // ground points measured from the point on the ground under the hip, tau after touchdown
  const heelAt = tau => a - strike[0] - v * tau, ballAt = tau => heelAt(tau) + sole;
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
  // the heel starts to lift at 45% of contact, or earlier if the flat-foot knee would otherwise straighten past 30 degrees
  const liftStart = (tc, mo) => {
    const t0 = tFlat * tc, t1 = REAL.rollStart * tc, lim = REAL.liftKnee * D, below = t => kneeFlat(mo, t) <= lim;
    if (below(t0)) return t0;
    for (let j = 1, prev = t0; j <= 40; prev = t0 + (t1 - t0) * j / 40, j++) {
      const t = t0 + (t1 - t0) * j / 40;
      if (below(t)) { let lo = prev, hi = t; for (let it = 0; it < 40; it++) { const m = (lo + hi) / 2; if (below(m)) hi = m; else lo = m; } return hi; }
    }
    return t1;
  };
  // centre of pressure: on the heel while the foot lowers, moving to the ball while it is flat, on the ball while the heel lifts
  const copAt = (tau, tc, tR) => tau < tFlat * tc ? heelAt(tau) : tau < tR ? heelAt(tau) + sole * (tau - tFlat * tc) / (tR - tFlat * tc) : ballAt(tau);
  // Force skew, body motion and heel-lift start for contact time tc. They depend on each other, so iterate. For a given body
  // motion the horizontal impulse, -integral F * cop / y, is linear in kappa.
  const contact = tc => {
    const k = Math.PI / tc, NS = 240; let kap = 0, mo = motion(tc, 0), tR = liftStart(tc, mo);
    for (let it = 0; it < 6; it++) {
      let i1 = 0, i2 = 0;
      for (let j = 0; j < NS; j++) { const tau = (j + 0.5) * tc / NS, w = copAt(tau, tc, tR) / mo.y(tau); i1 += Math.sin(k * tau) * w; i2 += Math.sin(2 * k * tau) * w; }
      kap = Math.max(-REAL.skewMax, Math.min(REAL.skewMax, i2 !== 0 ? -i1 / i2 : 0));
      mo = motion(tc, kap); tR = liftStart(tc, mo);
    }
    return { kap, mo, tR };
  };
  // can the leg still reach the ball at push-off with the heel lifted to its working limit and the knee bent at least 10 degrees?
  const blend = 2 * D, rollMax = REAL.rollMax * D, liftCap = rollMax - blend;   // the lift is eased off over its last 2 degrees
  const offCap = ballOff(liftCap), dMin = reach(REAL.kneeMinTO * D);
  const reaches = tc => Math.hypot(ballAt(tc) + offCap[0], contact(tc).mo.y(tc) - offCap[1]) <= dMin;
  // contact time: typical = foot leaves the usual distance behind the hip. A set time moves the push-off point instead.
  // Limits: some time in the air (one foot down at a time), the foot must at least pass under the hip, and the leg must reach.
  let tcMax = 0.95 * Ts;
  if (!reaches(tcMax)) {
    let lo = Math.min(tcMax, Math.max(0.06, (a + 0.05) / v)), hi = tcMax;
    for (let it = 0; it < 30; it++) { const mid = (lo + hi) / 2; if (reaches(mid)) lo = mid; else hi = mid; }
    tcMax = lo;
  }
  const tcMin = Math.min(Math.max(0.06, (a + 0.05) / v), tcMax);
  const tcTypical = Math.max(tcMin, Math.min((a + REAL.pushBehind * leg.legLen) / v, 0.9 * Ts, tcMax));
  const tcWanted = o.contactMs != null ? o.contactMs / 1000 : tcTypical;
  const tc = Math.max(tcMin, Math.min(tcMax, tcWanted));
  const tcLimited = Math.abs(tc - tcWanted) > 5e-4 ? (tc > tcWanted ? 'min' : 'max') : null;
  const tf = Ts - tc, Tsw = T - tc;
  const { kap, mo, tR } = contact(tc), yTO = mo.y(tc);
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
  // flat-foot knee at that moment to 22 degrees at toe-off (or to as bent as the leg allows with the heel at its working
  // limit, if it can't reach at 22). The heel lifts just enough to hold the knee on the floor: the ankle sits on a circle
  // of radius rho around the ball and on a circle of radius d around the hip. The lift is smoothly kept between 0 and
  // 65 degrees (no corners, so the joint speeds stay smooth).
  const rho = Math.hypot(ball, ha), del = Math.atan2(ha, ball);
  const liftFor = (tau, d) => {      // heel lift (rad) for hip-to-ankle distance d (the closest the ankle can get if d is too short)
    const Bx = ballAt(tau), y = mo.y(tau), L = Math.hypot(Bx, y), c = Math.max(-1, Math.min(1, (rho * rho + L * L - d * d) / (2 * rho * L)));
    return Math.PI - Math.atan2(y, -Bx) - Math.acos(c) - del;
  };
  const ramp = x => { if (x <= 0) return 0; if (x >= blend) return x; const s = x / blend; return blend * s * s * s * (6 - 8 * s + 3 * s * s); };  // smooth max(x, 0)
  const hR = 1e-4, kf = t => kneeFlat(mo, t);
  const kR = [kf(tR), (kf(tR + hR) - kf(tR - hR)) / (2 * hR), (kf(tR + hR) - 2 * kf(tR) + kf(tR - hR)) / (hR * hR)];
  // (if the flat-foot knee is still more bent than 22 degrees late in the contact, the floor ends at that bend instead,
  // so the heel still lifts and the knee still arrives at toe-off without extending)
  let kneeTO = Math.max(REAL.kneeTO * D, kf(REAL.holdAt * tc));
  if (liftFor(tc, reach(kneeTO)) > liftCap) kneeTO = kneeFor(Math.hypot(ballAt(tc) + offCap[0], yTO - offCap[1]));
  // Knee speed at toe-off: the one that leaves the thigh momentarily still (hip at its furthest extension, as in real
  // running) while the heel is still rising. With the thigh still, the ankle moves only by the knee bending, at
  // -kneeRate * Ls * (cos psi, sin psi); it must match the ankle's motion from the ground and the heel lift. When the heel
  // is already at its limit that can't be done, and the knee arrives at toe-off with zero speed instead.
  const kneeRateTO = (() => {
    const al = liftFor(tc, reach(kneeTO)), q = ballOff(al), X = ballAt(tc) + q[0], dz = yTO - q[1], d = Math.hypot(X, dz);
    const ps = Math.atan2(X, dz) + Math.acos(Math.max(-1, Math.min(1, (L1 * L1 + d * d - Ls * Ls) / (2 * L1 * d)))) - kneeFor(d);
    const qd = [ball * Math.sin(al) + ha * Math.cos(al), ball * Math.cos(al) - ha * Math.sin(al)];
    const a1 = -Ls * Math.cos(ps), b1 = -qd[0], e1 = -v, a2 = -Ls * Math.sin(ps), b2 = -qd[1], e2 = -mo.vTO, det = a1 * b2 - a2 * b1;
    const rate = (e1 * b2 - e2 * b1) / det, lift = (a1 * e2 - a2 * e1) / det;
    return al > blend && lift > 0 && rate > 0 && isFinite(rate) ? Math.min(rate, REAL.kneeRateTOMax * D) : 0;
  })();
  const kneeFloor = quintic(kR, [kneeTO, kneeRateTO, 0], tc - tR);
  // foot angle (anticlockwise +, so toes up +, heel up -) and the ankle position, tau after touchdown. Every piece is smooth
  // across touchdown, lift start and toe-off, so central differences there are accurate.
  function footAt(tau) {
    if (tau < tFlat * tc) {
      const T1 = tFlat * tc, u = tau / T1, B = 3 * b0 - Wtd, b = (1 - u) ** 3 * (b0 + B * u);
      const q = heelOff(b), H = heelAt(tau);
      return { X: H + q[0], Z: q[1], ang: b, heel: [H, 0], ball: [H + sole * Math.cos(b), sole * Math.sin(b)], toesDown: false };
    }
    if (tau < tR) { const H = heelAt(tau); return { X: H + hx, Z: ha, ang: 0, heel: [H, 0], ball: [H + sole, 0], toesDown: true }; }
    const al = rollMax - ramp(rollMax - ramp(liftFor(tau, reach(evalQ(kneeFloor, tau - tR)[0]))));
    const q = ballOff(al), Bx = ballAt(tau);
    return { X: Bx + q[0], Z: q[1], ang: -al, heel: [Bx - sole * Math.cos(al), sole * Math.sin(al)], ball: [Bx, 0], toesDown: true };
  }
  const dMax = L1 + Ls;
  function stanceAt(tau) {   // planted-leg geometry tau after touchdown
    const ft = footAt(tau), y = mo.y(tau), dz = y - ft.Z, d = Math.hypot(ft.X, dz);
    const gam = Math.acos(Math.max(-1, Math.min(1, (L1 * L1 + d * d - Ls * Ls) / (2 * L1 * d))));
    return { th: Math.atan2(ft.X, dz) + gam, ph: kneeFor(d), y, short: d - dMax, ...ft };
  }
  const roll = -footAt(tc).ang;
  // hip height, vertical velocity and ground force per unit mass at stride time t (t = 0 is this leg's toe-off), either foot down
  function hip(t) {
    const tl = (((t - Tsw) % Ts) + Ts) % Ts;
    if (tl < tc) return [mo.y(tl), mo.vy(tl), mo.f(tl)];
    const u = tl - tc; return [yTO + mo.vTO * u - G * u * u / 2, mo.vTO - G * u, 0];
  }
  let fMax = 0, fAt = 0;
  for (let j = 0; j <= 400; j++) { const s = j / 400, f = Math.sin(Math.PI * s) + kap * Math.sin(2 * Math.PI * s); if (f > fMax) { fMax = f; fAt = s; } }
  return { v, Ts, tc, tf, Tsw, a, ball, sole, landDeg, landTypical, tcTypical, tcMin, tcMax, tcLimited, kap, roll, strike: b0, yTD, yTO, vTD: mo.vTD, vTO: mo.vTO,
           liftStart: tR / tc, peakGRF: mo.Fpk * fMax / G, forcePeakAt: fAt, copAt: tau => copAt(tau, tc, tR), stanceAt, hip };
}

function realisticCore(o, leg, T, N, P, termKnee) {
  const { v, tc, tf, Tsw } = P, { L1, Ls, ha, footLen } = leg, D = Math.PI / 180;
  const h = 1e-4;
  // [value, rate, rate of rate] of the thigh and knee angles in contact, by central differences of the planted-leg geometry
  const derivS = tau => {
    const m = P.stanceAt(tau - h), c = P.stanceAt(tau), p = P.stanceAt(tau + h);
    const d = key => [c[key], (p[key] - m[key]) / (2 * h), (p[key] - 2 * c[key] + m[key]) / (h * h)];
    return { th: d('th'), ph: d('ph'), ang: d('ang'), g: c };
  };
  const TO = derivS(tc), TD = derivS(0);   // states at toe-off and touchdown (the geometry is smooth across both)
  // Swing curves match the angle and its rate at toe-off and touchdown. They need not match the angular acceleration:
  // the stance leg is treated as massless, so its torques jump at those instants anyway, and copying the large
  // accelerations of the planted-leg geometry would only make the swing leg lurch.
  // thigh in swing: from toe-off forward to its furthest point, then back slightly to touchdown
  const thMax = TD.th[0] + (4 + 1.5 * v) * D, sTh = 0.72;
  const thighSwing = swingCurve([0, sTh * Tsw, Tsw], [TO.th[0], thMax, TD.th[0]], TO.th[1], TD.th[1], 1);
  // knee in swing: from toe-off, through the peak bend and a terminal-swing bend, into touchdown
  const phTypical = Math.min(130, 60 + 12 * v);
  const phMax = o.realKnee != null ? o.realKnee : phTypical;
  const kneeSwing = swingCurve([0, 0.42 * Tsw, 0.88 * Tsw, Tsw], [TO.ph[0], phMax * D, termKnee * D, TD.ph[0]], TO.ph[1], TD.ph[1], 1);
  const K = blankK(N);
  let kneeLoad = 0, short = -1;
  for (let i = 0; i < N; i++) {
    const t = i * T / N, hp = P.hip(t);
    K.hipY[i] = hp[0]; K.hipVy[i] = hp[1]; K.gEff[i] = hp[2]; K.grf[i] = hp[2] / G;
    if (t < Tsw) {
      const th = thighSwing(t), ph = kneeSwing(t);
      K.th[i] = th[0]; K.thd[i] = th[1]; K.thdd[i] = th[2]; K.ph[i] = ph[0]; K.phd[i] = ph[1]; K.phdd[i] = ph[2]; K.count[i] = 1;
      // foot pitch (toe-down +): picks up from the heel-up angle at toe-off, square to the shin through mid-swing,
      // then toes up for the heel-first landing
      const w0 = sstep(t / (0.2 * Tsw)), wgt = sstep((t / Tsw - 0.45) / 0.45);
      K.footPitch[i] = (1 - w0) * P.roll + w0 * ((1 - wgt) * (-(th[0] - ph[0])) - wgt * P.strike);
    } else {
      const tau = t - Tsw, s = derivS(tau), g = s.g;
      K.th[i] = s.th[0]; K.thd[i] = s.th[1]; K.thdd[i] = s.th[2]; K.ph[i] = s.ph[0]; K.phd[i] = s.ph[1]; K.phdd[i] = s.ph[2];
      K.stance[i] = 1; K.footX[i] = g.X; K.ankleZ[i] = g.Z; K.footAng[i] = g.ang; K.footAngRate[i] = s.ang[1];
      K.heelX[i] = g.heel[0]; K.heelZ[i] = g.heel[1]; K.ballX[i] = g.ball[0]; K.ballZ[i] = g.ball[1]; K.toesDown[i] = g.toesDown ? 1 : 0;
      K.copX[i] = P.copAt(tau);
      kneeLoad = Math.max(kneeLoad, g.ph); short = Math.max(short, g.short);
    }
  }
  let yMax = -1e9, yMin = 1e9; for (let i = 0; i < N; i++) { yMax = Math.max(yMax, K.hipY[i]); yMin = Math.min(yMin, K.hipY[i]); }
  // foot clearance: lowest point of the swinging foot (heel or ball; toes can bend) from early swing until just before
  // touchdown. lateClear covers only the approach to landing, which the wrapper fixes by keeping the knee bent longer.
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
             kneeLoadDeg: kneeLoad / D, kneeTODeg: TO.ph[0] / D, kneeRateTD: TD.ph[1] / D, kneeRateTO: TO.ph[1] / D, reachShort: short,
             bob: yMax - yMin, hipMax: yMax, phMaxDeg: phMax, phTypical, minClear, lateClear, termKnee, a: P.a, ball: P.ball,
             landDeg: P.landDeg, landTypical: P.landTypical, contactLen: v * tc, tcTypical: P.tcTypical, tcMin: P.tcMin, tcMax: P.tcMax,
             tcLimited: P.tcLimited, pushBehind: v * tc - P.a };
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
// A runner landing close under the hip keeps the knee bent a little longer before touchdown so the foot doesn't hit
// the ground on the way in. If the landing approach would go more than DRAG_TOL into the ground (the same allowance the
// drag warning uses, since the foot shape is simplified), find (by bisection) the smallest
// late-swing knee bend between 14 and 50 degrees that keeps it within that.
function realisticKinematics(o, leg, T, N) {
  const P = realisticPlan(o, leg, T);
  const core = termKnee => realisticCore(o, leg, T, N, P, termKnee);
  const K = core(14);
  if (K.info.lateClear >= -DRAG_TOL) return K;
  const top = core(50);
  if (top.info.lateClear < -DRAG_TOL) return top;
  let lo = 14, hi = 50, best = top;
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
    S.Pa = new Float64Array(N).fill(NaN); S.tauA = new Float64Array(N).fill(NaN); S.Fx = new Float64Array(N); S.Fy = new Float64Array(N);
    // horizontal force: along the line from the centre of pressure toward the hip. The force timing (skew) already makes
    // braking and propulsion cancel; c tilts the line by whatever is left (zero unless the skew hit its limit).
    let num = 0, den = 0;
    for (let i = 0; i < N; i++) if (K.stance[i]) { const Fy = K.grf[i]; num += Fy * (-K.copX[i] / K.hipY[i]); den += Fy; }
    const c = den > 0 ? num / den : 0;
    let hp2 = 0, hn2 = 0, kp2 = 0, kn2 = 0, ap = 0, an = 0, sp2 = 0, sn2 = 0, ext = 0, sumCheck = 0;
    for (let i = 0; i < N; i++) {
      if (!K.stance[i]) continue;
      const Fy = o.M * G * K.grf[i], Fx = Fy * (-K.copX[i] / K.hipY[i] - c);
      S.Fx[i] = Fx / (o.M * G); S.Fy[i] = K.grf[i];
      const th = K.th[i], hipY = K.hipY[i];
      const kx = L1 * Math.sin(th), ky = hipY - L1 * Math.cos(th);
      const axp = K.footX[i], ayp = K.ankleZ[i], cx = K.copX[i];
      const cross = (rx, ry) => -(rx * Fy - ry * Fx);        // moment on the distal part from the proximal one: -(r_cop - r_joint) x F
      const Mh = cross(cx - 0, 0 - hipY), Mk = cross(cx - kx, 0 - ky), Ma = cross(cx - axp, 0 - ayp);
      const wT = K.thd[i], wS = K.thd[i] - K.phd[i], wF = K.footAngRate[i];
      const Ph = Mh * wT, Pk = Mk * (wS - wT), Pa = Ma * (wF - wS);
      S.tauH[i] = Mh; S.tauK[i] = Mk; S.tauA[i] = Ma; S.Ph[i] = Ph; S.Pk[i] = Pk; S.Pa[i] = Pa;
      if (Ph > 0) hp2 += Ph * dt; else hn2 -= Ph * dt;
      if (Pk > 0) kp2 += Pk * dt; else kn2 -= Pk * dt;
      if (Pa > 0) ap += Pa * dt; else an -= Pa * dt;
      const Ps = Ph + Pk + Pa; if (Ps > 0) sp2 += Ps * dt; else sn2 -= Ps * dt;
      // the foot point under the centre of pressure never moves, so joint power must equal the ground force's power on the hip
      const Pext = Fx * o.speed + Fy * K.hipVy[i];
      sumCheck += Math.abs(Ps - Pext) * dt; ext += Math.abs(Pext) * dt;
    }
    const r0 = o.elastic;
    const perStep = o.share ? Math.max(sp2 - r0 * sn2, 0)
                            : Math.max(hp2 - r0 * hn2, 0) + Math.max(kp2 - r0 * kn2, 0) + Math.max(ap - r0 * an, 0);
    const mechS = perStep * f * 2;
    stance = { hipPos: hp2, hipNeg: hn2, kneePos: kp2, kneeNeg: kn2, anklePos: ap, ankleNeg: an, pos: sp2, neg: sn2, net: sp2 - sn2, perStep,
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
