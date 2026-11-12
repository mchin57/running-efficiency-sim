// ---- Leg-swing physics: two-segment leg hung from the hip, prescribed motion, inverse dynamics ----
const G = 9.81;
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
      const Th = T / 2, s = t / Th, k = Math.PI / Th;
      K.ph[i] = phmax * Math.sin(Math.PI * s) ** 2;
      K.phd[i] = phmax * k * Math.sin(2 * Math.PI * s);
      K.phdd[i] = 2 * phmax * k * k * Math.cos(2 * Math.PI * s);
    }
    K.gEff[i] = o.gravity ? G : 0; K.count[i] = 1;
  }
  return K;
}
function blankK(N) {
  const z = () => new Float64Array(N);
  return { th: z(), thd: z(), thdd: z(), ph: z(), phd: z(), phdd: z(), gEff: z(), count: z(), stance: z(),
           hipY: z(), footX: z(), grf: z(), ankleLift: z(), footPitch: z(), copX: z(), ankleZ: z() };
}

// ---------- realistic stride: real ground, planted foot that rolls over the ball, bouncing body ----------
// Knee follows a smooth curve through typical running values. During contact the thigh angle and hip height
// are solved so the foot stays planted on the ground. The body's vertical acceleration (felt by the swinging
// leg as extra or reduced gravity) comes from a spring-mass ground force: F/m = Fpk sin(pi t/tc) in contact, 0 in flight.
// Typical landing: ankle touches down 0.256 leg lengths ahead of the hip and leaves 0.544 leg lengths behind it (0.8 leg lengths of contact).
const REAL = { landAhead: 0.256, pushBehind: 0.544, kneeTD: 20, rollStart: 0.45, rollDeg: 40 };
function realisticCore(o, leg, T, N, termKnee) {
  const v = o.speed, Ts = T / 2;
  const { L1, Ls, ha, footLen } = leg;
  // landing point: angle from vertical of the line from the hip to the ground under the ankle at touchdown
  const dTD = Math.sqrt(L1 * L1 + Ls * Ls + 2 * L1 * Ls * Math.cos(REAL.kneeTD * Math.PI / 180));
  const hipAtTD = x => Math.sqrt(Math.max(dTD * dTD - x * x, 1e-6)) + ha;
  const aTyp = REAL.landAhead * leg.legLen;
  const landTypical = Math.atan(aTyp / hipAtTD(aTyp)) * 180 / Math.PI;
  const landDeg = o.landDeg != null ? o.landDeg : landTypical;
  let a = aTyp;
  for (let it = 0; it < 60; it++) a = Math.min(0.9 * dTD, Math.tan(landDeg * Math.PI / 180) * hipAtTD(a));
  // contact time: typical = foot leaves the usual distance behind the hip. A set time moves the push-off point instead.
  // Limits: some time in the air (one foot down at a time), the foot must at least pass under the hip,
  // and it cannot push off more than 0.62 leg lengths behind the hip.
  const tcMax = Math.min(0.95 * Ts, (a + 0.62 * leg.legLen) / v), tcMin = Math.min(Math.max(0.06, (a + 0.05) / v), tcMax);
  const tcTypical = Math.max(tcMin, Math.min((a + REAL.pushBehind * leg.legLen) / v, 0.9 * Ts));
  const tcWanted = o.contactMs != null ? o.contactMs / 1000 : tcTypical;
  const tc = Math.max(tcMin, Math.min(tcMax, tcWanted));
  const tcLimited = Math.abs(tc - tcWanted) > 5e-4 ? (tc > tcWanted ? 'min' : 'max') : null;
  const tf = Ts - tc, Tsw = T - tc;
  const Fpk = Math.PI * G * Ts / (2 * tc), v0 = G * tf / 2, k = Math.PI / tc;
  const bodyAcc = tl => { tl = ((tl % Ts) + Ts) % Ts; return tl < tc ? Fpk * Math.sin(k * tl) - G : -G; };
  const grfAt = tl => { tl = ((tl % Ts) + Ts) % Ts; return tl < tc ? Fpk * Math.sin(k * tl) / G : 0; };
  const D = Math.PI / 180;
  // knee knots (time from toe-off, degrees): toe-off, swing peak, terminal swing, touchdown, stance peak, late stance
  const phTypical = Math.min(130, 60 + 12 * v);
  const phMax = o.realKnee != null ? o.realKnee : phTypical;
  const knots = [[0, 22], [0.42 * Tsw, phMax], [0.88 * Tsw, termKnee], [Tsw, REAL.kneeTD], [Tsw + 0.40 * tc, Math.min(52, 38 + 2.5 * v)], [Tsw + 0.80 * tc, 33]];
  const knee = periodicSpline(knots.map(p => p[0]), knots.map(p => p[1] * D), T);
  const ball = 0.45 * footLen;
  const sstep = s => { s = Math.min(1, Math.max(0, s)); return s * s * s * (s * (6 * s - 15) + 10); };
  function ankle(tau) {   // ankle position relative to the point on the ground below the hip
    const al = REAL.rollDeg * D * sstep((tau / tc - REAL.rollStart) / (1 - REAL.rollStart));
    return { X: a - v * tau + ball * (1 - Math.cos(al)), Z: ha + ball * Math.sin(al), al };
  }
  function stanceGeom(t) { // thigh angle and hip height that keep the foot planted
    const tau = t - Tsw, A = ankle(tau), ph = knee.at(t)[0];
    let th = Math.atan2(A.X, 0.8) + ph * 0.5;
    for (let it = 0; it < 30; it++) {
      const fx = L1 * Math.sin(th) + Ls * Math.sin(th - ph) - A.X;
      const dfx = L1 * Math.cos(th) + Ls * Math.cos(th - ph);
      th -= fx / dfx; if (Math.abs(fx) < 1e-12) break;
    }
    return { th, ph, hip: L1 * Math.cos(th) + Ls * Math.cos(th - ph) + A.Z, X: A.X, al: A.al };
  }
  const h = 1e-4;
  const deriv = (fn, t) => { const m = fn(t - h), c = fn(t), p = fn(t + h); return [c, (p - m) / (2 * h), (p - 2 * c + m) / (h * h)]; };
  const thS = t => stanceGeom(t).th;
  const TO = deriv(thS, T - h), TD = deriv(thS, Tsw + h);
  const thMax = TD[0] + (4 + 1.5 * v) * D, sTh = 0.72;
  const thVia = [thMax, 0, -2 * (thMax - TD[0]) / ((1 - sTh) * Tsw) ** 2];
  const thSegs = [quintic(TO, thVia, sTh * Tsw), quintic(thVia, TD, (1 - sTh) * Tsw)];
  // hip height for drawing: planted-leg geometry in contact, smooth bridge through flight
  const hipTD = stanceGeom(Tsw + h), hipTO = stanceGeom(T - h);
  const vUp = (hipTD.hip - hipTO.hip + G * tf * tf / 2) / tf;     // free-fall arc between toe-off and the next touchdown
  const hipAt = t => { const tl = (((t - Tsw) % Ts) + Ts) % Ts; if (tl < tc) return stanceGeom(Tsw + tl).hip;
    const u = tl - tc; return hipTO.hip + vUp * u - G * u * u / 2; };
  const K = blankK(N);
  for (let i = 0; i < N; i++) {
    const t = i * T / N, ph = knee.at(t);
    K.ph[i] = ph[0]; K.phd[i] = ph[1]; K.phdd[i] = ph[2];
    K.hipY[i] = hipAt(t); K.gEff[i] = G + bodyAcc(t - Tsw);
    if (t < Tsw) {
      const th = t < sTh * Tsw ? evalQ(thSegs[0], t) : evalQ(thSegs[1], t - sTh * Tsw);
      K.th[i] = th[0]; K.thd[i] = th[1]; K.thdd[i] = th[2]; K.count[i] = 1;
      // foot pitch (toe-down +): square to the shin early in swing, levelling out for landing late in swing
      const wgt = sstep((t / Tsw - 0.45) / 0.45);
      K.footPitch[i] = (1 - wgt) * (-(th[0] - K.ph[i]));
    } else {
      const th = deriv(thS, t), g = stanceGeom(t);
      K.th[i] = th[0]; K.thd[i] = th[1]; K.thdd[i] = th[2];
      K.stance[i] = 1; K.footX[i] = g.X; K.ankleLift[i] = g.al; K.ankleZ[i] = ha + ball * Math.sin(g.al);
      // centre of pressure: moves from just behind the ankle to the ball of the foot, then stays on the ball while the heel lifts
      const sc = (t - Tsw) / tc, x0 = g.X - ball * (1 - Math.cos(g.al));
      K.copX[i] = x0 + (sc < REAL.rollStart ? -0.10 * footLen + (ball + 0.10 * footLen) * (sc / REAL.rollStart) : ball);
    }
    K.grf[i] = grfAt(t - Tsw);
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
  K.info = { tc, tf, Tsw, stepLen: v * Ts, peakGRF: Fpk / G, bob: yMax - yMin, hipMax: yMax, phMaxDeg: phMax, phTypical, minClear, lateClear, termKnee, a, ball, landDeg, landTypical, contactLen: v * tc, tcTypical, tcMin, tcMax, tcLimited, pushBehind: v * tc - a };
  return K;
}
// periodic cubic spline through (x_i, y_i), period P; returns at(t) -> [y, y', y'']
function periodicSpline(xs, ys, P) {
  const n = xs.length, hh = [], M = new Array(n).fill(0);
  for (let i = 0; i < n; i++) hh.push(((xs[(i + 1) % n] - xs[i]) + P) % P || P);
  // solve cyclic tridiagonal system for second derivatives by Gauss elimination on a dense matrix (n is small)
  const A = [], b = [];
  for (let i = 0; i < n; i++) {
    const row = new Array(n).fill(0), hp = hh[(i - 1 + n) % n], hn = hh[i];
    row[(i - 1 + n) % n] += hp / 6; row[i] += (hp + hn) / 3; row[(i + 1) % n] += hn / 6;
    A.push(row); b.push((ys[(i + 1) % n] - ys[i]) / hn - (ys[i] - ys[(i - 1 + n) % n]) / hp);
  }
  const Ms = solve(A, b);
  return { at(t) {
    t = ((t % P) + P) % P;
    let i = n - 1; for (let j = 0; j < n; j++) { if (t >= xs[j] && (j === n - 1 || t < xs[j + 1])) { i = j; break; } }
    if (t < xs[0]) i = n - 1;
    const x0 = xs[i], hn = hh[i], y0 = ys[i], y1 = ys[(i + 1) % n], m0 = Ms[i], m1 = Ms[(i + 1) % n];
    let u = t - x0; if (u < 0) u += P;
    const w = hn - u;
    const y = m0 * w ** 3 / (6 * hn) + m1 * u ** 3 / (6 * hn) + (y0 - m0 * hn * hn / 6) * w / hn + (y1 - m1 * hn * hn / 6) * u / hn;
    const yd = -m0 * w * w / (2 * hn) + m1 * u * u / (2 * hn) - (y0 - m0 * hn * hn / 6) / hn + (y1 - m1 * hn * hn / 6) / hn;
    const ydd = m0 * w / hn + m1 * u / hn;
    return [y, yd, ydd];
  } };
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
// the ground on the way in. If the landing approach would go more than 1 cm into the ground (the same allowance the
// drag warning uses, since the foot shape is simplified), find (by bisection) the smallest
// late-swing knee bend between 14 and 50 degrees that keeps it within that.
function realisticKinematics(o, leg, T, N) {
  let K = realisticCore(o, leg, T, N, 14);
  if (K.info.lateClear >= -0.008) return K;
  const top = realisticCore(o, leg, T, N, 50);
  if (top.info.lateClear < -0.008) return top;
  let lo = 14, hi = 50, best = top;
  for (let it = 0; it < 6; it++) {
    const mid = (lo + hi) / 2, Km = realisticCore(o, leg, T, N, mid);
    if (Km.info.lateClear >= -0.008) { hi = mid; best = Km; } else lo = mid;
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
  let hp = 0, hn = 0, kp = 0, kn = 0, sp = 0, sn = 0, peakKnee = 0, peakHip = 0;
  for (let i = 0; i < N; i++) {
    if (!K.count[i]) { S.tauH[i] = NaN; S.tauK[i] = NaN; S.Ph[i] = NaN; S.Pk[i] = NaN; continue; }
    const th = K.th[i], thd = K.thd[i], thdd = K.thdd[i], g = K.gEff[i];
    const ps = th - K.ph[i], psd = thd - K.phd[i], psdd = thdd - K.phdd[i];
    const d = K.ph[i];
    const M12 = m2 * L1 * c2 * Math.cos(d), sd = Math.sin(d);
    const Q1 = M11 * thdd + M12 * psdd + m2 * L1 * c2 * sd * psd * psd + (m1 * c1 + m2 * L1) * g * Math.sin(th);
    const Q2 = M12 * thdd + M22 * psdd - m2 * L1 * c2 * sd * thd * thd + m2 * c2 * g * Math.sin(ps);
    const tauK = Q2, tauH = Q1 + Q2;
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
    const { Ls, ha } = leg;
    // horizontal force: along the line from the centre of pressure toward the hip, tilted by a constant so the
    // horizontal push sums to zero over the contact (steady running speed)
    let num = 0, den = 0;
    for (let i = 0; i < N; i++) if (K.stance[i]) { const Fy = K.grf[i]; num += Fy * (-K.copX[i] / K.hipY[i]); den += Fy; }
    const c = den > 0 ? num / den : 0;
    let hp2 = 0, hn2 = 0, kp2 = 0, kn2 = 0, ap = 0, an = 0, sp2 = 0, sn2 = 0, ext = 0, sumCheck = 0;
    const prevI = i => (i - 1 + N) % N, nextI = i => (i + 1) % N;
    for (let i = 0; i < N; i++) {
      if (!K.stance[i]) continue;
      const Fy = o.M * G * K.grf[i], Fx = Fy * (-K.copX[i] / K.hipY[i] - c);
      S.Fx[i] = Fx / (o.M * G); S.Fy[i] = K.grf[i];
      const th = K.th[i], ps = th - K.ph[i], hipY = K.hipY[i];
      const kx = L1 * Math.sin(th), ky = hipY - L1 * Math.cos(th);
      const axp = K.footX[i], ayp = K.ankleZ[i], cx = K.copX[i];
      const cross = (rx, ry) => -(rx * Fy - ry * Fx);        // moment on the distal part from the proximal one: -(r_cop - r_joint) x F
      const Mh = cross(cx - 0, 0 - hipY), Mk = cross(cx - kx, 0 - ky), Ma = cross(cx - axp, 0 - ayp);
      const wT = K.thd[i], wS = K.thd[i] - K.phd[i];
      const ip = K.stance[prevI(i)] ? prevI(i) : i, inx = K.stance[nextI(i)] ? nextI(i) : i;
      const wF = -(K.ankleLift[inx] - K.ankleLift[ip]) / (((inx - ip + N) % N || 1) * dt);
      const Ph = Mh * wT, Pk = Mk * (wS - wT), Pa = Ma * (wF - wS);
      S.tauH[i] = Mh; S.tauK[i] = Mk; S.tauA[i] = Ma; S.Ph[i] = Ph; S.Pk[i] = Pk; S.Pa[i] = Pa;
      if (Ph > 0) hp2 += Ph * dt; else hn2 -= Ph * dt;
      if (Pk > 0) kp2 += Pk * dt; else kn2 -= Pk * dt;
      if (Pa > 0) ap += Pa * dt; else an -= Pa * dt;
      const Ps = Ph + Pk + Pa; if (Ps > 0) sp2 += Ps * dt; else sn2 -= Ps * dt;
      const yd = (K.hipY[nextI(i)] - K.hipY[prevI(i)]) / (2 * dt);
      sumCheck += Math.abs(Ps - (Fx * o.speed + Fy * yd)) * dt; ext += Math.abs(Fx * o.speed + Fy * yd) * dt;
    }
    const r0 = o.elastic;
    const perStep = o.share ? Math.max(sp2 - r0 * sn2, 0)
                            : Math.max(hp2 - r0 * hn2, 0) + Math.max(kp2 - r0 * kn2, 0) + Math.max(ap - r0 * an, 0);
    const mechS = perStep * f * 2;
    stance = { hipPos: hp2, hipNeg: hn2, kneePos: kp2, kneeNeg: kn2, anklePos: ap, ankleNeg: an, pos: sp2, neg: sn2, perStep,
               mech: mechS, metab: mechS / o.eff, balanceErr: ext > 0 ? sumCheck / ext : 0 };
  }
  const r = o.elastic;
  const perStride = o.share ? Math.max(sp - r * sn, 0)
                            : Math.max(hp - r * hn, 0) + Math.max(kp - r * kn, 0);
  const mech = perStride * f * 2;        // both legs, W
  const metab = mech / o.eff;
  return { stance, leg, S, K, info: K.info || null, T, f, N, dt, hipPos: hp, kneePos: kp, hipNeg: hn, kneeNeg: kn, perStride, mech, metab,
           kcalh: metab * 3600 / 4184, peakKneeDeg: peakKnee * 180 / Math.PI, peakHipDeg: peakHip * 180 / Math.PI };
}
