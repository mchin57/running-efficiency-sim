# Running Efficiency Sim

An interactive simulator for the mechanical work and energy cost of running. It animates a two-segment leg and works out the hip, knee and ankle torques at every instant of the stride, then adds up the work the muscles do. Swing (leg in the air) and stance (foot on the ground) are measured separately.

Open `index.html` in a browser. It is a single self-contained file with no build step and no dependencies apart from Google Fonts.

## What you can change

- **Runner:** average US man or woman (CDC 2021–2023), or any height and body mass
- **Leg motion:**
  - *Realistic stride:* the foot plants on the ground, the body bounces, and the knee follows typical running motion
  - *Straight*, *Bend, then straighten*, *Stay bent:* the hip hangs in the air and the leg swings like a pendulum
- **Realistic stride controls:** running speed, cadence, peak knee bend in swing, landing point (degrees in front of the hip), and ground contact time and ground contact length (one setting shown both ways: length = speed × time; set either, and that one stays put when speed changes)
- **Typical values for the speed:** an arrow under the knee-bend, landing, contact and cadence sliders marks the typical value for the current speed (from measured runners); click it to use that value
- **Physics:** gravity, whether hip and knee can share energy, tendon energy return, muscle efficiency

The page shows swing power, stance power, the whole-stride total and food energy. It also has live joint power readouts, power, torque and angle charts over one stride, and a chart of swing cost against knee bend that marks where the foot would drag.

## How it works

**Swing (inverse dynamics).** The leg is a two-segment pendulum: thigh, then shank plus foot. Segment masses, centres of mass and inertias come from de Leva (1996); segment lengths come from Winter's ratios of height. The motion is prescribed, and the equations of motion give the hip and knee torques:

```
Q1 = M11·θ̈ + M12·ψ̈ + m2·L1·c2·sin(θ−ψ)·ψ̇² + (m1·c1 + m2·L1)·g·sin θ
Q2 = M12·θ̈ + M22·ψ̈ − m2·L1·c2·sin(θ−ψ)·θ̇² + m2·c2·g·sin ψ
knee torque = Q2,  hip torque = Q1 + Q2,  power = torque × joint angular velocity
```

In the realistic stride, the hip's vertical acceleration acts on the swinging leg as extra or reduced gravity. That is about 2–3 g while the other foot is pushing and zero in the air.

**Body motion.** The vertical ground force is a spring-mass pulse that supports body weight over each step, `F/m = Fpk·[sin(πs) + κ·sin(2πs)]` for `s` from 0 to 1 through the contact and zero in the air. The `κ` term shifts the force earlier or later without changing its impulse. The body's height comes from that force in closed form and repeats every step, so the drawn hip, the swinging leg and the planted leg all share one motion that obeys F = ma.

**Stance (quasi-static).** The ground force acts at the centre of pressure: under the heel while the foot lowers, moving to the ball while it is flat, out along the flat toes to the tip while the heel lifts, and on the tip while the toes roll. It points from there at the hip, and `κ` is chosen so braking and propulsion cancel over each contact (steady speed). `κ` can't go past ±0.5 (the force would go negative), so contact times where no `κ` balances are left out of the allowed contact-time range. Likewise, if the leg can't reach the ground at push-off even with the shortest contact, the landing point moves forward until it can. If no contact time balances with the chosen landing point, the force is tilted off the hip instead and the page flags it as not steady running. Hip, knee, ankle and toe-joint moments come from that force with the leg treated as weightless. Because the force passes through the hip, the hip does no work in stance (except in flagged settings); all modelled hip work is in the swing. Summed joint power equals the power the force delivers to the body, and the body ends each step with the energy it started with.

**Realistic stride geometry.** The foot (a rigid rearfoot plus toes, the front 30%) never slides. It lands on the heel with the toes 12° up and the knee at 10°. If the leg would lock straight just after landing (a foot landing almost under the hip, or a very short contact), it lands with more knee bend instead. The foot lowers flat, then the rearfoot turns about the ball as the heel lifts (up to 65°) while the pressure moves out to the toe tip. If the leg still can't reach, the toes roll up about their tip (up to 45°), as far as needed. Rolling off the toes is what lets contact reach realistic lengths. The toe joint mostly absorbs energy as it bends (about 25 J per step at 3 m/s). During contact the knee and thigh angles are whatever the leg needs to reach from the hip to the ankle. The foot lowers at the rate that gives a typical knee speed (300°/s) just after landing. The heel starts to lift at 55% of contact, or earlier if, after its deepest bend, the knee would straighten past 30°; it lifts just enough to keep the knee from straightening past a floor that ends near 22° at toe-off. Swing curves match each joint's angle and speed at toe-off and touchdown, the toe-off speeds taken from inside the contact; in between, the thigh swings forward to a measured, speed-dependent peak (see below) and the knee to a peak bend of 60° + 12°·v (capped at 130°). If the foot would scuff the ground on its way in to landing, the thigh swings up to 15° further forward, which lifts the foot.

**Typical values.** The defaults for landing point and contact, and the arrows on the sliders, follow measured runners. From Clark, Meng & Stearne 2020 (3.1–10 m/s): thigh angle at touchdown 0.033·v + 0.335 rad (which sets the landing point, about 18–23° from 3 to 6 m/s), contact time 0.78·v^−0.89 s scaled with height from a 1.75 m runner, and peak forward thigh swing 0.102·v + 0.273 rad. Cadence: 133.9 + 10.9·v steps/min, a straight-line fit to 28 recreational runners at 2.5–4.5 m/s (Fukuchi et al. 2017, Table 4), about 167 steps/min at 3 m/s; it defaults to 180 rather than the typical value, and its arrow shows the typical one. At typical cadences the model reaches the typical contact time from 2.5 m/s up (at 2 m/s it reaches 366 ms of 421). At faster-than-typical cadences at slow speeds there isn't room between steps, so it uses the nearest it can.

## How it compares with lab data (average man, 3 m/s, 180 spm)

**Being re-derived.** The numbers in this section and in *Read with care* were computed before the touchdown and pressure-path changes of November 2026 and have not all been recomputed yet. For current, like-for-like comparisons against measured runners, run `node tests/validate.js`.

| Measure | Model | Published |
|---|---|---|
| Total mechanical work | 4.7 W/kg | 0.40–0.50 kcal/kg/km ≈ 5.0–6.3 W/kg (Cavagna) |
| Work moving the body | 2.9 W/kg (stance) | 0.25 kcal/kg/km ≈ 3.1 W/kg (Cavagna) |
| Work moving the limbs | 1.8 W/kg (legs only) | ≈ 2.4 W/kg incl. arms (Cavagna et al. 1997 equation) |
| Joint shares of positive work | ankle and toes 41%, hip 29%, knee 30% | ankle 47%, hip 32% at 2 m/s (Farris & Sawicki 2012) |
| Net metabolic cost | 18.8 W/kg with no tendon return | ≈ 11–12 W/kg (≈ 3.7–4 J/kg/m) |

The mechanical work lines up well, a little under Cavagna's total. Food energy is well above measured costs with no tendon return, because real running recycles energy elastically. About 40% tendon return matches measured costs. Leg swing's share of the cost is 38% at 3 m/s and 32–33% from 4 to 6 m/s at 180 steps/min, above the high end of lab estimates (7–26%); at the cadences runners actually use at 3 m/s it is far lower. The hip does none of its work in stance (see above).

**Swing and speed.** With the thigh's forward swing from the measured fit and contact close to measured, the thigh's motion matches the same study well: at its step rates, the peak backward swing is 19–28° against 24–28° measured (matching from 4 m/s up), the total thigh range is within 2° of measured from 3.5 m/s up, and the average thigh angular speed is within 7%. At those step rates swing cost rises steeply with speed: about 60 W at 3 m/s, 125 W at 4, 240 W at 5 and 420 W at 6. At a fixed 180 steps/min it rises less (about 160 W at 3 m/s to 240 W at 6), because shorter ground contact at speed leaves the leg more time to swing. The forward swing is always kept at least 4° past the thigh's touchdown angle, which only matters when the landing point is set well forward.

Per step, the split between joints is off in a consistent way. At 3.5 m/s (typical cadence) the model's ankle does +1.33 / −0.29 J/kg against +0.97 / −0.48 measured, and the knee +0.73 / −1.48 against +0.45 / −0.75 (Schache et al. 2015). The toe joint absorbs 0.31 J/kg against about 0.38 measured (Day & Hahn 2020).

## Read with care

The physics is checked (see `tests/check.js`), but three simplifications skew particular conclusions:

- **Joint shares.** The ground force is aimed through the hip, so the hip does no stance work and the knee and ankle do it all. Real hip muscles do a good share of stance work, so the knee and ankle look more important here than they are, with peak knee moments about 1.4–1.7 times measured values and peak ankle moments 25–55% above them. Whole-stride totals are less affected than the split between joints.
- **Cadence and swing.** The swing leg follows a set path, so a higher cadence covers the same path in less time and swing work climbs very steeply. Measured limb-moving work rises much more gently with step rate: the direction is right, the size is exaggerated.
- **Food energy.** One muscle efficiency (25% by default) converts all positive work, and braking work is free. Measured running looks far more efficient (apparent efficiencies of about 45–70% have been reported, Cavagna & Kaneko 1977, largely from tendon recoil), and braking does cost some energy. Hence the ~40% tendon return needed to match measured costs. Comparisons between settings are more trustworthy than absolute calories.

Also known and not yet fixed:

- **Force timing.** The vertical ground force peaks at 33–37% of contact, against about 43% measured (its force-weighted centre is about 45%). At typical settings the timing shift sits at or near its limit.
- **Push-off.** The heel lifts about 63° at toe-off, against about 40–49° measured. The final pressure point is the toe tip, whereas measured pressure ends under the big-toe pad, about 85% of the way along the foot.

## Possible model updates

Not started; listed for a future round. Each would move the default numbers.

1. **Push off from the big-toe pad.** Move the final pressure point and the toe pivot from the toe tip to about 85% of the way along the foot, where measured pressure ends. Expected: less heel lift at toe-off (now about 63°, measured about 40–49°) and a shorter lever at the toe joint.
2. **Let the ground force pass in front of and behind the hip.** For example, aim it at a centre of mass above the hip, or use a measured horizontal-force shape. Expected: the hip does some stance work, and knee moments and knee braking work drop (now about 1.4–1.7× and 2× published).
3. **Peak the force later.** Change the vertical force shape so it peaks nearer 43% of contact, as measured, rather than 33–37%.

Kept as they are, by choice: cadence defaults to 180 steps/min (its arrow shows the typical value), and settings that can't be steady running stay reachable, with a warning on the page.

## Development

The page is built from two source files:

- `src/physics.js`: the simulation engine (no DOM)
- `src/template.html`: the page, with a `/*PHYSICS*/` placeholder

```
python build.py       # writes index.html
node tests/check.js   # sanity checks: known results, energy balance, geometry, all settings

node tools/fetch-fukuchi.js    # downloads the measured running data into data/ (git-ignored), checking size and md5
node tools/reduce-fukuchi.js   # reduces it to tests/data/fukuchi2017-summary.json (means and SDs across subjects)
node tests/validate.js         # compares the model with the measurements (--cohort, --speed, --set key=value, --json)
```

The measured data are from Fukuchi, Fukuchi & Duarte (2017), PeerJ 5:e3298, doi 10.7717/peerj.3298 (Figshare doi 10.6084/m9.figshare.4543435), licensed CC BY 4.0.

## Sources

- CDC, Anthropometric Reference Data for Children and Adults: United States, 2021–2023
- de Leva P. (1996). Adjustments to Zatsiorsky-Seluyanov's segment inertia parameters. *J Biomech* 29(9).
- Winter D.A. *Biomechanics and Motor Control of Human Movement.*
- Cavagna G.A. et al. Mechanical work in running; Cavagna, Mantovani, Willems & Musch (1997), *Pflügers Arch* 434.
- Cavagna G.A., Kaneko M. (1977). Mechanical work and efficiency in level walking and running. *J Physiol* 268(2).
- Day E.M., Hahn M.E. (2020). *Gait Posture* 84:187–191.
- Fukuchi R.K., Fukuchi C.A., Duarte M. (2017). A public dataset of running biomechanics and the effects of running speed on lower extremity kinematics and kinetics. *PeerJ* 5:e3298.
- Clark K.P., Meng C.R., Stearne D.J. (2020). 'Whip from the hip': thigh angular motion, ground contact mechanics, and running speed. *Biology Open* 9(10), bio053546.
- Schache A.G., Brown N.A.T., Pandy M.G. (2015). *J Exp Biol* 218(15).
- Farris D.J., Sawicki G.S. (2012). *J R Soc Interface* 9(66).
- Modica J.R., Kram R. (2005). *J Appl Physiol* 98(6); Arellano C.J., Kram R. (2014). *Integr Comp Biol* 54(6).
