# Running Efficiency Sim

An interactive simulator for the mechanical work and energy cost of running. It animates a two-segment leg and works out the hip, knee and ankle torques at every instant of the stride, then adds up the work the muscles do. Swing (leg in the air) and stance (foot on the ground) are measured separately.

Open `index.html` in a browser. It is a single self-contained file with no build step and no dependencies apart from Google Fonts.

## What you can change

- **Runner:** average US man or woman (CDC 2021–2023), or any height and body mass
- **Leg motion:**
  - *Realistic stride:* the foot plants on the ground, the body bounces, and the knee follows typical running motion
  - *Straight*, *Bend, then straighten*, *Stay bent:* the hip hangs in the air and the leg swings like a pendulum
- **Realistic stride controls:** running speed, cadence, peak knee bend in swing, landing point (degrees in front of the hip) and ground contact time
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

**Stance (quasi-static).** The ground force acts at the centre of pressure: under the heel while the foot lowers, moving to the ball while it is flat, and on the ball while the heel lifts. It points from there at the hip, and `κ` is chosen so braking and propulsion cancel over each contact (steady speed). Hip, knee and ankle moments come from that force with the leg treated as weightless. Because the force passes through the hip, the hip does no work in stance; all modelled hip work is in the swing. Summed joint power equals the power the force delivers to the body, and the body ends each step with the energy it started with.

**Realistic stride geometry.** The foot is rigid and never slides. It lands on the heel with the toes 12° up and the knee at 20°, lowers flat, then rolls over the ball as the heel lifts (up to 65°). During contact the knee and thigh angles are whatever the leg needs to reach from the hip to the ankle. The foot lowers at the rate that gives a typical knee speed (300°/s) just after landing. The heel lifts just enough to keep the knee from straightening past a floor that ends near 22° at toe-off. Swing curves match each joint's angle and speed at toe-off and touchdown. If the foot would hit the ground on its way in to landing, the knee stays bent a little longer before touchdown.

## How it compares with lab data (average man, 3 m/s, 180 spm)

| Measure | Model | Published |
|---|---|---|
| Total mechanical work | 4.8 W/kg | 0.40–0.50 kcal/kg/km ≈ 5.0–6.3 W/kg (Cavagna) |
| Work moving the body | 3.2 W/kg (stance) | 0.25 kcal/kg/km ≈ 3.1 W/kg (Cavagna) |
| Work moving the limbs | 1.6 W/kg (legs only) | ≈ 2.4 W/kg incl. arms (Cavagna et al. 1997 equation) |
| Joint shares of positive work | ankle 33%, hip 30%, knee 37% | ankle 47%, hip 32% at 2 m/s (Farris & Sawicki 2012) |
| Net metabolic cost | 19.2 W/kg with no tendon return | ≈ 11–12 W/kg (≈ 3.7–4 J/kg/m) |

The mechanical work lines up well, a little under Cavagna's total. Food energy is well above measured costs with no tendon return, because real running recycles energy elastically. About 40% tendon return matches measured costs. Leg swing's share of the cost is 33% at 3 m/s and 20–24% from 4 to 6 m/s, at or above the high end of lab estimates (7–26%). The knee does more of the work than in the lab data, and the hip does none of its work in stance (see above).

Two known limits. Swing cost barely grows above 3 m/s (about 115–145 W), whereas real swing cost rises with speed: the thigh's range in swing is set by the landing and push-off positions, which scale with leg length rather than speed. Ankle work matches measured values at 5 m/s (1.3 J/kg) but keeps rising beyond that, whereas measured ankle work levels off.

## Development

The page is built from two source files:

- `src/physics.js`: the simulation engine (no DOM)
- `src/template.html`: the page, with a `/*PHYSICS*/` placeholder

```
python build.py       # writes index.html
node tests/check.js   # sanity checks: known results, energy balance, geometry, all settings
```

## Sources

- CDC, Anthropometric Reference Data for Children and Adults: United States, 2021–2023
- de Leva P. (1996). Adjustments to Zatsiorsky-Seluyanov's segment inertia parameters. *J Biomech* 29(9).
- Winter D.A. *Biomechanics and Motor Control of Human Movement.*
- Cavagna G.A. et al. Mechanical work in running; Cavagna, Mantovani, Willems & Musch (1997), *Pflügers Arch* 434.
- Schache A.G., Brown N.A.T., Pandy M.G. (2015). *J Exp Biol* 218(15).
- Farris D.J., Sawicki G.S. (2012). *J R Soc Interface* 9(66).
- Modica J.R., Kram R. (2005). *J Appl Physiol* 98(6); Arellano C.J., Kram R. (2014). *Integr Comp Biol* 54(6).
