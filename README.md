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

**Stance (quasi-static).** A spring-mass vertical ground force supports body weight over each step. It acts at the centre of pressure, which moves from just behind the ankle to the ball of the foot. Its horizontal part is set so braking and propulsion cancel over each contact. Hip, knee and ankle moments come from that force with the leg treated as weightless. Summed joint power matches the power delivered to the body to within about 4%.

**Realistic stride geometry.** During contact, the thigh angle and hip height are solved so the planted foot never slides, and the foot rolls 40° over the ball before push-off. Between steps the body follows a free-fall arc. If the foot would hit the ground on its way in to landing, the knee stays bent a little longer before touchdown.

## How it compares with lab data (average man, 3 m/s, 180 spm)

| Measure | Model | Published |
|---|---|---|
| Total mechanical work | 5.2 W/kg | 0.40–0.50 kcal/kg/km ≈ 5.0–6.3 W/kg (Cavagna) |
| Work moving the body | 3.6 W/kg (stance) | 0.25 kcal/kg/km ≈ 3.1 W/kg (Cavagna) |
| Work moving the limbs | 1.6 W/kg (legs only) | ≈ 2.4 W/kg incl. arms (Cavagna et al. 1997 equation) |
| Joint shares of positive work | ankle 38%, hip 41%, knee 20% | ankle 47%, hip 32% at 2 m/s (Farris & Sawicki 2012) |
| Net metabolic cost | 20.9 W/kg with no tendon return | ≈ 11–12 W/kg (≈ 3.7–4 J/kg/m) |

The mechanical work lines up well. Food energy is about twice too high with no tendon return, because real running recycles energy elastically. About 35% tendon return matches measured costs. Leg swing's share of the cost (about 25–31%) is at the high end of lab estimates (7–26%). Ankle work keeps rising at high speed, whereas measured ankle work levels off near 1.3 J/kg by about 5 m/s.

## Development

The page is built from two source files:

- `src/physics.js`: the simulation engine (no DOM)
- `src/template.html`: the page, with a `/*PHYSICS*/` placeholder

```
python3 build.py      # writes index.html
node tests/check.js   # sanity checks against known results
```

## Sources

- CDC, Anthropometric Reference Data for Children and Adults: United States, 2021–2023
- de Leva P. (1996). Adjustments to Zatsiorsky-Seluyanov's segment inertia parameters. *J Biomech* 29(9).
- Winter D.A. *Biomechanics and Motor Control of Human Movement.*
- Cavagna G.A. et al. Mechanical work in running; Cavagna, Mantovani, Willems & Musch (1997), *Pflügers Arch* 434.
- Schache A.G., Brown N.A.T., Pandy M.G. (2015). *J Exp Biol* 218(15).
- Farris D.J., Sawicki G.S. (2012). *J R Soc Interface* 9(66).
- Modica J.R., Kram R. (2005). *J Appl Physiol* 98(6); Arellano C.J., Kram R. (2014). *Integr Comp Biol* 54(6).
