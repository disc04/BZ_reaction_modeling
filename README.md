# Belousov–Zhabotinsky Reaction Modeling

Numerical simulation of the **Belousov–Zhabotinsky (BZ) reaction**, the classic
oscillating chemical reaction, using the **Oregonator** kinetic model.

<p align="center">
  <img src="figures/bz_spiral.gif" width="300" alt="Spiral wave">
  <img src="figures/bz_targets.gif" width="300" alt="Target waves">
</p>
<p align="center"><sub>Simulated spiral and target waves in an excitable BZ medium (one-minute loops).</sub></p>

![Oregonator oscillations](figures/oregonator_timeseries.png)

## Background

In the BZ reaction, bromate oxidizes an organic substrate (malonic acid) in acid,
catalysed by a redox couple (Ce³⁺/Ce⁴⁺ or ferroin/ferriin). Instead of relaxing
monotonically to equilibrium, the intermediates oscillate for tens of minutes
(visible as red ↔ blue colour switching with ferroin).

The **Oregonator** (Field & Noyes, 1974) reduces the full FKN mechanism to five
steps and three intermediates:

| Step | Reaction | Rate |
|------|----------|------|
| O1 | A + Y → X + P | k₁·A·Y |
| O2 | X + Y → 2P | k₂·X·Y |
| O3 | A + X → 2X + 2Z | k₃·A·X |
| O4 | 2X → A + P | k₄·X² |
| O5 | B + Z → ½f·Y | k_c·B·Z |

A = BrO₃⁻, B = malonic acid, **X = HBrO₂** (activator), **Y = Br⁻** (inhibitor),
**Z = Ce⁴⁺** (oxidised catalyst). A and B are held constant.

Rate equations (mass-action kinetics):

```
dx/dt = k1·A·y − k2·x·y + k3·A·x − 2·k4·x²
dy/dt = −k1·A·y − k2·x·y + ½·f·kc·B·z
dz/dt = 2·k3·A·x − kc·B·z
```

**How it oscillates:** autocatalytic growth of HBrO₂ (O3) is suppressed while
Br⁻ is high. Once Br⁻ falls below a threshold, HBrO₂ explodes and oxidises the
catalyst (Z). Z then slowly regenerates Br⁻ (O5), shutting autocatalysis off
again. The delay introduced by Z is essential: eliminating it removes the
oscillations.

Conventions follow Field & Noyes (1974): step O3 is the net of
BrO₃⁻ + HBrO₂ → 2BrO₂• and 2BrO₂• + 2Ce³⁺ → 2HBrO₂ + 2Ce⁴⁺, so bromate is the
reactant and two Ce⁴⁺ are formed per event.

### Why the Oregonator?

| Model | Variables | Derived from BZ chemistry? | Stable limit cycle? | Main drawback |
|---|---|---|---|---|
| Lotka–Volterra | 2 | No | No (conservative centres) | Each perturbation gives a new orbit; not realistic chemistry |
| Brusselator | 2 | No (abstract; built for the Bray–Liebhafsky reaction) | Yes | Trimolecular step that is chemically improbable |
| Berlinator | 3 | Partly | No (stable nodes and a saddle) | Does not oscillate on its own |
| Full FKN mechanism | ~10 species | Yes | Yes | Many uncertain rate constants, very stiff, hard to analyse |
| Extended Oregonator (Showalter) | 4+ | Yes | Yes | More parameters for little qualitative gain |
| **Oregonator** | **3 (2 after reduction)** | **Yes (reduction of FKN)** | **Yes** | **Lumped parameters f and k_c** |

The Oregonator is the smallest model derived directly from the real BZ
mechanism. Its rate constants are measured, and its variables map onto actual
species (HBrO₂, Br⁻, Ce⁴⁺). It is also the only model above that covers this
project's whole scope in one framework: sustained oscillations with a
realistic period, reduction to two variables for nullcline and bifurcation
analysis, and excitability, which gives target and spiral waves. It has long
been the standard reference model for BZ dynamics, so results can be checked
against the literature.

**Limitations:** f and k_c are lumped, semi-empirical parameters; the total
amount of catalyst is not conserved; and the model does not reproduce the
aperiodic (chaotic) behaviour seen in some BZ experiments.

### Oscillation period

With the default Field–Noyes parameters (f = 1) the period is ≈ 390 s. The
period is set mainly by the slow catalyst-reduction step O5: after each spike,
Ce⁴⁺ must decay with rate constant k_c·B before Br⁻ falls low enough to let
autocatalysis restart, giving period ≈ 8 / (k_c·B). k_c is a lumped, empirical
constant, so the absolute period is a tunable quantity; experimental periods
(tens of seconds to minutes) depend on [BrO₃⁻], [malonic acid], [H⁺] and the
catalyst.

### The stoichiometric factor f

f is the number of Br⁻ ions returned per two Ce⁴⁺ reduced in O5. It lumps the
complex organic chemistry of Process C into one number and sets the **strength
of the negative feedback**:

- **f < ≈ 0.5**: too little Br⁻ is regenerated to switch autocatalysis off. The
  system sits in an **oxidised steady state** (high HBrO₂ and Ce⁴⁺, low Br⁻).
- **≈ 0.5 < f < ≈ 2.4**: oscillations. The theoretical limit is
  0.5 < f < 1 + √2 as k_c → 0. The period is shortest near f ≈ 1 and grows
  toward both edges.
- **f > ≈ 2.4**: Br⁻ is regenerated so strongly that autocatalysis never
  restarts. The system sits in a **reduced steady state** (low HBrO₂ and Ce⁴⁺, Br⁻ held above the switching threshold).

Both edges are Hopf bifurcations.

![Time series for different f](figures/oregonator_f_timeseries.png)

![Period vs f](figures/oregonator_f_period.png)

### Phase portrait: the limit cycle

![Phase portrait](figures/oregonator_phase_portrait.png)

In (x, y, z) phase space the oscillation is a closed orbit, a **limit cycle**.
Trajectories started from very different concentrations (grey), including a
tiny kick away from the steady state, all converge onto the same orbit. It is an
**attractor**: amplitude and period are fixed by the chemistry, not by the
initial conditions. This is the key difference from the Lotka–Volterra model,
whose closed orbits are conservative centres that shift with every perturbation.

The steady state at the centre is **unstable**. Its Jacobian eigenvalues at
f = 1 are λ ≈ −14.2, +1.69 and +5.1·10⁻⁴ s⁻¹; the two positive values push
trajectories away from it and onto the cycle.

The Br⁻ vs Ce⁴⁺ projection maps the cycle onto the FKN processes. Br⁻ drops
below its threshold (A), autocatalysis oxidises the catalyst (B), and Ce⁴⁺
reduction releases Br⁻, followed by a slow joint decay that resets the clock (C).
The very fast jumps and the slow drift are the signature of a **relaxation
oscillator**.

### Reduced two-variable Oregonator: nullclines and Hopf bifurcations

After scaling (Tyson, 1985), Br⁻ turns out to be by far the fastest variable
(ε′ ≈ 2·10⁻⁵ vs ε ≈ 10⁻²), so it can be set to its quasi-steady value
y = f·z/(q + x). This leaves two variables (Tyson & Fife, 1980):

```
ε dx/dτ = x(1 − x) − f·z·(x − q)/(x + q)
  dz/dτ = x − z
```

with x = scaled [HBrO₂], z = scaled [Ce⁴⁺], τ = k_c·B·t, ε = k_c·B/(k₃·A) and
q = 2k₁k₄/(k₂k₃). The reduction reproduces the full model: the steady state is
identical, and the period is 372 s vs 387 s.

![Nullclines](figures/reduced_nullclines.png)

The x-nullcline is **N-shaped**; the z-nullcline is the line z = x. Their
intersection is the steady state, and f moves it along the N:
- **f = 0.4:** the intersection sits on the right (outer) branch, a stable
  oxidised state (blue for ferroin).
- **f = 1:** it sits on the middle branch, which repels. The trajectory creeps
  along the slow outer branches and jumps quickly between them. This is the
  geometric picture of a relaxation oscillator.
- **f = 2.6:** the intersection sits on the left branch, a stable reduced state (red for ferroin).

![Bifurcation diagram](figures/reduced_bifurcation.png)

The trace of the 2×2 Jacobian changes sign at two **Hopf points**,
**f = 0.505 and f = 2.397**, matching the window found by brute-force
simulation of the full model (0.525–2.375). The two edges behave differently:
- **Lower edge (subcritical Hopf):** full-size oscillations appear abruptly.
  For 0.503 < f < 0.505 the stable steady state and the large cycle coexist
  (**bistability**), and the outcome depends on how strongly the system is
  perturbed.
- **Upper edge (supercritical Hopf + canard explosion):** a small cycle is
  born at f = 2.397 and grows smoothly. Near f ≈ 2.386 its amplitude jumps by
  about four orders of magnitude within Δf ≈ 0.001. This is a **canard
  explosion**, typical of slow–fast systems.

### Spatial patterns: target and spiral waves

In an unstirred layer, each point runs its own BZ kinetics and diffusion of
HBrO₂ couples neighbouring points. Adding diffusion to the reduced Oregonator
gives a reaction–diffusion system:

```
∂x/∂t = [x(1 − x) − f·z·(x − q)/(x + q)] / ε + Dₓ ∇²x
∂z/∂t = x − z                                + D_z ∇²z
```

Parameters (ε = 0.05, q = 0.002, f = 2.5, Dₓ = 1, D_z = 0.6) put the medium in
the **excitable** regime: f is just above the upper Hopf point, so the rest
state is stable, but a supra-threshold kick triggers a full excursion.
- **Trigger waves:** autocatalysis at the front plus diffusion of HBrO₂ ignite
  the neighbouring medium. The refractory tail (high Ce⁴⁺ and Br⁻) behind the
  front prevents backward propagation, so colliding waves annihilate.
- **Targets:** small discs with oscillatory kinetics (f = 1) act as pacemakers,
  like dust or bubbles in the dish, and emit concentric rings.
- **Spirals:** a planar wave cut in half leaves a free end, which curls up into
  a self-sustained rotating spiral. It needs no pacemaker.

![Target waves](figures/bz_targets.png)
![Spiral wave](figures/bz_spiral.png)

Numerics: a 400 × 400 grid, an isotropic 9-point Laplacian with no-flux
boundaries, and explicit Euler stepping in pure NumPy. Colours mimic the
ferroin indicator (red = reduced, blue = oxidised).

## Part 2: validation on experimental droplet data

**The study.** Chang, de Planque & Zauner (2018) recorded ferroin-catalysed BZ
oscillations in droplets in oil by video and quantified the colour waves. The
open dataset ([doi:10.5258/SOTON/D0363](https://doi.org/10.5258/SOTON/D0363),
CC-BY) contains:
- concentration series in which one component is varied at a time (H₂SO₄,
  NaBrO₃, malonic acid, ferroin), with oscillation frequency, lifetime, wave
  count and amplitude for 4–7 repeat droplets each;
- the full time series of one droplet over its oscillatory lifetime
  (289 curated peaks).

**What the data means for our model.** Acid and bromate enter the Oregonator
directly through the rate constants (k₁, k₂, k₃ ∝ [H⁺]; A = [BrO₃⁻]), so their
effect on frequency is a direct prediction. Malonic acid enters through B and,
indirectly, through the bromide yield f. Measured frequency falls with
malonic acid, the opposite of the simple prediction period ≈ 8/(k_c·B), which
makes it a discriminating test. Ferroin has no explicit counterpart in the
model, because the Oregonator does not conserve catalyst, so it is out of scope.

**Rationale.** Part 1 shows that the model is internally consistent; Part 2
asks whether it is right. Because f and k_c are lumped parameters, the aim is to
reproduce **trends and orders of magnitude**, not exact waveforms, and to show
clearly where the model's simplifications break down.

### Test 1: acid and bromate

![Acid and bromate test](figures/validation_acid_bromate.png)

The rate constants carry their acid orders (k₁ ∝ h², k₂ ∝ h, k₃ ∝ h; [H⁺] from
H₂SO₄ including HSO₄⁻ dissociation). The only free parameter, k_c, was
calibrated once on the base recipe (0.5 M H₂SO₄, 0.47 M NaBrO₃, 0.18 M malonic
acid; 0.106 Hz → k_c ≈ 34 M⁻¹ s⁻¹). Everything else is a prediction.
- **Model A (k_c constant) fails.** Below about 0.33 M acid or bromate it
  predicts no oscillations at all, and above that the frequency hardly
  changes, whereas the droplets oscillate across the whole range and speed up
  about 20×.
- **Model B (k_c ∝ [H⁺][BrO₃⁻]) captures the acid series** (log-log slope 1.47
  vs 1.66 measured) **but underestimates the bromate effect** (slope 1.00 vs
  1.73).

The catalyst-reduction rate must therefore rise with acid and bromate, which
is consistent with its dependence on bromomalonic acid, whose formation is
acid- and bromate-driven. The steeper-than-linear bromate dependence is
something the lumped Oregonator does not capture. Caveat: the data are median
frequencies over a droplet's lifetime as reactants deplete, while the model
uses the initial recipe.

### Test 2: malonic acid

![Malonic acid test](figures/validation_malonic_acid.png)

In the Oregonator, malonic acid (B) enters **only through the catalyst-reset
rate k_c·B** (step O5), so the whole question is how that product changes.
- **(a) Predictions.** With k_c constant, more malonic acid speeds the
  oscillations up, and above 0.28 M they stop altogether. With k_c·B constant,
  malonic acid has no effect. The droplets do neither: their frequency falls
  2.5× between 0.18 and 0.54 M.
- **(b) What the data would require.** Matching each measured frequency
  needs k_c·B to *drop* about 14× (log-log slope −2.3) as malonic acid triples,
  the opposite of what a reset rate driven by the organic substrate should do.
- **(c) The role of f.** The frequency map shows contours that run almost
  parallel to the f axis: f mainly decides *whether* the system oscillates,
  not how fast. No value of f rescues model A above 0.28 M, so the bromide
  yield cannot explain the slowing either.

**Conclusion:** the malonic acid trend lies outside what the Oregonator can
represent. The likely missing chemistry is malonic acid's role as a bromine
scavenger and source of bromomalonic acid, which controls bromide production
over the whole cycle. Extended mechanisms that keep bromomalonic acid as a
separate variable (e.g. Györgyi & Field, 1992) would be needed. Effects specific
to droplets, such as bromine partitioning into the oil phase, may also contribute.

### Test 3: time series of one droplet

![Time-series test](figures/validation_timeseries.png)

The droplet (base recipe) oscillates for 70 min in three phases: **initial**
(small, fast, ~5 s period), **main** (large, period lengthening from 12.7 to
about 17 s) and **late** (period suddenly shortening to about 12 s while the
amplitude collapses) before the medium is exhausted. Because the classic
Oregonator holds the reactants constant and cannot age, bromate and malonic
acid are made dynamic here (`bz_models/depletion.py`). Two parameters are
calibrated: k_c0 on the period at 10–20 min, and a depletion scale on the
70-min lifetime.
- **Main phase: reproduced.** The slow lengthening of the period from 10 to
  58 min follows the data closely, a genuine prediction since only the start
  point and the end time were fitted.
- **Initial phase: missed.** The model oscillates at full size from the first
  second. The induction chemistry (bromination of malonic acid, which first
  raises and then lowers bromide) is not in the model.
- **Late phase: wrong route out.** The model ends when malonic acid runs out,
  with the period diverging and the amplitude growing. The droplet ends the
  opposite way: shorter period and collapsing amplitude, the signature of
  approaching a supercritical Hopf bifurcation (the small-cycle regime seen in
  Part 1).
- **Amplitude and waveform: limited by catalyst conservation.** Model ferriin
  peaks grow far beyond the 2 mM of ferroin present, and its oxidation spikes
  (about 1 s) are much narrower than the droplet's (fast rise, roughly 8–10 s
  decay; the 2.5 s video frames blur only part of this).

### Part 2 summary

The Oregonator captures the core oscillator, the acid dependence (once k_c
scales with acid and bromate) and the slow ageing of the main phase. It fails
where its lumped chemistry matters: the malonic acid trend, the induction and
exhaustion phases, and absolute amplitudes. The two most direct improvements
would be a catalyst-conserving variant and a model that tracks bromomalonic
acid explicitly (e.g. Györgyi & Field, 1992).

## Usage

Install the package in editable mode (ideally inside a virtual environment).
This makes `bz_models` importable from anywhere, including notebooks, and code
changes take effect without reinstalling:

```bash
pip install -e .        # or: pip install -r requirements.txt
```

Then run any script from the repository root:

```bash
python scripts/run_oregonator.py   # single run, f = 1
python scripts/explore_f.py        # scan over f
python scripts/phase_portrait.py   # limit cycle in phase space
python scripts/reduced_oregonator.py  # nullclines and Hopf bifurcations
python scripts/spatial_patterns.py    # target and spiral waves (~2–3 min)
python scripts/validate_acid_bromate.py  # Part 2, test 1 (needs Chang_BZ_data/, ~1 min)
python scripts/validate_malonic_acid.py  # Part 2, test 2 (~5 min)
python scripts/validate_timeseries.py    # Part 2, test 3 (~3 min)
```

```python
from bz_models import OregonatorParams, plot_timeseries, simulate

t, c = simulate(OregonatorParams(f=1.5), t_end=2000)
plot_timeseries(t, c)
```

The system is stiff (rate constants span ~8 orders of magnitude), so it is
integrated with LSODA and tight tolerances.

## Project structure

```
pyproject.toml        package metadata and dependencies
bz_models/            model code (installable package)
  oregonator.py       rate equations, solver, steady state, Jacobian, period
  reduced.py          two-variable Oregonator, nullclines, Hopf points
  spatial.py          2D reaction-diffusion: targets and spirals
  recipe.py           rate constants from a BZ recipe ([H+], bromate, malonic acid)
  data.py             loader for the Chang et al. droplet dataset
  validation.py       shared calibration of k_c on the base recipe
  depletion.py        Oregonator with slow consumption of bromate and malonic acid
  plotting.py         plotting utilities
scripts/
  run_oregonator.py   simulate and plot time series
  explore_f.py        effect of f on dynamics and period
  phase_portrait.py   limit cycle in (x, y, z) space
  reduced_oregonator.py  nullclines and bifurcation diagram
  spatial_patterns.py    target and spiral waves, GIF animations
  validate_acid_bromate.py  Part 2, test 1: frequency vs acid and bromate
  validate_malonic_acid.py  Part 2, test 2: malonic acid and the role of f
  validate_timeseries.py    Part 2, test 3: one droplet over its lifetime
figures/              generated plots
```

## Roadmap

**Part 1: model and dynamics**
- [x] Three-variable Oregonator: sustained oscillations
- [x] Effect of the stoichiometric factor f
- [x] Phase portrait / limit cycle
- [x] Two-variable reduced Oregonator: nullclines, Hopf bifurcation in f
- [x] Reaction–diffusion: target and spiral waves in 2D

**Part 2: validation on droplet data** ([doi:10.5258/SOTON/D0363](https://doi.org/10.5258/SOTON/D0363))
- [x] Acid and bromate test: predicted vs measured frequency as [H₂SO₄] and [NaBrO₃] vary
- [x] Malonic acid test: frequency vs [malonic acid], and the role of f
- [x] Time-series test: waveform and slowing of oscillations as reactants deplete

## References

1. Field, R. J., Kőrös, E. & Noyes, R. M. (1972). *J. Am. Chem. Soc.* 94, 8649–8664.
2. Field, R. J. & Noyes, R. M. (1974). *J. Chem. Phys.* 60, 1877–1884.
3. Tyson, J. J. & Fife, P. C. (1980). Target patterns in a realistic model of
   the Belousov–Zhabotinskii reaction. *J. Chem. Phys.* 73, 2224–2237.
4. Tyson, J. J. (1985). A quantitative account of oscillations, bistability, and
   travelling waves in the Belousov–Zhabotinskii reaction. In *Oscillations and
   Traveling Waves in Chemical Systems* (Field & Burger, eds.), Wiley.
5. Field, R. J. (2007). Oregonator. *Scholarpedia* 2(5), 1386.
6. Jahnke, W. & Winfree, A. T. (1991). A survey of spiral-wave behaviors in the
   Oregonator model. *Int. J. Bifurcation Chaos* 1, 445–466.
7. Vîlcu, R. & Bala, D. (2004). Models of oscillating chemical reactions.
   *Analele Universității din București – Chimie* XIII, 277–286.
8. Chang, K. M., de Planque, M. R. R. & Zauner, K.-P. (2018). Towards functional
   droplet architectures: a Belousov–Zhabotinsky medium for networks.
   *Sci. Rep.* 8, 12656. Data: [doi:10.5258/SOTON/D0363](https://doi.org/10.5258/SOTON/D0363).
9. Györgyi, L. & Field, R. J. (1992). A three-variable model of deterministic
   chaos in the Belousov–Zhabotinsky reaction. *Nature* 355, 808–810.
