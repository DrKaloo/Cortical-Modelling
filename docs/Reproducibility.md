# Reproducibility

## Software

The project is written in MATLAB.

The source uses:

- forward Euler integration;
- `rng` for deterministic seeding;
- `randn` for stochastic input;
- `contour` for nullclines;
- `quiver` for the phase-plane vector field;
- numerical Newton refinement for fixed points;
- `eig` for Jacobian stability;
- `.mat` files to pass generated results between scripts.

## Execution order

Run from a directory in which `src/` is on the MATLAB path:

```text
Q1
Q2
Q3
Q4
Q5
Q6
Q7_phase_plane
Q7_jacobian
Q8
```

`Q4.m` creates `Q4_results.mat`, which is loaded by Q5 and Q6.

`Q7_phase_plane.m` creates `Q7_results.mat`, which is loaded by `Q7_jacobian.m`.

## Randomness

The scripts call `rng(1)`, so a fixed MATLAB release and compatible random-number-generator behaviour should give deterministic repeated runs.

Exact cross-version reproducibility is not guaranteed because:

- random-number algorithms can differ across environments or settings;
- numerical library details can affect Newton refinement and eigenvalue reporting;
- plotting features can differ by MATLAB release.

There is also a known Q1 issue: a preliminary call to `simulate_trial` consumes one trial's random numbers before the 100 displayed trajectories are stored. This means Q1's stored 100-trial sequence is offset by one trial relative to Q2/Q3, even though all scripts use `rng(1)`.

## Generated files

The following are generated outputs and are ignored by `.gitignore`:

- `Q4_results.mat`
- `Q7_results.mat`
- `Q8_results.mat`
- MATLAB `.fig` files
- autosave files

The PNG figures in `figures/` are retained because they document the graded submission.

## Numerical interpretation

The model uses Gaussian input samples at every discrete timestep with the amplitude given in the original exercise. This should be interpreted as the discrete noise model used for this project rather than as a claim that the implementation is a mathematically normalized continuous-time stochastic differential equation.

