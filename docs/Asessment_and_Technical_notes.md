# Assessment and technical notes

## Assessment outcome

Original assessment mark: **87/100**.

The visible per-question marks in the returned document are:

| Question | Mark |
|---|---:|
| Q1 | 1.0 / 1.0 |
| Q2 | 0.5 / 1.0 |
| Q3 | 0.5 / 1.0 |
| Q4 | 1.0 / 1.0 |
| Q5 | 0.8 / 1.0 |
| Q6 | 1.0 / 1.0 |
| Q7 | 1.7 / 2.0 |
| Q8 | 1.8 / 2.0 |

These visible marks sum to 8.3/10 before any extra-credit contribution. The optional Jacobian section was marked "Perfect!", so the final 87/100 is consistent with extra credit contributing beyond the base question total. This is an inference from the marked document, not an independently supplied grading formula.

## Marker feedback

### Q2/Q3: threshold detection

The marker identified the `if` / `elseif` threshold-crossing implementation as biased.

Why this matters: if both populations cross the threshold during the same discrete Euler step, `S1` is checked first, so Choice 1 wins automatically. A more defensible implementation should detect crossing in both populations independently and compare estimated within-step crossing times, for example by linear interpolation between the previous and current values.

### Methods beyond taught material

The marker also noted that methods not covered in class should be explained or referenced. The phase-plane code uses numerical zero-contour extraction, a coarse fixed-point search and Newton refinement. These are legitimate methods, but a public scientific repository should state clearly which parts are core taught methods and which are extensions.

### AI-reference comment

The returned document contains the comment that there were no references to AI. The document itself does not establish whether AI was used in the assessed work. Do not add a retroactive AI-use declaration unless it is factually true. If AI was used and university rules required disclosure, any public archival note should be accurate about what was used and when.

### Q5: presentation

The marker considered the plot to contain too much unused space. Because the selected probabilities lie roughly between 0.54 and 1.0, a tighter y-range would make the four-point relationship easier to inspect while retaining a scientifically honest scale.

### Q7: numerical stability

The marker praised the analysis but wanted stability tested by initializing the dynamical system close to each fixed point and observing whether trajectories are attracted or repelled. The main missing test concerned FP2, the saddle.

A stronger numerical test would:

1. perturb each fixed point in several directions;
2. integrate the deterministic system from each perturbed initial condition;
3. compare distance from the fixed point over time;
4. show contraction in all local directions around stable nodes;
5. show contraction along one local direction and expansion along another around the saddle.

### Q7: vector field

The quiver arrows are useful, but the report should explain that each arrow represents the local deterministic derivative vector

```math
\mathbf{v}(S_1,S_2)=
\begin{bmatrix}
\dfrac{dS_1}{dt} \\
\dfrac{dS_2}{dt}
\end{bmatrix}.
```
at that phase-plane location. Arrow direction shows local flow direction; arrow magnitude reflects local speed unless the vectors are explicitly normalized.

### Q7 Jacobian

The optional analytical Jacobian analysis received the strongest feedback ("Perfect!"). It is a major portfolio strength because it connects numerical phase-plane structure to local linear stability through eigenvalues.

### Q8: model fit

The marker praised the figure and presentation but noted that the fit could improve.

The largest mismatch is psychometric rather than chronometric: the model becomes too accurate at intermediate coherences. The submitted discussion already recognizes this rather than hiding it, which is scientifically preferable.

## Additional technical review for GitHub

These points were not all explicit in the marker comments but matter for a public research-code repository.

### 1. Q1 reproducibility mismatch

`Q1.m` resets the random seed and then calls `simulate_trial` once only to obtain the time-vector length:

```matlab
[~, ~, ~, ~, ~, ~, t] = simulate_trial(I1, I2, T, dt);
```

That call consumes the random numbers for an entire trial. The subsequent 100 stored trials therefore begin from the *second* random trial in the seeded stream.

`Q2.m` and `Q3.m`, by contrast, reset the same seed and immediately begin collecting trials. Their first collected trial is therefore the trial discarded by Q1.

Consequence: the statement that Q1, Q2 and Q3 use exactly the same 100 trials is not strictly correct.

Clean fix: construct the time vector directly:

```matlab
t = (0:dt:T)';
Nt = numel(t);
```

or store the first simulation instead of discarding it.

### 2. Anonymous transfer-function guard

`Q7_phase_plane.m`, `Q7_jacobian.m` and `Q8.m` use an anonymous function of the form

```matlab
(abs(d*x) < 1e-6) .* taylor_term + ...
(abs(d*x) >= 1e-6) .* exact_term
```

This looks like a branch but both algebraic terms may still be evaluated before masking. At an exact singular point, the nominally "masked" exact expression can generate `NaN`, and `0 * NaN` remains `NaN`.

The safer implementation is an ordinary helper function that evaluates the Taylor approximation and exact formula on separate indexed subsets.

### 3. Q8 zero-coherence plotting

The Q8 figure sets `XScale` to `log` while `coherences` contains zero. MATLAB cannot place $x=0$ on a logarithmic axis. The script then supplies a tick at $x=1$ labelled `"0"`, but the actual 0%-coherence observation is still absent from its true coordinate.

Use a linear x-axis for the public figure, or explicitly transform the coordinates and document the transformation.

### 4. Q8 zero-coherence accuracy

The script fixes 0%-coherence accuracy to 0.5 by convention rather than estimating it from the simulated choices. This is defensible for a "probability correct" curve when no motion direction is objectively correct at zero coherence, but it should be stated as a convention.

A useful diagnostic is to report the simulated fraction choosing population 1 at 0% coherence separately; large deviation from 0.5 would reveal implementation or symmetry problems.

### 5. Fit metric interpretation

The reported combined error adds accuracy MSE to reaction-time MSE after dividing RT residuals by 200 ms. This is a custom objective, not a standard goodness-of-fit statistic. The README labels it accordingly.

### 6. Parameter fitting

The submitted Q8 parameters are a manually selected, biologically motivated set. The code does not perform formal parameter optimization. A future revision should call this a "manual fit" or "parameter tuning", not an optimized fit.

## Public-repository recommendation

For a portfolio repository, retain the scientific results and code, but do not upload:

- the marked DOCX;
- the student number;
- verbatim assignment questions;
- tutor annotations as screenshots;
- university-owned diagrams or brief text;
- generated `.mat` files unless they are genuinely needed as source data.

Keep the repository private if the assignment is still reused for assessment or if university policy restricts public release of solutions.

