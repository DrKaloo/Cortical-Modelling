%% Q7 - phase plane analysis (nullclines, fixed points, sample trajectories)
% Kaloyan Todorov
%
% Steps:
%   1. Evaluate dS1/dt and dS2/dt on a 400x400 grid
%   2. Draw both nullclines as the zero level sets (contour, level 0)
%   3. Locate fixed points: coarse grid scan + Newton refinement
%   4. Overlay one noisy trial per choice, cut off after decision + 800 ms
%
% Saves fp, I1, I2 to Q7_results.mat for the Jacobian script.

clc, clearvars, close all
rng(1)

% Parameters (same as simulate_trial.m)
tau   = 0.100;
gamma = 0.641;
a     = 270e9;
b     = 108;
d     = 0.15;
W11   = 0.2609e-9;
W22   = 0.2609e-9;
W12   = 0.0497e-9;
W21   = 0.0497e-9;

I1 = 0.3310e-9;
I2 = 0.3304e-9;

% Deterministic flow (no noise, for the structural analysis)
resp = @(x) (abs(d*x) < 1e-6) .* (1/d + x/2) + ...
            (abs(d*x) >= 1e-6) .* ( x ./ (1 - exp(-d*x)) );

dS1 = @(S1, S2) -S1/tau + (1 - S1) .* gamma .* ...
                resp(a*(W11*S1 - W12*S2 + I1) - b);
dS2 = @(S1, S2) -S2/tau + (1 - S2) .* gamma .* ...
                resp(a*(W22*S2 - W21*S1 + I2) - b);

% 400x400 grid over [0,1]^2
N = 400;
s = linspace(0, 1, N);
[S1g, S2g] = meshgrid(s, s);
dS1_grid = dS1(S1g, S2g);
dS2_grid = dS2(S1g, S2g);

% Fixed-point search: coarse scan for local minima of |dS1|+|dS2|
Ncoarse = 60;
sc = linspace(0.01, 0.99, Ncoarse);
[S1c, S2c] = meshgrid(sc, sc);
mag_coarse = abs(dS1(S1c, S2c)) + abs(dS2(S1c, S2c));

candidates = [];
for i = 2:Ncoarse-1
    for j = 2:Ncoarse-1
        local = mag_coarse(i-1:i+1, j-1:j+1);
        if mag_coarse(i, j) == min(local(:))
            candidates = [candidates; S1c(i, j), S2c(i, j)];
        end
    end
end

% Newton refined on each candidate
fp = [];
for k = 1:size(candidates, 1)
    x = candidates(k, :)';
    for iter = 1:60
        F = [dS1(x(1), x(2)); dS2(x(1), x(2))];
        if norm(F) < 1e-10, break; end

        h  = 1e-6;
        Jn = [(dS1(x(1)+h, x(2)) - dS1(x(1)-h, x(2)))/(2*h), ...
              (dS1(x(1), x(2)+h) - dS1(x(1), x(2)-h))/(2*h); ...
              (dS2(x(1)+h, x(2)) - dS2(x(1)-h, x(2)))/(2*h), ...
              (dS2(x(1), x(2)+h) - dS2(x(1), x(2)-h))/(2*h)];

        x = x - Jn \ F;
        if any(x < -0.1) || any(x > 1.1), break; end
    end

    converged = norm([dS1(x(1), x(2)); dS2(x(1), x(2))]) < 1e-6;
    inside    = all(x >= 0) && all(x <= 1);
    if converged && inside

        % Skip duplicates
        is_new = true;
        for m = 1:size(fp, 1)
            if norm(fp(m, :) - x') < 1e-3, is_new = false; break; end
        end
        if is_new, fp = [fp; x']; end
    end
end
fp = sortrows(fp);

fprintf('\n Q7: Fixed points \n')
for k = 1:size(fp, 1)
    fprintf('  FP%d: (S1, S2) = (%.4f, %.4f)\n', k, fp(k,1), fp(k,2))
end


% Sampling one noisy trial per choice
T  = 3;
dt = 0.1e-3;

S1_ch1 = []; S2_ch1 = [];
S1_ch2 = []; S2_ch2 = [];
got1 = false; got2 = false;
attempt = 0;

while (~got1 || ~got2) && attempt < 50
    attempt = attempt + 1;
    [S1, S2, ~, ~, decision, dt_dec, tt] = simulate_trial(I1, I2, T, dt);
    if isnan(dt_dec), continue; end

    % Keep trajectory up to decision + 800 ms so it settles at its attractor
    t_end   = min(dt_dec + 0.80, T);
    idx_end = find(tt <= t_end, 1, 'last');

    if decision == 1 && ~got1
        S1_ch1 = S1(1:idx_end);  S2_ch1 = S2(1:idx_end);  got1 = true;
    elseif decision == 2 && ~got2
        S1_ch2 = S1(1:idx_end);  S2_ch2 = S2(1:idx_end);  got2 = true;
    end
end
fprintf('Sampled %d trials to obtain one trajectory per choice\n', attempt)


% Plot
figure('Color', 'w', 'Position', [100 100 780 680]);
hold on


% Sparse quiver for direction field
qstep = 26;
qx = S1g(1:qstep:end, 1:qstep:end);
qy = S2g(1:qstep:end, 1:qstep:end);
qu = dS1_grid(1:qstep:end, 1:qstep:end);
qv = dS2_grid(1:qstep:end, 1:qstep:end);
quiver(qx, qy, qu, qv, 1.3, 'Color', [0.82 0.82 0.82], 'LineWidth', 0.6)

% Nullclines
contour(S1g, S2g, dS1_grid, [0 0], 'b', 'LineWidth', 2);
contour(S1g, S2g, dS2_grid, [0 0], 'r', 'LineWidth', 2);


% Trajectories
if got1
    plot(S1_ch1, S2_ch1, 'Color', [0 0.5 0], 'LineWidth', 1.3)
    plot(S1_ch1(1), S2_ch1(1), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0 0.5 0], 'MarkerEdgeColor', 'k')
end
if got2
    plot(S1_ch2, S2_ch2, 'Color', [0.6 0 0.6], 'LineWidth', 1.3)
    plot(S1_ch2(1), S2_ch2(1), 'o', 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.6 0 0.6], 'MarkerEdgeColor', 'k')
end

% Fixed points with labels
for k = 1:size(fp, 1)
    plot(fp(k,1), fp(k,2), 'ks', 'MarkerSize', 11, ...
         'MarkerFaceColor', 'k', 'LineWidth', 1)
    text(fp(k,1)+0.02, fp(k,2)+0.03, sprintf('FP%d', k), ...
         'FontSize', 11, 'FontWeight', 'bold')
end


% Handles for a clean legend
h_ds1 = plot(NaN, NaN, 'b', 'LineWidth', 2);
h_ds2 = plot(NaN, NaN, 'r', 'LineWidth', 2);
h_tr1 = plot(NaN, NaN, 'Color', [0 0.5 0],   'LineWidth', 1.3);
h_tr2 = plot(NaN, NaN, 'Color', [0.6 0 0.6], 'LineWidth', 1.3);
h_fp  = plot(NaN, NaN, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
legend([h_ds1, h_ds2, h_tr1, h_tr2, h_fp], ...
       {'dS_1/dt = 0  (S_1 nullcline)', ...
        'dS_2/dt = 0  (S_2 nullcline)', ...
        'Noisy trial -> Choice 1', ...
        'Noisy trial -> Choice 2', ...
        'Fixed points'}, ...
       'Location', 'northeast', 'FontSize', 9, 'AutoUpdate', 'off')

xlabel("S_1  (choice 1 gating variable)")
ylabel("S_2  (choice 2 gating variable)")
title("Phase plane analysis - Kaloyan")
xlim([0 1]); ylim([0 1])
axis square
grid on; box on
set(gca, 'FontSize', 11, 'GridAlpha', 0.15)

saveas(gcf, 'Q7_phase_plane_kaloyan.png')
save('Q7_results.mat', 'fp', 'I1', 'I2')
