%% Q4 - find stimulus strengths giving ~50, 60, 80, 100% choice-1 rates
% Kaloyan Todorov
%
% Stimulus strength is defined (per the brief) as s = I1 - I2. I kept the
% sum (I1+I2)/2 fixed at I_mean so only the difference varies, which
% isolates stimulus strength as a single control variable.
%
% A 13-point grid of s values is scanned with 300 trials each, and for
% each target percentage I pick the grid point with the closest observed
% rate. Results are saved to Q4_results.mat for use by Q5 and Q6.

clc, clearvars, close all
rng(1)

I_mean = (0.3310e-9 + 0.3304e-9) / 2;   % 0.3307e-9 A

T        = 3;
dt       = 0.1e-3;
N_trials = 300;

% Grid concentrated where the psychometric transition happens
stim_grid = [ -1e-13, linspace(0, 1.2e-12, 10), 1.6e-12, 2.5e-12 ];
n_grid = length(stim_grid);

pct_choice1      = zeros(n_grid, 1);
mean_dt_c1       = nan(n_grid, 1);
all_dts_per_grid = cell(n_grid, 1);   % kept DTs for SEM in Q6

fprintf('\nScanning %d stimulus strengths, %d trials each: \n', n_grid, N_trials)

for g = 1:n_grid
    stim = stim_grid(g);
    I1   = I_mean + stim/2;
    I2   = I_mean - stim/2;

    decisions = zeros(N_trials, 1);
    dts       = nan(N_trials, 1);

    for trial = 1:N_trials
        [~, ~, ~, ~, decisions(trial), dts(trial), ~] = ...
            simulate_trial(I1, I2, T, dt);
    end

    pct_choice1(g) = 100 * sum(decisions == 1) / N_trials;
    mask_c1        = (decisions == 1);
    if any(mask_c1)
        mean_dt_c1(g) = mean(dts(mask_c1));
    end
    all_dts_per_grid{g} = dts(mask_c1);

    fprintf('  stim = %+.2e A   -> choice 1: %5.1f %%   mean DT (c1): %.0f ms\n', ...
            stim, pct_choice1(g), 1000 * mean_dt_c1(g))
end

% Pick grid point closest to each target
targets       = [50, 60, 80, 100];
picked_stim   = zeros(size(targets));
picked_pct    = zeros(size(targets));
picked_dt     = zeros(size(targets));
picked_dt_sem = zeros(size(targets));

fprintf('\n Q4 picks \n')
for i = 1:length(targets)
    [~, idx] = min(abs(pct_choice1 - targets(i)));
    picked_stim(i) = stim_grid(idx);
    picked_pct(i)  = pct_choice1(idx);
    picked_dt(i)   = mean_dt_c1(idx);

    dts_picked = all_dts_per_grid{idx};
    if length(dts_picked) > 1
        picked_dt_sem(i) = std(dts_picked) / sqrt(length(dts_picked));
    end

    fprintf('Q4_%d  target %3d%%  ->  stim = %+.3e A   (%.1f %%, mean DT = %.0f ms)\n', ...
            i, targets(i), picked_stim(i), picked_pct(i), 1000*picked_dt(i))
end

% Save everything Q5 and Q6 will need
picked_n_trials = N_trials * ones(size(targets));
picked_n_wins   = round(picked_pct * N_trials / 100);
save('Q4_results.mat', 'stim_grid', 'pct_choice1', 'mean_dt_c1', ...
                       'targets', 'picked_stim', 'picked_pct', 'picked_dt', ...
                       'picked_n_trials', 'picked_n_wins', 'picked_dt_sem')

fprintf('\nResults saved as Q4_results.mat\n')
