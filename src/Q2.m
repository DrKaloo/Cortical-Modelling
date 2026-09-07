%% Q2 - percentage of trials where population 1 reaches threshold first
% Kaloyan Todorov

clc, clearvars
rng(1)   % same method as Q1, so these are the same 100 trials

I1 = 0.3310e-9;
I2 = 0.3304e-9;

T        = 3;
dt       = 0.1e-3;
N_trials = 100;

decisions = zeros(N_trials, 1);

for trial = 1:N_trials
    [~, ~, ~, ~, decisions(trial), ~, ~] = simulate_trial(I1, I2, T, dt);
end

n1   = sum(decisions == 1);
n2   = sum(decisions == 2);
nnd  = sum(decisions == 0);
pct1 = 100 * n1 / N_trials;

fprintf('\n Q2 results \n')
fprintf('Population 1 wins : %d / %d  (%.1f %%)\n', n1,  N_trials, pct1)
fprintf('Population 2 wins : %d / %d  (%.1f %%)\n', n2,  N_trials, 100*n2/N_trials)
fprintf('Undecided         : %d / %d  (%.1f %%)\n', nnd, N_trials, 100*nnd/N_trials)
