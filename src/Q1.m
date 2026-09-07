%% Q1 - plot trajectories of 100 trials
% Kaloyan Todorov

clc, clearvars, close all
rng(1)

I1 = 0.3310e-9;
I2 = 0.3304e-9;

T        = 3;
dt       = 0.1e-3;
N_trials = 100;

% One run just to get the time vector size
[~, ~, ~, ~, ~, ~, t] = simulate_trial(I1, I2, T, dt);
Nt = length(t);

S1_all = zeros(Nt, N_trials);
S2_all = zeros(Nt, N_trials);

for trial = 1:N_trials
    [S1_all(:,trial), S2_all(:,trial), ~, ~, ~, ~, ~] = ...
        simulate_trial(I1, I2, T, dt);
end

% Plot - S2 first (red), then S1 (blue) so choice 1 sits on top
figure('Color', 'w', 'Position', [100 100 800 450]); hold on
for trial = 1:N_trials
    plot(t, S2_all(:,trial), 'Color', [1 0 0 0.25], 'LineWidth', 0.5);
end
for trial = 1:N_trials
    plot(t, S1_all(:,trial), 'Color', [0 0 1 0.25], 'LineWidth', 0.5);
end
yline(0.5, '--k', 'LineWidth', 1.2, 'Label', 'threshold')

% Handles so the legend picks the right colours
h1 = plot(NaN, NaN, 'Color', 'b', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'Color', 'r', 'LineWidth', 1.5);
legend([h1 h2], {'S_1 (choice 1)', 'S_2 (choice 2)'}, ...
       'Location', 'southeast', 'AutoUpdate', 'off')

xlabel("Time (s)")
ylabel("Synaptic gating S")
title("100-trial decision network trajectories - Kaloyan")
ylim([0 1])
box on
set(gca, 'FontSize', 11)

saveas(gcf, 'Q1_trajectories_kaloyan.png')
