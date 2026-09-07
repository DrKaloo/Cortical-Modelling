%% Q6 - chronometric curve using the same 4 picks from Q4
% Kaloyan Todorov
% Error bars are the SEM of mean DT across choice-1 trials at each point.

clc, clearvars, close all

% For marking and loading purposes (check)
if ~exist('Q4_results.mat', 'file')
    error('Q4_results.mat not found. Run Q4.m first.')
end


load('Q4_results.mat', 'picked_stim', 'picked_dt', 'picked_dt_sem')

figure('Color', 'w', 'Position', [100 100 700 450]); hold on
errorbar(picked_stim, 1000*picked_dt, 1000*picked_dt_sem, 'o-', ...
         'LineWidth', 1.5, 'MarkerSize', 9, ...
         'MarkerFaceColor', 'r', 'Color', 'r', 'CapSize', 8)

xlabel("Stimulus strength (I_1 - I_2) [A]")
ylabel("Mean decision time for choice 1 (ms)")
title("Chronometric curve - Kaloyan")
grid on
set(gca, 'FontSize', 11)

saveas(gcf, 'Q6_chronometric_kaloyan.png')

fprintf('\n Q6: Chronometric curve \n')
fprintf('Stimulus strength (A)   Mean DT (ms)   SEM (ms)\n')
for i = 1:length(picked_stim)
    fprintf('  %+.3e           %4.0f           %.1f\n', ...
            picked_stim(i), 1000*picked_dt(i), 1000*picked_dt_sem(i))
end
