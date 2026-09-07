%% Q5 - psychometric curve using the 4 picks from Q4
% Kaloyan Todorov
% Error bars are binomial SE = sqrt(p(1-p)/n), with n = 300.

clc, clearvars, close all

% For marking and loading purposes (check)
if ~exist('Q4_results.mat', 'file')
    error('Q4_results.mat not found. Run Q4.m first.')
end



load('Q4_results.mat', 'picked_stim', 'picked_pct', 'picked_n_trials')

p  = picked_pct / 100;
n  = picked_n_trials;
se = sqrt(p .* (1 - p) ./ n);

figure('Color', 'w', 'Position', [100 100 700 450]); hold on
errorbar(picked_stim, p, se, 'o-', ...
         'LineWidth', 1.5, 'MarkerSize', 9, ...
         'MarkerFaceColor', 'b', 'Color', 'b', 'CapSize', 8)

xlabel("Stimulus strength (I_1 - I_2) [A]")
ylabel("Probability of choice 1")
title("Psychometric curve - Kaloyan")
ylim([0 1.05])
grid on
set(gca, 'FontSize', 11)

saveas(gcf, 'Q5_psychometric_kaloyan.png')

fprintf('\n Q5: Psychometric curve \n')
fprintf('Stimulus strength (A)   P(choice 1)   binomial SE\n')
for i = 1:length(picked_stim)
    fprintf('  %+.3e           %.3f         %.3f\n', ...
            picked_stim(i), p(i), se(i))
end
