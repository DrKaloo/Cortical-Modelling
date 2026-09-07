%% Q3 - mean decision time on choice-1 trials
% Kaloyan Todorov

clc, clearvars
rng(1)

I1 = 0.3310e-9;
I2 = 0.3304e-9;

T        = 3;
dt       = 0.1e-3;
N_trials = 100;

decisions = zeros(N_trials, 1);
dts       = nan(N_trials, 1);

for trial = 1:N_trials
    [~, ~, ~, ~, decisions(trial), dts(trial), ~] = ...
        simulate_trial(I1, I2, T, dt);
end

% Mean crossing time on pop-1 winning trials
mask = (decisions == 1);
mdt  = mean(dts(mask));
medt = median(dts(mask));
sdt  = std(dts(mask));

fprintf('\n Q3 results \n')
fprintf('Trials where population 1 won : %d\n', sum(mask))
fprintf('Mean decision time            : %.3f s  (%.0f ms)\n', mdt,  1000*mdt)
fprintf('Median decision time          : %.3f s  (%.0f ms)\n', medt, 1000*medt)
fprintf('SD decision time              : %.3f s  (%.0f ms)\n', sdt,  1000*sdt)
