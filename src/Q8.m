%% Q8 - fit to Roitman & Shadlen (2002)
% Kaloyan Todorov
%
% Four biologically-based extensions to the base model (standard
% Wong & Wang 2006):
%
%   (1) Coherence-dependent inputs:  I1(c') = I0*(1 + kappa*c'/100),
%       I2(c') = I0*(1 - kappa*c'/100). Reflects the approximately linear
%       coherence tuning of MT neurons (Britten et al., 1993).
%
%   (2) Lower decision threshold S_thr = 0.45 (was 0.5). Trained monkeys
%       adopt reduced thresholds to balance speed vs accuracy
%       (Hanks et al., 2014).
%
%   (3) Non-decision time t_nd = 300 ms added to RT. Sensory + MT latency
%       ~100 ms (Britten et al., 1993) plus motor execution ~200 ms
%       (Huk & Shadlen, 2005).
%
%   (4) Doubled input noise sigma = 0.06 nA (was 0.03). Attentional
%       fluctuations + MT Poisson variability (Cohen & Kohn, 2011);
%       needed to reproduce the graded accuracy at mid coherences.

clc, clearvars, close all
rng(1)

% Fitted parameters
kappa = 0.03;
S_thr = 0.45;
t_nd  = 0.300;
sigma = 0.06e-9;

% Base model (unchanged from previous questions)
tau   = 0.100;
gamma = 0.641;
a     = 270e9;
b     = 108;
d     = 0.15;
W11   = 0.2609e-9;     W22 = 0.2609e-9;
W12   = 0.0497e-9;     W21 = 0.0497e-9;
I0    = 0.3307e-9;


% Roitman & Shadlen 2002 Fig 3 targets
coherences = [0, 3.2, 6.4, 12.8, 25.6, 51.2];
target_acc = [0.50, 0.58, 0.68, 0.83, 0.95, 0.99];
target_rt  = [860,  850,  790,  720,  620,  530];

n_coh    = length(coherences);
T        = 3;
dt       = 0.1e-3;
N_trials = 500;


% Response function
resp = @(x) (abs(d*x) < 1e-6) .* (1/d + x/2) + ...
            (abs(d*x) >= 1e-6) .* ( x ./ (1 - exp(-d*x)) );

model_acc = zeros(1, n_coh);
model_rt  = zeros(1, n_coh);

fprintf('\n Q8: Roitman & Shadlen fit (%d trials per coherence) \n', N_trials)

for ci = 1:n_coh
    c  = coherences(ci);
    I1 = I0 * (1 + kappa*c/100);
    I2 = I0 * (1 - kappa*c/100);

    decisions = zeros(N_trials, 1);
    dts       = nan(N_trials, 1);

    for trial = 1:N_trials
        S1 = 0; S2 = 0;
        for k = 1:round(T/dt)
            Isyn1 = W11*S1 - W12*S2 + I1 + sigma*randn;
            Isyn2 = W22*S2 - W21*S1 + I2 + sigma*randn;

            r1 = resp(a*Isyn1 - b);
            r2 = resp(a*Isyn2 - b);

            S1 = S1 + dt*(-S1/tau + (1 - S1)*gamma*r1);
            S2 = S2 + dt*(-S2/tau + (1 - S2)*gamma*r2);

            if S1 >= S_thr
                decisions(trial) = 1;
                dts(trial)       = k * dt;
                break
            elseif S2 >= S_thr
                decisions(trial) = 2;
                dts(trial)       = k * dt;
                break
            end
        end
    end

    % Accuracy: 0.5 by convention at c = 0, else P(correct) on decided trials
    if c == 0
        model_acc(ci) = 0.5;
        mean_dt       = mean(dts(~isnan(dts)));
    else
        correct       = (decisions == 1);
        n_decided     = sum(~isnan(dts));
        model_acc(ci) = sum(correct) / n_decided;
        mean_dt       = mean(dts(correct));
    end


    % Total RT includes non-decision time
    model_rt(ci) = 1000 * (mean_dt + t_nd);

    fprintf('  c = %5.1f%%   acc = %.3f  (target %.2f)   RT = %4.0f ms  (target %d ms)\n', ...
            c, model_acc(ci), target_acc(ci), model_rt(ci), target_rt(ci))
end


% Fit quality
err_acc   = mean( (model_acc - target_acc).^2 );
err_rt    = mean( ((model_rt - target_rt) / 200).^2 );
total_err = err_acc + err_rt;

fprintf('\nFit error:  accuracy MSE = %.4f,  RT MSE (normalised) = %.4f,  combined = %.4f\n', ...
        err_acc, err_rt, total_err)

% Side-by-side plot
figure('Color', 'w', 'Position', [100 100 1200 460]);

subplot(1, 2, 1); hold on
plot(coherences, target_acc, 'ko-', 'MarkerSize', 9, ...
     'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
     'DisplayName', 'Roitman & Shadlen 2002 (target)')
plot(coherences, model_acc, 'bs-', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b', 'LineWidth', 1.5, ...
     'DisplayName', 'Model fit')
xlabel("Motion coherence (%)")
ylabel("Probability correct")
title("Psychometric function - Kaloyan")
ylim([0.4 1.05]); xlim([0.9 60])
set(gca, 'XScale', 'log', 'FontSize', 11, 'GridAlpha', 0.15, ...
         'XTick',      [1, 3.2, 6.4, 12.8, 25.6, 51.2], ...
         'XTickLabel', {'0', '3.2', '6.4', '12.8', '25.6', '51.2'})
grid on; box on
legend('Location', 'southeast', 'FontSize', 9)

subplot(1, 2, 2); hold on
plot(coherences, target_rt, 'ko-', 'MarkerSize', 9, ...
     'MarkerFaceColor', 'w', 'LineWidth', 1.5, ...
     'DisplayName', 'Roitman & Shadlen 2002 (target)')
plot(coherences, model_rt, 'rs-', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'r', 'LineWidth', 1.5, ...
     'DisplayName', 'Model fit')
xlabel("Motion coherence (%)")
ylabel("Mean RT for correct trials (ms)")
title("Chronometric function - Kaloyan")
xlim([0.9 60])
set(gca, 'XScale', 'log', 'FontSize', 11, 'GridAlpha', 0.15, ...
         'XTick',      [1, 3.2, 6.4, 12.8, 25.6, 51.2], ...
         'XTickLabel', {'0', '3.2', '6.4', '12.8', '25.6', '51.2'})
grid on; box on
legend('Location', 'northeast', 'FontSize', 9)

saveas(gcf, 'Q8_roitman_shadlen_kaloyan.png')
save('Q8_results.mat', 'coherences', 'model_acc', 'model_rt', ...
     'target_acc', 'target_rt', 'kappa', 'S_thr', 't_nd', 'sigma', 'total_err')

fprintf('\nFigure saved: Q8_roitman_shadlen_kaloyan.png\n')
