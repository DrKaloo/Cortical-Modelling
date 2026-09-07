function [S1, S2, r1, r2, decision, dt_decision, t] = simulate_trial(I1, I2, T, dt)
% Single trial of the two-population rate decision network (Wong & Wang 2006
% reduced). Integrates the noisy dynamics with forward Euler and returns
% the first threshold crossing (S >= 0.5).
%
% Inputs : I1, I2  external input currents for pop 1 and pop 2 (A)
%          T       total simulation time (s)          [default 3]
%          dt      integration step (s)               [default 1e-4]
%
% Outputs: S1, S2       synaptic gating trajectories
%          r1, r2       firing rate trajectories
%          decision     1 if S1 crossed first, 2 if S2, 0 if neither
%          dt_decision  time of first crossing (s), NaN if none
%          t            time vector
%
% Kaloyan Todorov

if nargin < 3, T  = 3;
end

if nargin < 4, dt = 0.1e-3;
end

% Parameters (from the brief)
tau   = 0.100;
gamma = 0.641;
a     = 270e9;
b     = 108;
d     = 0.15;
W11   = 0.2609e-9;
W22   = 0.2609e-9;
W12   = 0.0497e-9;
W21   = 0.0497e-9;
sigma = 0.03e-9;
thr   = 0.5;

t  = 0:dt:T;
N  = length(t);
S1 = zeros(N, 1);
S2 = zeros(N, 1);
r1 = zeros(N, 1);
r2 = zeros(N, 1);

decision    = 0;
dt_decision = NaN;

% Euler integration
for k = 2:N

    % Synaptic currents (Eq 5, 6)
    Isyn1 = W11*S1(k-1) - W12*S2(k-1) + I1 + sigma*randn;
    Isyn2 = W22*S2(k-1) - W21*S1(k-1) + I2 + sigma*randn;

    % Firing rates (Eq 3, 4); Taylor guard near x = 0 to avoid 0/0
    x1 = a*Isyn1 - b;
    x2 = a*Isyn2 - b;
    if abs(d*x1) < 1e-6
        r1(k) = 1/d + x1/2;
    else
        r1(k) = x1 / (1 - exp(-d*x1));
    end
    if abs(d*x2) < 1e-6
        r2(k) = 1/d + x2/2;
    else
        r2(k) = x2 / (1 - exp(-d*x2));
    end

    % Gating update (Eq 1, 2)
    S1(k) = S1(k-1) + dt*(-S1(k-1)/tau + (1 - S1(k-1))*gamma*r1(k));
    S2(k) = S2(k-1) + dt*(-S2(k-1)/tau + (1 - S2(k-1))*gamma*r2(k));

    % First threshold crossing
    if decision == 0
        if     S1(k) >= thr
            decision    = 1;
            dt_decision = t(k);
        elseif S2(k) >= thr
            decision    = 2;
            dt_decision = t(k);
        end
    end
end
end
