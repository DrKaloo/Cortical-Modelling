%% Q7 Jacobian stability analysis
% Kaloyan Todorov
%
% For each fixed point from Q7_phase_plane.m:
%   - build the 2x2 Jacobian of the deterministic flow
%   - compute its eigenvalues
%   - classify: all Re(lambda) < 0 -> stable, mixed signs -> saddle,
%     all Re > 0 -> unstable
%
% F_i = -S_i/tau + (1 - S_i)*gamma*r(x_i),  x_i = a*Isyn_i - b
% r(x)  = x / (1 - exp(-d*x))
% r'(x) = [1 - exp(-d*x) - x*d*exp(-d*x)] / [1 - exp(-d*x)]^2

clc, clearvars, close all

% For marking and loading purposes (check)
if ~exist('Q7_results.mat', 'file')
    error('Q7_results.mat not found. Run Q7_phase_plane.m first.')
end


load('Q7_results.mat', 'fp', 'I1', 'I2')

tau   = 0.100;
gamma = 0.641;
a     = 270e9;
b     = 108;
d     = 0.15;
W11   = 0.2609e-9;
W22   = 0.2609e-9;
W12   = 0.0497e-9;
W21   = 0.0497e-9;

% Response function and its analytical derivative
r  = @(x) (abs(d*x) < 1e-6) .* (1/d + x/2) + ...
          (abs(d*x) >= 1e-6) .* ( x ./ (1 - exp(-d*x)) );

rp = @(x) (abs(d*x) < 1e-6) .* 0.5 + ...
          (abs(d*x) >= 1e-6) .* ...
          ( (1 - exp(-d*x) - x.*d.*exp(-d*x)) ./ (1 - exp(-d*x)).^2 );

fprintf('\n Q7: Jacobian analysis \n')
fprintf('Fixed points sorted by S1 (ascending):\n\n')

for k = 1:size(fp, 1)
    S1 = fp(k, 1);
    S2 = fp(k, 2);
    J  = build_jacobian(S1, S2, tau, gamma, a, b, ...
                        W11, W22, W12, W21, I1, I2, r, rp);
    lambda = eig(J);

    fprintf('FP%d  at  (S1, S2) = (%.4f, %.4f)\n', k, S1, S2)
    fprintf('  Jacobian:\n')
    fprintf('    [ %+9.3f   %+9.3f ]\n',   J(1,1), J(1,2))
    fprintf('    [ %+9.3f   %+9.3f ]\n\n', J(2,1), J(2,2))

    re = real(lambda);
    im = imag(lambda);
    if all(abs(im) < 1e-9)
        fprintf('  Eigenvalues:  %+.3f,  %+.3f   (real)\n', re(1), re(2))
    else
        fprintf('  Eigenvalues:  %+.3f %+.3fi,  %+.3f %+.3fi\n', ...
                re(1), im(1), re(2), im(2))
    end

    if all(re < 0)
        cls = 'Stable (attractor)';
    elseif all(re > 0)
        cls = 'Unstable (repeller)';
    else
        cls = 'Saddle';
    end
    fprintf('  Classification: %s\n\n', cls)
end


% 2x2 analytical Jacobian at (S1, S2)
function J = build_jacobian(S1, S2, tau, gamma, a, b, ...
                            W11, W22, W12, W21, I1, I2, r, rp)

    Isyn1 = W11*S1 - W12*S2 + I1;
    Isyn2 = W22*S2 - W21*S1 + I2;
    x1 = a*Isyn1 - b;
    x2 = a*Isyn2 - b;

    % d(x_i)/d(S_j)
    dx1_dS1 =  a * W11;    dx1_dS2 = -a * W12;
    dx2_dS1 = -a * W21;    dx2_dS2 =  a * W22;

    % Chain rule
    dr1_dS1 = rp(x1) * dx1_dS1;
    dr1_dS2 = rp(x1) * dx1_dS2;
    dr2_dS1 = rp(x2) * dx2_dS1;
    dr2_dS2 = rp(x2) * dx2_dS2;

    r1 = r(x1);
    r2 = r(x2);

    % J_ij = d(F_i)/d(S_j)
    J11 = -1/tau - gamma*r1 + (1 - S1)*gamma*dr1_dS1;
    J12 =                     (1 - S1)*gamma*dr1_dS2;
    J21 =                     (1 - S2)*gamma*dr2_dS1;
    J22 = -1/tau - gamma*r2 + (1 - S2)*gamma*dr2_dS2;

    J = [J11 J12; J21 J22];
end
