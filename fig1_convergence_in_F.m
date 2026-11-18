% fig1_convergence_in_F.m
%
% Figure 1: exponential convergence of the full Fourier-Galerkin BEM in the
% truncation parameter F (Proposition 1), N = 2 circles, delta = 1e-8.
%   (a) |omega_F - omega_ref| for F = 0..5, reference F_ref = 6
%   (b) Fourier coefficients |phi_n| of the eigenmode on circle 1 (F = 25)

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

geom   = make_circles([0 1], [0.25 0.25]);
delta  = 1e-8;
F_list = 0:5;
F_ref  = 6;
F_coef = 25;

[~, K10, K20] = eff_matrices_circles(geom, 0);
seeds = theorem_seeds(K10, K20, delta);
asm   = @(F) @(o) assemble_A_circles(geom, F, delta, o);

om_ref = bem_solve(asm(F_ref), seeds, F_ref + 1);
err = nan(2, numel(F_list));
for k = 1:numel(F_list)
    err(:, k) = abs(bem_solve(asm(F_list(k)), om_ref, F_list(k) + 1, 1e-6) - om_ref);
end
[~, vec] = bem_solve(asm(F_coef), om_ref, F_coef + 1);

fs = 20;  sty = {'-o', '-.^'};
fig1a = figure('Name', 'Figure 1a', 'Color', 'w');  hold on;
for br = 1:2
    semilogy(F_list, err(br, :), sty{br}, 'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 8);
end
set(gca, 'YScale', 'log', 'FontSize', fs - 4, 'TickLabelInterpreter', 'latex');  grid on;
xlabel('$F$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('Absolute error', 'Interpreter', 'latex', 'FontSize', fs);
save_figure(fig1a, 'fig1a_frequency_error_vs_F');

fig1b = figure('Name', 'Figure 1b', 'Color', 'w');  hold on;
n = -F_coef:F_coef;  sty = {'-', '--'};
for br = 1:2
    c = abs(vec(1:2*F_coef+1, br));               % interior density on circle 1
    semilogy(n, c / max(c), sty{br}, 'Color', 'k', 'LineWidth', 1.5);
end
set(gca, 'YScale', 'log', 'FontSize', fs - 4, 'TickLabelInterpreter', 'latex');  grid on;
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('Eigenmode Fourier coefficients', 'Interpreter', 'latex', 'FontSize', fs);
save_figure(fig1b, 'fig1b_mode_coefficient_decay');
fprintf('omega_ref = %s\n', mat2str(om_ref, 10));
disp(err)
