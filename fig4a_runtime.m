% fig4a_runtime.m
%
% Figure 4a: wall-clock time to compute all N resonances of N random circles
% with (i) the asymptotic formulas of Theorem 2, (ii) the effective matrix of
% Theorem 1 (F = 3) and (iii) the full BEM (F = 5) seeded by (ii).
% Set N_list = 2.^(1:9) and N_bem_max = 16 to probe large N (Theorem 1 cost
% is dominated by the O((N(2F+1))^3) generalised eigensolve).

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

N_list    = 2.^(1:5);
N_bem_max = Inf;
F_bem = 5;  F_eff = 3;  delta = 1e-5;
rng(0);

t = nan(3, numel(N_list));
for k = 1:numel(N_list)
    N = N_list(k);
    geom = make_circles(random_centers(N, 1.5, max(8, 3*sqrt(N))), 0.3*ones(1, N));

    tic;                                                    % Theorem 2
    [~, K10, K20] = eff_matrices_circles(geom, 0);
    om_th = theorem_seeds(K10, K20, delta);
    t(1, k) = toc;

    tic;                                                    % Theorem 1
    [C0, K1, K2] = eff_matrices_circles(geom, F_eff);
    om_eff = effective_solve(C0, K1, K2, [], delta, om_th, [], 'cluster');
    t(2, k) = toc;

    if N <= N_bem_max                                       % full BEM
        tic;
        om_bem = bem_solve(@(o) assemble_A_circles(geom, F_bem, delta, o), om_eff, F_bem + 1);
        t(3, k) = toc;
        fprintf('N = %4d: %8.3f %8.3f %8.3f s   (%d/%d BEM, max|eff-bem| = %.1e)\n', N, t(:, k), ...
            nnz(~isnan(om_bem)), N, max(abs(om_eff - om_bem)));
    else
        fprintf('N = %4d: %8.3f %8.3f      -   s\n', N, t(1:2, k));
    end
end

fs = 30;  lw = 3;
fig = figure('Name', 'Figure 4a', 'Position', [50 50 1000 800], 'Color', 'w');
h(1) = loglog(N_list, t(1, :), '-ks',  'LineWidth', lw, 'MarkerSize', 10, 'MarkerFaceColor', 'k');  hold on;
h(2) = loglog(N_list, t(2, :), '--ko', 'LineWidth', lw, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
h(3) = loglog(N_list, t(3, :), '-.k^', 'LineWidth', lw, 'MarkerSize', 12);
set(gca, 'FontSize', fs, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex', 'XTick', N_list);
grid on;  xlim(N_list([1 end]));
xlabel('$N$ (number of circles)', 'Interpreter', 'latex');
ylabel('Runtime (in seconds)', 'Interpreter', 'latex');
legend(h, {'Theorem 2', sprintf('Theorem 1, $F=%d$', F_eff), 'Full BEM'}, ...
    'Interpreter', 'latex', 'Location', 'northwest');
save_figure(fig, 'fig4a_runtime');
