% fig4b_fig5_circles_N25.m
%
% N = 25 random circles, delta = 1e-5.
% Figure 4b: resonances in the complex plane (BEM, Theorem 1, Theorem 2).
% Figure 5 : Re(u) of branches 1, 6, 10 from (a) the full BEM, (b) the
%            effective matrix of Theorem 1 and (c) the asymptotic modes of
%            Theorem 2.

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

N = 25;  F = 3;  F_eff = 3;  delta = 1e-5;  a = 0.3;
branches = [1 6 10];
rng(1);
geom = make_circles(random_centers(N, 1.5, 10), a*ones(1, N));

%% resonances ------------------------------------------------------------
[~, K10, K20] = eff_matrices_circles(geom, 0);
[om_th, X] = theorem_seeds(K10, K20, delta);
[C0, K1, K2] = eff_matrices_circles(geom, F_eff);
K = 2*F_eff + 1;
Vs = zeros(N*K, N);  Vs((0:N-1)*K + F_eff + 1, :) = X;
[om_eff, V_eff] = effective_solve(C0, K1, K2, [], delta, om_th, Vs);
[om_bem, V_bem] = bem_solve(@(o) assemble_A_circles(geom, F, delta, o), om_eff, F + 1);
fprintf('BEM branches found: %d/%d\n', nnz(~isnan(om_bem)), N);
fprintf('max |Thm1 - BEM| = %.2e,  max |Thm2 - BEM| = %.2e\n', ...
    max(abs(om_eff - om_bem)), max(abs(om_th - om_bem)));

%% Figure 4b -------------------------------------------------------------
fig4b = figure('Name', 'Figure 4b', 'Position', [50 50 800 600], 'Color', 'w');
h(1) = plot(real(om_bem), imag(om_bem), 'ro', 'MarkerSize', 10, 'LineWidth', 2);  hold on;
h(2) = plot(real(om_eff), imag(om_eff), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
h(3) = plot(real(om_th),  imag(om_th),  'ks', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('$\mathrm{Re}(\omega)$', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('$\mathrm{Im}(\omega)$', 'Interpreter', 'latex', 'FontSize', 30);
legend(h, {'BEM', sprintf('Theorem 1 ($F=%d$)', F_eff), 'Theorem 2'}, ...
    'Interpreter', 'latex', 'Location', 'best', 'FontSize', 24);
set(gca, 'FontSize', 24, 'TickLabelInterpreter', 'latex');  grid on;
save_figure(fig4b, 'fig4b_resonances_N25');

%% Figure 5 --------------------------------------------------------------
circ = struct('type', 'star', 'r0', num2cell(geom.a.'), 'kc', [], 'ac', [], 'ks', [], 'bs', [], ...
    'c', num2cell(geom.c.'), 'rot', 0, 'nsign', 1, 'owner', num2cell(1:N));
qp = make_curves(circ, 256);
G  = field_grid(qp, 0.6, 200, 0.02);
D  = N*(2*F + 1);
sq = sqrt(qp.perim(:));                     % normalised -> arclength basis
U  = cell(3, numel(branches));
for k = 1:numel(branches)
    br = branches(k);
    u = nan(numel(G.pts), 1, 'like', 1i);    % (a) BEM: interior / exterior densities
    u(G.int) = layer_field(G.pts(G.int), qp, V_bem(1:D, br),     F, om_bem(br));
    u(G.ext) = layer_field(G.pts(G.ext), qp, V_bem(D+1:end, br), F, om_bem(br));
    U{1, k} = normalise_field(u);
    c = reshape(V_eff(:, br), [], N) ./ sq.';            % (b) Theorem 1
    u(:) = NaN;  m = G.int | G.ext;
    u(m) = layer_field(G.pts(m), qp, c(:), F_eff, om_eff(br));
    U{2, k} = normalise_field(u, U{1, k});
    u(m) = layer_field(G.pts(m), qp, X(:, br) ./ sq, 0, om_th(br));   % (c) Theorem 2
    U{3, k} = normalise_field(u, U{1, k});
end
lab = {'', sprintf(' (Theorem 1, $F=%d$)', F_eff), ' (Theorem 2)'};
name = {'fig5a_modes_BEM', 'fig5b_modes_theorem1', 'fig5c_modes_theorem2'};
for r = 1:3
    t = arrayfun(@(b) sprintf('Branch %d%s', b, lab{r}), branches, 'UniformOutput', false);
    save_figure(plot_mode_strip(G, U(r, :), t, name{r}), name{r});
end
