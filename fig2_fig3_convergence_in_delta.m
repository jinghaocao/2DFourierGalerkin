% fig2_fig3_convergence_in_delta.m
%
% Figure 2: error of the asymptotic resonances of Theorem 2 vs delta.
% Figure 3: error of the effective-matrix resonances of Theorem 1 vs delta
%           for F = 0, 3, 6.
% Reference: full BEM (F = 5).  N = 8 random circles of radius 0.3.

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));

delta_list = 10.^(-5:-1:-8);
F_bem      = 5;
F_list     = [0 3 6];
N          = 8;
rng(0);
geom = make_circles(random_centers(N, 1.5, max(8, 3*sqrt(N))), 0.3*ones(1, N));

nd = numel(delta_list);  nF = numel(F_list);
err_th  = nan(N, nd);            % Theorem 2, per branch
err_eff = nan(N, nF, nd);        % Theorem 1, per branch and F

[~, K10, K20] = eff_matrices_circles(geom, 0);
for iF = 1:nF
    [C0{iF}, K1{iF}, K2{iF}] = eff_matrices_circles(geom, F_list(iF)); %#ok<SAGROW>
end

for id = 1:nd
    delta = delta_list(id);
    [om_th, X] = theorem_seeds(K10, K20, delta);
    om_eff = nan(N, nF);
    for iF = 1:nF
        K = 2*F_list(iF) + 1;
        Vs = zeros(N*K, N);  Vs((0:N-1)*K + F_list(iF) + 1, :) = X;   % embed at m = 0
        om_eff(:, iF) = effective_solve(C0{iF}, K1{iF}, K2{iF}, [], delta, om_th, Vs);
    end
    om_bem = bem_solve(@(o) assemble_A_circles(geom, F_bem, delta, o), ...
                       om_eff(:, nF), F_bem + 1, 1e-12);
    err_th(:, id)     = abs(om_th - om_bem);
    err_eff(:, :, id) = abs(om_eff - om_bem);
    fprintf('delta = %.0e : %d/%d BEM branches,  Thm2 log %.2e,  Thm1(F=%d) log %.2e\n', ...
        delta, nnz(~isnan(om_bem)), N, err_th(1, id), F_list(end), err_eff(1, end, id));
end

%% Figure 2 --------------------------------------------------------------
fs = 30;  lw = 3;
ad = abs(log(delta_list));
anchor = @(e, r) r * (e(end) / r(end));          % reference slope through last point

fig2 = figure('Name', 'Figure 2', 'Position', [50 50 1900 760], 'Color', 'w');
ax = subplot(1, 2, 1); hold on;
e = err_th(1, :);
h1 = loglog(delta_list, e, '-ks', 'LineWidth', lw, 'MarkerFaceColor', 'k', 'MarkerSize', 9);
h2 = loglog(delta_list, anchor(e, sqrt(delta_list)./ad), '--r', 'LineWidth', lw + 1);
style_axes(ax, fs, 'Log branch');
legend([h1 h2], {'Log branch', '$O(\sqrt{\delta}/\log\delta)$'}, 'Interpreter', 'latex', 'Location', 'northeast');

ax = subplot(1, 2, 2); hold on;
grey = [0.45 0.45 0.45];
for br = 2:N
    hb = loglog(delta_list, err_th(br, :), '-o', 'Color', grey, 'LineWidth', 2, ...
        'MarkerSize', 6, 'MarkerFaceColor', 'w');
end
e = mean(err_th(2:end, :), 1, 'omitnan');
hm = loglog(delta_list, e, '-ks', 'LineWidth', lw, 'MarkerFaceColor', 'k', 'MarkerSize', 9);
hr = loglog(delta_list, anchor(e, sqrt(delta_list./ad)), '--r', 'LineWidth', lw + 1);
style_axes(ax, fs, 'Regular branches');
legend([hb hm hr], {'Individual branches', 'Mean', '$O(\sqrt{\delta/\log\delta})$'}, ...
    'Interpreter', 'latex', 'Location', 'northeast');
save_figure(fig2, 'fig2_theorem2_error_vs_delta');

%% Figure 3 --------------------------------------------------------------
% F = 3 and F = 6 coincide to plotting accuracy (the Fourier truncation error
% is far below the O(delta^{3/2}) asymptotic remainder), so they are drawn
% with distinct markers, F = 3 on top.
sty = {{[0.6 0.6 0.6], '-',  'o', 9,  [0.6 0.6 0.6]}, ...     % F = 0
       {[0.3 0.3 0.3], '--', 's', 17, 'none'}, ...            % F = 3
       {[0 0 0],       '-',  '.', 26, [0 0 0]}};              % F = 6
order = [1 3 2];
data  = {squeeze(err_eff(1, :, :)), squeeze(mean(err_eff(2:end, :, :), 1, 'omitnan'))};
refs  = {(delta_list./ad).^1.5, delta_list.^1.5 .* ad.^2};
rlab  = {'$O((\delta/|\log\delta|)^{3/2})$', '$O(\delta^{3/2}|\log\delta|^2)$'};
ttl   = {'Log branch', 'Regular branches (mean)'};
fig3 = figure('Name', 'Figure 3', 'Position', [50 50 1900 760], 'Color', 'w');
for p = 1:2
    ax = subplot(1, 2, p); hold on;
    h = gobjects(nF, 1);
    for iF = order
        s = sty{iF};
        h(iF) = loglog(delta_list, data{p}(iF, :), s{2}, 'Color', s{1}, 'LineWidth', lw, ...
            'Marker', s{3}, 'MarkerSize', s{4}, 'MarkerEdgeColor', s{1}, 'MarkerFaceColor', s{5});
    end
    hr = loglog(delta_list, anchor(data{p}(end, :), refs{p}), '--r', 'LineWidth', lw + 1);
    style_axes(ax, fs, ttl{p});
    y = [data{p}(:); hr.YData(:)];  ylim([min(y)/2.5, max(y)*2.5]);
    legend([h; hr], [arrayfun(@(f) sprintf('$F=%d$', f), F_list, 'UniformOutput', false), rlab(p)], ...
        'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', fs - 4);
end
save_figure(fig3, 'fig3_theorem1_error_vs_delta');

function style_axes(ax, fs, ttl)
set(ax, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse', 'FontSize', fs, ...
    'LineWidth', 1.2, 'TickLabelInterpreter', 'latex');
grid(ax, 'on');
xlabel(ax, '$\delta$', 'Interpreter', 'latex');
ylabel(ax, 'Absolute error', 'Interpreter', 'latex');
title(ax, ttl, 'Interpreter', 'latex');
end
