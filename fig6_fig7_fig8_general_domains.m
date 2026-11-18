% fig6_fig7_fig8_general_domains.m
%
% Resonant modes for general smooth resonators (Section 4.2), delta = 1e-5:
%   Figure 6: 5 rings                    (two boundary curves per resonator)
%   Figure 7: 8 rotated ellipses
%   Figure 8: 4 curved double rings      (non-symmetric star-shaped curves)
% Pipeline per case: Theorem 2 seeds (N x N reduction) -> Theorem 1
% effective matrix -> full BEM; Re(u) of branches 1, ceil(N/2), N is plotted.

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'lib'));
delta = 1e-5;

for c = {'rings', 'ellipses', 'curved_rings'}
    [curves, P] = geometry(c{1});
    qd   = make_curves(curves, P.Q);
    Nres = max(qd.owner);

    [C0, K1, K2, M]  = eff_matrices_curves(qd, P.F_eff);
    [K1r, K2r, Pnull] = reduce_to_resonators(qd, C0, K1, K2, M, P.F_eff, Nres);
    [om_th, X] = theorem_seeds(K1r, K2r, delta);
    [om_eff, ~] = effective_solve(C0, K1, K2, M, delta, om_th, Pnull*X);
    [om_bem, V] = bem_solve(@(o) assemble_A_curves(qd, P.F, delta, o), om_eff, P.F + 1);

    fprintf('\n%s: %d/%d BEM branches\n  %-4s %-28s %-10s %-10s\n', c{1}, ...
        nnz(~isnan(om_bem)), Nres, 'br', 'omega_BEM', '|Thm2-BEM|', '|Thm1-BEM|');
    fprintf('  %-4d %+.6e %+.6ei  %.2e   %.2e\n', ...
        [1:Nres; real(om_bem.'); imag(om_bem.'); abs(om_th - om_bem).'; abs(om_eff - om_bem).']);

    br = unique([1, ceil(Nres/2), Nres]);
    qp = make_curves(curves, 256);
    G  = field_grid(qp, P.pad, 200, P.excl);
    D  = qd.N*(2*P.F + 1);
    U  = cell(1, numel(br));
    for k = 1:numel(br)
        u = nan(numel(G.pts), 1, 'like', 1i);
        u(G.int) = layer_field(G.pts(G.int), qp, V(1:D, br(k)),     P.F, om_bem(br(k)));
        u(G.ext) = layer_field(G.pts(G.ext), qp, V(D+1:end, br(k)), P.F, om_bem(br(k)));
        U{k} = normalise_field(u);
    end
    fig = plot_mode_strip(G, U, arrayfun(@(b) sprintf('Branch %d', b), br, 'UniformOutput', false), c{1});
    save_figure(fig, P.name);
end

%% ------------------------------------------------------------------------
function [curves, P] = geometry(name)
% curves: see make_curves.  P: truncations F (BEM) / F_eff (Theorem 1),
% quadrature Q, plot padding and near-boundary exclusion distance.
star = @(r0, kc, ac, ks, bs) struct('r0', r0, 'kc', kc, 'ac', ac, 'ks', ks, 'bs', bs, 'a', [], 'b', []);
switch name
    case 'rings'
        rng(2);  N = 5;
        ctr = random_centers(N, 2.5, 4*sqrt(N));
        for k = N:-1:1
            curves(2*k-1) = curve('star', ctr(k), 0, +1, k, star(1.0, [], [], [], []));
            curves(2*k)   = curve('star', ctr(k), 0, -1, k, star(0.8, [], [], [], []));
        end
        P = struct('F', 8, 'F_eff', 4, 'Q', 128, 'pad', 0.6, 'excl', 0.03, 'name', 'fig6_rings');
    case 'ellipses'
        rng(3);  N = 8;
        ctr = random_centers(N, 0.875, 4*0.35*sqrt(N));
        rot = pi*rand(1, N);
        for k = N:-1:1
            s = star([], [], [], [], []);  s.a = 0.35;  s.b = 0.20;
            curves(k) = curve('ellipse', ctr(k), rot(k), +1, k, s);
        end
        P = struct('F', 8, 'F_eff', 3, 'Q', 128, 'pad', 0.25, 'excl', 0.01, 'name', 'fig7_ellipses');
    case 'curved_rings'
        rng(7);  N = 4;
        ctr = random_centers(N, 0.8, 1.6);
        rot = 2*pi*rand(1, N);
        so = star(0.25,  [5 6 8 10], [0.02 0.01 0.01 0.01], 3, 0.01);
        si = star(0.05,  [2 5 7], [0.005 0.005 0.005],       3, 0.005);
        for k = N:-1:1
            curves(2*k-1) = curve('star', ctr(k), rot(k), +1, k, so);
            curves(2*k)   = curve('star', ctr(k), rot(k), -1, k, si);
        end
        P = struct('F', 28, 'F_eff', 12, 'Q', 256, 'pad', 0.12, 'excl', 0.006, 'name', 'fig8_curved_rings');
end
end

function s = curve(type, ctr, rot, nsign, owner, s)
s.type = type;  s.c = ctr;  s.rot = rot;  s.nsign = nsign;  s.owner = owner;
end
