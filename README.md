# Fourier–Galerkin methods for subwavelength resonances in 2D

MATLAB code reproducing the figures of

> J. Cao, *Fourier–Galerkin Methods for Subwavelength Resonances in
> two-dimensional Acoustic Metamaterials*.

Each script is self-contained, runs from this folder, and writes `.eps` and
`.png` files to `figures/`.

| Script | Paper figure | Content |
|---|---|---|
| `fig1_convergence_in_F.m` | Fig. 1 | Exponential convergence of the full BEM in the truncation `F` (Prop. 1) |
| `fig2_fig3_convergence_in_delta.m` | Figs. 2, 3 | Error vs. δ of Theorem 2 (asymptotic) and Theorem 1 (effective matrix, F = 0, 3, 6) |
| `fig4a_runtime.m` | Fig. 4a | Runtime vs. number of resonators N |
| `fig4b_fig5_circles_N25.m` | Figs. 4b, 5 | N = 25 circles: spectrum, and modes from BEM / Theorem 1 / Theorem 2 |
| `fig6_fig7_fig8_general_domains.m` | Figs. 6, 7, 8 | Modes of rings, ellipses and curved double rings |

Tested with MATLAB R2022a; no toolboxes are required.

## Library (`lib/`)

| Function | Purpose |
|---|---|
| `make_circles`, `random_centers` | Circle geometries |
| `make_curves` | Quadrature on general smooth curves (star-shaped or elliptic; inner boundaries of rings via `nsign = -1`) |
| `assemble_A_circles` | Full BEM matrix `A_F(ω,δ)` (eq. 25) for circles, closed-form multipole entries (§4.1.1) |
| `assemble_A_curves` | Full BEM matrix on general curves, with singularity splitting on the diagonal |
| `eff_matrices_circles` | `C_0^F, K_1^F, K_2^F` (eq. 46) for circles, closed-form diagonal symbols (eq. 72) |
| `eff_matrices_curves` | Same for general curves, Kussmaul–Martensen splitting for `K_2^F` (eqs. 74–76) |
| `reduce_to_resonators` | `N × N` leading-order matrices for resonators with several boundary curves |
| `theorem_seeds` | Asymptotic resonances of Theorem 2 (eqs. 56–57) |
| `effective_solve` | Resonances of the effective matrix `R_F(ω,δ)` of Theorem 1 (eq. 47) |
| `bem_solve` | Newton/Müller refinement on the full BEM (Algorithm 1) |
| `layer_field`, `field_grid`, `normalise_field`, `plot_mode_strip`, `save_figure` | Mode evaluation and plotting |

## Conventions

- Unknowns of the full BEM are ordered `[interior density; exterior density]`,
  with Fourier modes `-F..F` per boundary curve, expanded w.r.t. arclength.
- The effective matrices use the orthonormal basis `e^{imθ}/sqrt|∂D_j|` (eq. 44).
- Random configurations use fixed seeds, so every figure is reproducible.
