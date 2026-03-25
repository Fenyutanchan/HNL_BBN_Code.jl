# HNL-BBN Reproduction

Reproducing all derivations and figures in [note](refs/note_from_Xun-Jie.pdf) with a clear, modular, and reusable codebase.

## Stack

- **Julia** with CairoMakie for plotting
- **Jupyter** (julia-vscode) notebooks for interactive exploration
- **DifferentialEquations.jl** for ODE solving

## Shared module

`src/HNL_BBN_Code.jl` — physical constants, cosmology helpers, thermal rates, and HNL physics functions used across notebooks.

## Notebooks

| # | Notebook | PDF Section | Equations | Figure | Status |
|---|----------|-------------|-----------|--------|--------|
| 01 | [thermal-params](notebooks/01-thermal-params.ipynb) | §2.1 | 2.3–2.5 | Fig 1 | [ ] |
| 02 | [gamma-regimes](notebooks/02-gamma-regimes.ipynb) | §2.1 | 2.6–2.8 | Fig 2 | [ ] |
| 03 | [freeze-in](notebooks/03-freeze-in.ipynb) | §2.1.1, App A | 2.9–2.12 | — | [ ] |
| 04 | [freeze-out](notebooks/04-freeze-out.ipynb) | §2.1.2 | 2.13 | — | [ ] |
| 05 | [abundance-contours](notebooks/05-abundance-contours.ipynb) | §2.1.3 | 2.14–2.17 | Fig 3 | [ ] |
| 06 | [nu-gamma-splitting](notebooks/06-nu-gamma-splitting.ipynb) | App B | B.1–B.4 | — | [ ] |
| 07 | [np-conversion](notebooks/07-np-conversion.ipynb) | App C, D | C.1–C.3, D.1–D.6 | — | [ ] |
| 08 | [dark-decay](notebooks/08-dark-decay.ipynb) | §2.2 | TBD | — | [ ] |

## References

- [`refs/code-2507.12270`](refs/code-2507.12270) — code from our previous work [arXiv:2507.12270](https://arxiv.org/abs/2507.12270).

## Setup

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```
