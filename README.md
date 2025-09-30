# Richardson–BCS Solver & Correlators (Python)

A robust Python implementation of the **reduced BCS (Richardson) model** that:

* solves the **Bethe–Richardson equations** for the rapidities $E_\alpha$ (exact, non-mean-field),
* constructs the **Gaudin (Jacobian) matrix**,
* and computes **pair–pair correlators**  $\langle A^\dagger_{\mathbf{k}} A_{\mathbf{k}'} \rangle$
  on square lattices up to at least (L=20) with strong numerical safeguards.

This repository is designed to be **plug-and-play** for generating trustworthy benchmarks and datasets for many-body physics.

---

## Intuition & Background (for newcomers)

### Model

We study the reduced BCS Hamiltonian with time-reversed pairing:

$$
H=\sum_{\mathbf{k},s}\varepsilon_{\mathbf{k}},c^\dagger_{\mathbf{k}s}c_{\mathbf{k}s}
;+; g\sum_{\mathbf{k}}A^\dagger_{\mathbf{k}}\sum_{\mathbf{k}'}A_{\mathbf{k}'},,
\qquad
A^\dagger_{\mathbf{k}}=c^\dagger_{\mathbf{k}\uparrow}c^\dagger_{-\mathbf{k}\downarrow}.
$$

Here $\varepsilon_{\mathbf{k}}$ is the single-particle dispersion (default: $\varepsilon_{\mathbf{k}}=\cos k_x+\cos k_y$ on an $L\times L$ grid), and (g) is the attractive pairing strength.

### Exact (non-BCS) eigenstates

Unlike mean-field BCS, this model is exactly solvable. The (2M)-particle eigenstates take a **Richardson (Bethe) ansatz**:

$$
\big|\Psi_{2M}\big\rangle=\prod_{\alpha=1}^{M} B^\dagger_\alpha\big|0\rangle,
\qquad
B^\dagger_\alpha=\sum_{\mathbf{k}} \frac{A^\dagger_{\mathbf{k}}}{2\varepsilon_{\mathbf{k}}+2g-E_\alpha},,
$$

where the complex parameters $E_\alpha$ (called **rapidities**) satisfy a coupled nonlinear system.

### Richardson equations (key equation #1)

For each $\alpha=1,\dots,M$,

$$
\boxed{
\frac{1}{g};+; \sum_{\mathbf{k}} \frac{1}{,2\varepsilon_{\mathbf{k}}+2g-E_\alpha,}
;-; \sum_{\beta\neq\alpha} \frac{2}{,E_\beta-E_\alpha,}
;=; 0
}
\tag{R}
$$

### Gaudin matrix (key equation #2)

Define the **Gaudin (Jacobian) matrix** $G\in\mathbb{C}^{M\times M}$:

$$
\boxed{
\begin{aligned}
G_{\alpha\alpha}&=\sum_{\mathbf{k}}\frac{1}{\big(2\varepsilon_{\mathbf{k}}+2g-E_\alpha\big)^2}
;-;\sum_{\beta\neq\alpha}\frac{2}{\big(E_\beta-E_\alpha\big)^2},\
G_{\alpha\mu}&=\frac{2}{\big(E_\mu-E_\alpha\big)^2}\quad(\mu\neq\alpha).
\end{aligned}}
\tag{G}
$$

### Pair–pair correlator (key equation #3)

Let $\varphi_\alpha(\mathbf{k})=\big(2\varepsilon_{\mathbf{k}}+2g-E_\alpha\big)^{-1}$ and
$\boldsymbol{\varepsilon}(\mathbf{k})=[\varphi_1(\mathbf{k}),\dots,\varphi_M(\mathbf{k})]^T$.
Then the normalized pair–pair correlator is

$$
\boxed{
\frac{\langle \Psi_{2M}\lvert A^\dagger_{\mathbf{k}}A_{\mathbf{k}'}\rvert \Psi_{2M}\rangle}
{\langle \Psi_{2M}\vert \Psi_{2M}\rangle}
;=;
\boldsymbol{\varepsilon}(\mathbf{k})^{!T},G^{-1},\boldsymbol{\varepsilon}(\mathbf{k}')
}
\tag{C}
$$

This is what the code evaluates on both the representative ${\mathbf{k},-\mathbf{k}}$ set and the full grid.

---

## What this code implements

* **Lattice & symmetry:** builds the (k)-mesh, tracks time-reversed pairs ${\mathbf{k},-\mathbf{k}}$, handles TRIM points.
* **Root solver for (R):**

  * **Multiplicity-aware seeding:** puts initial roots across the distinct gaps of $A_{\mathbf{k}}=2\varepsilon_{\mathbf{k}}+2g$, allowing multiple seeds per gap when degeneracies occur.
  * **Complex-safe residuals/Jacobians:** a tiny imaginary regulator $+i,\eta$ prevents divisions by (near) zero.
  * **Adaptive continuation in (g):** start from small (g), step to $g_{\rm target}$ with automatic step shrinking in hard regions; quick **polish** at each step.
  * **Damped Newton:** Tikhonov ladder $(G+\lambda I)\delta=-R$and pseudoinverse fallback stabilize updates when (G) is ill-conditioned.
* **Correlators:** builds $\Phi$ with entries $\varphi_\alpha(\mathbf{k})$, forms $C=\Phi^{!T}G^{-1}\Phi$, and returns the **real-Hermitian** version $\tfrac12(C+C^\dagger)$.

---

## Install

```bash
git clone https:https://github.com/LexuZhao/richardson-bcs-exact/tree/main.git

# optional virtualenv
python -m venv .venv && source .venv/bin/activate

pip install -r requirements.txt
# or simply: pip install numpy
```

---

## Quick start

```python
from richardson_bcs import correlators_for_L

# Square lattice L x L, pairing density alpha = M/L^2, coupling g
out = correlators_for_L(L=12, alpha=0.25, g=1.0, verbose=True, save_npz=True)

E = out["E_roots"]   # Richardson rapidities (complex)
G = out["G"]         # Gaudin matrix (complex)
C = out["C_full"]    # real-Hermitian correlator on full k-grid (L^2 x L^2)

# sanity diagnostic (weighted diagonal sum on full grid)
print(out["M_from_full"])
```

---

## API (most-used functions)

### `correlators_for_L(L, alpha=0.25, g=1.0, verbose=True, save_npz=True) -> dict`

Builds the grid, solves (R) via continuation, constructs (G) and correlators. Returns a dictionary with:

* `E_roots` ((M,)): rapidities $E_\alpha$
* `G` ((M,M)): Gaudin matrix
* `C_reps` $(\Omega,\Omega)$: correlator on representatives ${\mathbf{k},-\mathbf{k}}$
* `C_full` $(L^2,L^2)$: correlator on full grid
* plus `k_reps`, `k_full`, `weights_full`, `M_from_full`, etc.

### `continuation_solve_richardson(eps, g_target, M, ...) -> np.ndarray`

Adaptive homotopy in (g) with:

* **ETA plateau → cool-down** (big $\eta$ early for stability; small $\eta$ near the end for accuracy),
* **tiered polish** at each (g),
* **step growth/bisection** based on residual and conditioning.

### `solve_richardson_once(eps, g, M, E_init=None, tol=1e-10, maxiter=300, verbose=False)`

One Newton solve at fixed (g) with:

* complex-safe residuals,
* Tikhonov ladder and `pinv` fallback,
* pole-aware backtracking and per-iteration step capping.

---

## Numerical knobs (how to tune)

Defined at the top of the file:

* `ETA` — tiny imaginary regulator $ \eta $:
  keeps denominators like $2\varepsilon+2g-E_\alpha$ and $E_\beta-E_\alpha$ off real poles.
  **Raise** (e.g., $10^{-7})–(10^{-6})$ for big (L)/tight gaps; **lower** (e.g., $10^{-12}$) for final polishing.

* `PINV_RCOND` — SVD cutoff for pseudoinverse (G^\dagger):
  **larger** = stronger damping when (G) is ill-conditioned (safer, possibly slower).

* `STEP_LIMIT_FAC` — caps the Newton step as a fraction of a typical (A)-gap:
  **smaller** to prevent overshoot near poles (0.1–0.3), **larger** to speed up when smooth (0.6–0.9).

* `BT_MAX` — backtracking tries per Newton iteration:
  increase if steps are often rejected.

**Typical symptoms → quick fixes**

* `||R||_∞` stalls + `cond(G)` huge: ↑ `ETA`, ↑ `PINV_RCOND`, ↓ `STEP_LIMIT_FAC`, ↑ continuation steps.
* Good at small (g), bad later: keep a **larger ETA** for longer fraction of the path; do **polish tiers** at each (g).
* Exploding diagnostics in (C): switch earlier to `pinv` and use a larger `rcond`.

---

## Reproducible equations used in code

* Richardson equations (R) for rapidities $E_\alpha$.
* Gaudin matrix (G) as the Jacobian of (R).
* Correlator formula (C) $C_{\mathbf{k},\mathbf{k}'}=\Phi^{!T}G^{-1}\Phi$ with $\Phi_{\alpha\mathbf{k}}=(2\varepsilon_\mathbf{k}+2g-E_\alpha)^{-1}$.

These three boxed formulas are exactly what the implementation evaluates.

---

## Performance & scaling tips

* Cost is dominated by solving linear systems with (G): roughly (O(M^3)) per Newton step.
* For larger (L) (dense/degenerate spectra):

  * increase continuation `n_steps`,
  * keep the ETA **plateau** longer before cooling,
  * tighten step limiter (`STEP_LIMIT_FAC`),
  * and prefer `pinv` when `cond(G)` is large.

After convergence, reduce `ETA` and do a **final polish** at target (g) for best accuracy.

---

## Repository layout

```
richardson_bcs.py      # main module (solver + correlators)
examples/
  quickstart.ipynb     # end-to-end run and plots
LICENSE
README.md
```

---

## How to cite

If this code helps your work, please cite the repository:

```bibtex
@misc{LexuZhao_RichardsonBCS_2025,
  author       = {Lexu Zhao},
  title        = {Richardson–BCS Solver \& Correlators (Python)},
  year         = {2025},
  howpublished = {\url{https://github.com/LexuZhao/richardson-bcs-exact/tree/main}},
  note         = {Exact Richardson roots, Gaudin matrix, and pair--pair correlators with robust continuation}
}
```
