# Richardson–BCS Solver & Correlators (Python)

A robust Python implementation of the **reduced BCS (Richardson) model** that
(i) solves the **Bethe–Richardson equations** for the rapidities \(E_\alpha\) (exact, non-mean-field),
(ii) constructs the **Gaudin (Jacobian) matrix**, and
(iii) computes **pair–pair correlators** \(\langle A^\dagger_{\mathbf{k}} A_{\mathbf{k}'} \rangle\)
on square lattices (tested up to \(L=20\)) with strong numerical safeguards.

---

## Intuition & Background (for newcomers)

**Model.** We study the reduced BCS Hamiltonian with time-reversed pairing
$$
H=\sum_{\mathbf{k},s}\varepsilon_{\mathbf{k}}\,c^\dagger_{\mathbf{k}s}c_{\mathbf{k}s}
\;+\; g\sum_{\mathbf{k}}A^\dagger_{\mathbf{k}}\sum_{\mathbf{k}'}A_{\mathbf{k}'},
\qquad
A^\dagger_{\mathbf{k}}=c^\dagger_{\mathbf{k}\uparrow}c^\dagger_{-\mathbf{k}\downarrow}.
$$
Here \(\varepsilon_{\mathbf{k}}\) is the single-particle dispersion (default
\(\varepsilon_{\mathbf{k}}=\cos k_x+\cos k_y\) on an \(L\times L\) grid), and \(g\) is the attractive pairing strength.

**Exact eigenstates.** Unlike mean-field BCS, this model is **exactly solvable**. The \(2M\)-particle eigenstates take a Richardson (Bethe) ansatz
$$
\big|\Psi_{2M}\big\rangle=\prod_{\alpha=1}^{M} B^\dagger_\alpha\big|0\rangle,
\qquad
B^\dagger_\alpha=\sum_{\mathbf{k}} \frac{A^\dagger_{\mathbf{k}}}{\,2\varepsilon_{\mathbf{k}}+2g-E_\alpha\,},
$$
where complex parameters \(E_\alpha\) (the **rapidities**) solve the coupled non-linear **Richardson equations**.

**Richardson equations (key equation #1).** For each \(\alpha=1,\dots,M\),
$$
\boxed{
\frac{1}{g}
+ \sum_{\mathbf{k}} \frac{1}{\,2\varepsilon_{\mathbf{k}}+2g-E_\alpha\,}
- \sum_{\beta\neq\alpha} \frac{2}{\,E_\beta-E_\alpha\,}
= 0
}
\tag{R}
$$

**Gaudin matrix (key equation #2).** Define \(G\in\mathbb{C}^{M\times M}\):
$$
\boxed{
\begin{aligned}
G_{\alpha\alpha}&=\sum_{\mathbf{k}}\frac{1}{\big(2\varepsilon_{\mathbf{k}}+2g-E_\alpha\big)^2}
-\sum_{\beta\neq\alpha}\frac{2}{\big(E_\beta-E_\alpha\big)^2},\\
G_{\alpha\mu}&=\frac{2}{\big(E_\mu-E_\alpha\big)^2}\qquad(\mu\neq\alpha).
\end{aligned}}
\tag{G}
$$

**Pair–pair correlator (key equation #3).** Let \(\varphi_\alpha(\mathbf{k})=(2\varepsilon_{\mathbf{k}}+2g-E_\alpha)^{-1}\) and
\(\boldsymbol{\varepsilon}(\mathbf{k})=[\varphi_1(\mathbf{k}),\dots,\varphi_M(\mathbf{k})]^T\). Then
$$
\boxed{
\frac{\langle \Psi_{2M}\lvert A^\dagger_{\mathbf{k}}A_{\mathbf{k}'}\rvert \Psi_{2M}\rangle}
{\langle \Psi_{2M}\vert \Psi_{2M}\rangle}
=\boldsymbol{\varepsilon}(\mathbf{k})^{T}G^{-1}\boldsymbol{\varepsilon}(\mathbf{k}')
}
\tag{C}
$$
which the code evaluates on both the representative set \(\{\mathbf{k},-\mathbf{k}\}\) and the full grid.

---

## What this code implements

- **Lattice & symmetry:** builds the \(k\)-mesh, groups \(\{\mathbf{k},-\mathbf{k}\}\), and handles TRIM points.
- **Root solver for (R):**
  - *Multiplicity-aware seeding:* initial roots distributed across distinct gaps of \(A_{\mathbf{k}}=2\varepsilon_{\mathbf{k}}+2g\), allowing multiple seeds per gap when degeneracies occur.
  - *Complex-safe residuals/Jacobians:* a tiny imaginary regulator \(+i\eta\) prevents division by (near) zero.
  - *Adaptive continuation in \(g\):* start at small \(g\), step to \(g_{\rm target}\) with automatic bisection in hard regions plus quick **polish** at each step.
  - *Damped Newton:* Tikhonov ladder \((G+\lambda I)\delta=-R\) and `pinv` fallback stabilize updates when \(G\) is ill-conditioned.
- **Correlators:** build \(\Phi\) with entries \(\varphi_\alpha(\mathbf{k})\), form \(C=\Phi^{T}G^{-1}\Phi\), and return the **real-Hermitian** version \(\tfrac12(C+C^\dagger)\).

---

## Install

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
python -m venv .venv && source .venv/bin/activate  # optional
pip install -r requirements.txt  # or: pip install numpy
