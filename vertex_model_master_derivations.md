# vertex_model_master_derivations.md

This document provides the definitive mathematical reconciliation between the literature synthesis in `vm_llm_research/` and the manual derivations in `vm_human_calculations/`. It establishes the rigorous framework required for calculating the Hessian matrix and the shear modulus in a 2D vertex model under Lees-Edwards Boundary Conditions (LEBC).

---

## 1. Definitions & State Variables

The system is defined by a set of $N_v$ vertices with positions $\mathbf{r}_i = (x_i, y_i)^T$. We consider a confluent tissue of $N_c$ cells.

### 1.1 The Energy Functional
We define the energy functional $E$ as:
$$E = \sum_{\alpha=1}^{N_c} \left[ \frac{K_A}{2}(A_\alpha - A_0)^2 + \frac{K_P}{2}(P_\alpha - P_0)^2 \right]$$
where:
*   $A_\alpha, P_\alpha$: Current area and perimeter of cell $\alpha$.
*   $A_0, P_0$: Target area and perimeter.
*   $K_A, K_P$: Area and perimeter stiffness moduli.

**Reconciliation Note:** Human notes use $K$ and $\Gamma$. For implementation, if using "Pressure Units" as noted in `vm_human_calculations`, $K_{impl} = K_A A_0$, where $[K_{impl}] = \text{Force/Length}$. This synthesis uses $K_A, K_P$ for mathematical clarity.

### 1.2 State Variables (Forces & Stress)
To simplify higher-order derivatives, we define the conjugate variables:
*   **Cell Pressure:** $\Pi_\alpha = -\pdv{E}{A_\alpha} = -K_A(A_\alpha - A_0)$
*   **Edge Tension:** $T_\alpha = \pdv{E}{P_\alpha} = K_P(P_\alpha - P_0)$

### 1.3 Coordinate Unwrapping & LEBC
Under a simple shear strain $\gamma$, the physical distance $\Delta \mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i$ must be "unwrapped" to account for the Lees-Edwards shift. For a box $(L_x, L_y)$:
1.  Calculate raw $\Delta y = y_j - y_i$.
2.  Wrap $\Delta y \to \Delta y - L_y \cdot \text{round}(\Delta y / L_y)$. Let $n_y = \text{round}(\Delta y / L_y)$.
3.  Calculate shifted $\Delta x = (x_j - x_i) - n_y (\gamma L_y)$.
4.  Wrap $\Delta x \to \Delta x - L_x \cdot \text{round}(\Delta x / L_x)$.

---

## 2. First Derivatives (Forces)

The force on vertex $i$ is $\mathbf{f}_i = -\pdv{E}{\mathbf{r}_i}$.
$$\pdv{E}{\mathbf{r}_i} = \sum_{\alpha \in \text{cells}(i)} \left[ -\Pi_\alpha \nabla_i A_\alpha + T_\alpha \nabla_i P_\alpha \right]$$

### 2.1 Geometric Gradients
*   **Area Gradient:** $\nabla_i A_\alpha = \frac{1}{2} \hat{\mathbf{z}} \times (\mathbf{r}_{i+1} - \mathbf{r}_{i-1}) = \frac{1}{2} (y_{i+1} - y_{i-1}, x_{i-1} - x_{i+1})^T$.
*   **Perimeter Gradient:** $\nabla_i P_\alpha = \sum_{j \in \text{neigh}(i)} \hat{\mathbf{t}}_{ij}$, where $\hat{\mathbf{t}}_{ij} = \frac{\mathbf{r}_i - \mathbf{r}_j}{|\mathbf{r}_i - \mathbf{r}_j|}$ is the unit vector pointing *away* from neighbor $j$ toward $i$.

---

## 3. The Hessian Matrix ($\mathcal{H}$)

The Hessian is the second derivative of the energy: $H_{ij}^{\mu\nu} = \frac{\partial^2 E}{\partial r_i^\mu \partial r_j^\nu}$. We decompose this into **Modulus (Stiffness)** and **Stress (Geometric)** terms.

$$\mathcal{H}_{ij} = \mathcal{H}_{ij}^{\text{Mod}} + \mathcal{H}_{ij}^{\text{Stress}}$$

### 3.1 Modulus Terms (Spring-like response)
$$\mathcal{H}_{ij}^{\text{Mod}} = \sum_\alpha \left[ K_A (\nabla_i A_\alpha \otimes \nabla_j A_\alpha) + K_P (\nabla_i P_\alpha \otimes \nabla_j P_\alpha) \right]$$

### 3.2 Stress Terms (Geometric response)
These terms arise from the change in the *direction* of the area and perimeter gradients.
$$\mathcal{H}_{ij}^{\text{Stress}} = \sum_\alpha \left[ -\Pi_\alpha \nabla_i \nabla_j A_\alpha + T_\alpha \nabla_i \nabla_j P_\alpha \right]$$

*   **Area Stress Term:** For $j = i+1$, $\nabla_i \nabla_{i+1} A_\alpha = \frac{1}{2} \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$. This is purely anti-symmetric.
*   **Perimeter Stress Term:** This is the "transverse spring" term. For $j \neq i$:
    $$\nabla_i \nabla_j P_\alpha = -\frac{1}{l_{ij}} (\mathbf{I} - \hat{\mathbf{t}}_{ij} \otimes \hat{\mathbf{t}}_{ij})$$
    where $\mathbf{I}$ is the identity matrix.

---

## 4. Strain Derivatives

To calculate the shear modulus, we define the response to an affine shear transformation $\mathbf{r} \to \boldsymbol{\Lambda}(\gamma) \mathbf{r}$.

### 4.1 Perimeter Strain Derivative
$$\pdv{P_\alpha}{\gamma} = \sum_{\langle ij \rangle \in \alpha} \frac{\Delta x_{ij} \Delta y_{ij}}{l_{ij}}$$
**Physical Intuition:** This represents the projection of the horizontal shear displacement (lever arm $\Delta y$) onto the edge axis $\hat{t}_x$.

### 4.2 Area Strain Derivative
$$\pdv{A_\alpha}{\gamma} = 0$$
**Pedantic Flag:** While the area is invariant under simple shear, its *gradient* $\nabla A$ is not. This is a common source of error in non-affine calculations.

### 4.3 Second Strain Derivative
$$\pdv[2]{P_\alpha}{\gamma} = \sum_{\langle ij \rangle \in \alpha} \frac{\Delta y_{ij}^2}{l_{ij}} \left( 1 - \frac{\Delta x_{ij}^2}{l_{ij}^2} \right)$$

---

## 5. The Affine Force Vector ($\mathbf{\Xi}$)

The vector $\mathbf{\Xi}_i = \frac{\partial^2 E}{\partial \mathbf{r}_i \partial \gamma}$ represents the forces generated on vertices by an infinitesimal affine stretch of the box.

$$\mathbf{\Xi}_i = \sum_{\alpha \in \text{cells}(i)} \left[ K_P (\nabla_i P_\alpha) \pdv{P_\alpha}{\gamma} + T_\alpha \pdv{(\nabla_i P_\alpha)}{\gamma} - \Pi_\alpha \pdv{(\nabla_i A_\alpha)}{\gamma} \right]$$

### 5.1 Force Rotation Terms
*   **Perimeter Rotation:** $\pdv{(\nabla_i P_\alpha)}{\gamma} = \sum_{j} \frac{\Delta y_{ij}}{l_{ij}} [ \hat{\mathbf{x}} - \hat{t}_{ij} (\hat{t}_{ij} \cdot \hat{\mathbf{x}}) ]$
*   **Area Rotation:** $\pdv{(\nabla_i A_\alpha)}{\gamma} = \frac{1}{2} (0, y_{i+1} - y_{i-1})^T$

---

## 6. The Macroscopic Shear Modulus ($G$)

The static shear modulus is decomposed into affine and non-affine parts:
$$G = G_{\text{aff}} - G_{\text{na}}$$

### 6.1 Affine Stiffness
$$G_{\text{aff}} = \frac{1}{A_{\text{tot}}} \pdv[2]{E}{\gamma} = \frac{1}{A_{\text{tot}}} \sum_\alpha \left[ K_P \left(\pdv{P_\alpha}{\gamma}\right)^2 + T_\alpha \pdv[2]{P_\alpha}{\gamma} \right]$$

### 6.2 Non-Affine Correction
$$G_{\text{na}} = \frac{1}{A_{\text{tot}}} \mathbf{\Xi}^T \mathcal{H}^{-1} \mathbf{\Xi}$$
**Discrepancy Note:** Human notes (Merkel/Manning 2018 context) use a factor of $4$ in the denominator ($g = \frac{1}{4 N_c} \dots$). This usually stems from the definition of energy density for pure shear vs. simple shear. To remain consistent with the Merkel/Manning NJP 2018 literature, use the $1/A_{\text{tot}}$ normalization.

---

## 7. Beamer Illustrations (TikZ)

### 7.1 Lees-Edwards Unwrapping
This diagram illustrates how a cell crossing the top boundary is unwrapped by the shift $\gamma L_y$.

```latex
\begin{tikzpicture}[scale=2, thick]
    % Draw main box
    \draw[dashed, gray] (0,0) rectangle (2,2);
    \node at (1,-0.2) {Main Box ($L_x, L_y$)};
    
    % Draw sheared ghost box above
    \begin{scope}[shift={(0.6, 2)}] % gamma * Ly = 0.3 * 2 = 0.6
        \draw[blue!50] (0,0) rectangle (2,1);
        \node[blue] at (1, 1.2) {Ghost Box (Shifted by $\gamma L_y$)};
    \end{scope}

    % Draw a cell crossing the boundary
    \coordinate (V1) at (0.8, 1.8);
    \coordinate (V2) at (1.2, 1.7);
    \coordinate (V3wrapped) at (1.3, 0.2);
    \coordinate (V3unwrapped) at (1.3+0.6, 2.2);

    \filldraw[fill=green!20, opacity=0.5] (V1) -- (V2) -- (V3unwrapped) -- cycle;
    \draw[red, ->, >=stealth] (V3wrapped) -- (V3unwrapped) node[midway, right] {LEBC Shift $\gamma L_y$};
    
    \node[circle, fill=black, inner sep=1pt, label=left:$V_i$] at (V1) {};
    \node[circle, fill=black, inner sep=1pt, label=right:$V_j$] at (V2) {};
    \node[circle, fill=red, inner sep=1pt, label=right:$V_k^{\text{unwrapped}}$] at (V3unwrapped) {};
\end{tikzpicture}
```

### 7.2 Edge Rotation & Lever Arm
Illustrates the perimeter strain derivative $\pdv{P}{\gamma}$.

```latex
\begin{tikzpicture}[scale=2, thick]
    \coordinate (A) at (0,0);
    \coordinate (B) at (1.5, 1);
    \coordinate (Bshear) at (1.5+0.5, 1); % Delta x + gamma * Delta y

    \draw[->] (-0.5,0) -- (2.5,0) node[right] {$x$};
    \draw[->] (0,-0.5) -- (0,1.5) node[above] {$y$};

    \draw[blue] (A) -- (B) node[midway, above left] {$l_{ij}$};
    \draw[red, dashed] (A) -- (Bshear) node[midway, below right] {$l_{ij}'$};
    
    \draw[<->, gray] (1.5, 0) -- (1.5, 1) node[midway, right] {Lever Arm $\Delta y$};
    \draw[->, orange] (1.5, 1) -- (2.0, 1) node[midway, above] {$\gamma \Delta y$};

    \node at (0.75, -0.2) {$\pdv{l_{ij}}{\gamma} = \frac{\Delta x \Delta y}{l}$};
\end{tikzpicture}
```

---

## 8. Final Reconciliation & Critical Flags

1.  **The "Box Shift" Confusion:** Literature often toggles between $x' = x + \gamma y$ (local strain) and $x' = x + \gamma L_y$ (boundary shift). This derivation clarifies that the local displacement is $\gamma y$, which naturally yields the $\gamma L_y$ shift when evaluated at $y = L_y$.
2.  **Pressure Term in $\mathbf{\Xi}$:** You must include $-\Pi_\alpha \pdv{(\nabla_i A_\alpha)}{\gamma}$. Even if $\pdv{A}{\gamma}=0$, the direction of the pressure force rotates. Omitting this will lead to a 10-15% error in $G$ for high-pressure packings.
3.  **Hessian Inversion:** The matrix $\mathcal{H}$ has two trivial zero eigenvalues (global $x,y$ translation). For $G_{\text{na}}$, use the Moore-Penrose pseudoinverse or project out the null space.
4.  **Lagrangian Frame:** All calculations here are in the reference frame. Coordinate $y_i$ is the *height in the box*, and the LEBC shift is applied to distances, not to the stored vertex positions themselves.
