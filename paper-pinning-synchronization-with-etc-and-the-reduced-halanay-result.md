***Pinning Synchronization With Alternative Block Form Halanay Result***

---

Closed form conditions (as opposed to an LMI that needs a solver)

## Problem Formulation
Consider the synchronization network:

```math
\dot X(t) =F(t,X(t)) + c (A \otimes \Gamma) X(t-\tau(t)) + U(t_k), \quad t \ge 0  \quad (1)
```

- $`X(t) \in \mathbb{R}^{Nn}`$ represents the state of the network
- $`F: \mathbb R_+ \times \mathbb R^{nN} \to\mathbb R^{nN}`$ such that $`F := [f^\top, f^\top ... \text{N times}]^\top`$ where $`f: \mathbb{R_+ \times R^n \to R^n}`$ is some nonlinear function
- $`c>0`$ denotes the global coupling strenght
- $`A\in\mathbb R^{N \times N}`$ and $`\Gamma \in \mathbb R^{n\times n}`$ are the outer and inner coupling matrices respectivley
- $`A`$'s diagonal is defined such that for any row $`i`$, $`\sum_{j = 1}^N a_{ij} = 0`$
- $`\tau(t)`$ is the bounded time varying continuous delay such that $`0 \le \tau(t) \le \tau_0`$
- The network initial condition is given by $`X(s) = \Phi(s)`$, $`s \in [-\tau_0, 0]`$ on which $`\Phi`$ is continuous
- Let $`Z \in \mathbb{R}^{Nn}`$ such that $`Z:= [z^\top, z^{\top} \text{ ... N times}]^\top`$ where $`z(t) \in \mathbb R^n`$ is the solution of $`\dot z = f(t,z(t))`$.

Let order the vectors such that the first $`l`$ nodes be the ones that are controlled and the remaining $`N - l`$ are uncontrolled. The system can then be written nodewise as:

```math
\begin{aligned}
&\dot x_i(t) = f(t,x_i(t)) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau(t)) + u_i(t_k), \quad i= 1, ... l \\
&\dot x_i(t) = f(t,x_i(t)) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau(t)), \quad i= l+1... N
\end{aligned}\quad (2)
```

With $`x_i \in \mathbb{R}^{n}`$, $`u_i(t_k) : \mathbb R^n \times \mathbb R\to \mathbb {R}^n`$, with the desired state of this node being $`z(t)`$.  

Define the synchronization error as $`e_i = x_i - z`$, with $`e = [e_1^\top, e_2^\top...e_N^\top]^\top`$, $`\hat e = [e_1^\top, e_2^\top...e_l^\top ]^\top`$ and $`\bar e = [e_{l+1}^\top... e_N^\top]^{\top}`$. We can rewrite the system in terms of its error$`^{[Ap.1]}`$, with $`\bar f_i(t) := f(t,x_i(t)) - f(t,z(t))`$.

```math
\begin{aligned}
&\dot e_i(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau(t)) + u_i(t_k) , \quad i=1, ... l \\
&\dot e_i(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau(t)), \quad i=l+1... N
\end{aligned} \quad (2)
```

Defining the control function:

```math
u_i(t_k) = -d_ie_i(t_k), \quad i=1,...l
```

Let $`\bar d:= \min\{d_1, ... d_l\}`$. Define engineering parameters $`\alpha, \beta >0`$ and $`\sigma\ge 0`$. 

Define the event triggered control as follows:

```math
\begin{aligned}
&\text{ETC Sampling Error:} \quad \epsilon(t) := \hat e(t) - \hat e(t_k), \quad t \in [t_{k}, t_{k+1}), \quad \text{for } k\in \mathbb{N}  \\
& \text{Error Measure:} \quad \mathcal C (e(t), \epsilon(t)) := \bar d \|\epsilon(t)\| \| \hat e(t)\|  \\
& \text{Event Threshold:} \quad \zeta := \alpha \exp({-\beta t}) - \sigma\|\hat e(t)\|^2  \\
& \text{Event Subsequence:} \quad t_{k+1} := \inf\{t \gt t_{k} : \mathcal C(x(t),\epsilon(t)) \ge \zeta(t)\}
\end{aligned}
```

We can thus write our control as:

```math
u(t_k) = d_i \epsilon_i(t) - d_i  e_i(t), \quad i = 1,...l
```

Our complete system model is:

```math
\begin{aligned}
&\dot e_i(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau(t)) + d_i \epsilon_i(t) - d_i  e_i(t), \quad i=1, ... l \\
&\dot e_i(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau(t)), \quad i=l+1... N
\end{aligned} \quad (3)
```

### Assumption 1:
Assume that there exists a constant $`\gamma`$ such that: 

```math
(x-y)^\top [f(t,x) - f(t,y)] \le \gamma \|x-y\|^2 \quad (As.1)
```

For any $`x,y \in \mathbb{R}^n`$ and $`t\ge0`$

### Assumption 2:
Assume that there exists a nonnegative constant $`L`$ such that:

```math
\|f(t,x) - f(t,y) \| \le  L\|x-y\| \quad (As.2)
```

For any $`x,y \in \mathbb{R}^n`$ and $`t\ge0`$

### Lemma: A Halanay Inequality

---

For $`a_i,b_i,c_i > 0`$ and if 

```math
a_1 > b_1 + c_1, \quad a_2 > b_2 + c_2
```

And if the continuous function $`V_i : [\tau_0, \infty) \to \mathbb{R}_+`$ satisfies:

```math
\begin{aligned}
&\dot V_1(t) \le -a_1V(t) + b_1\|V_1\|_\tau + c_1\|V_2\|_\tau + \alpha e^{-\beta t} \\
&\dot V_2(t) \le -a_2V(t) + b_2\|V_2\|_\tau + c_2\|V_1\|_\tau \\
\end{aligned} \quad \forall t \ge t_0
```

then: $`V_i(t) \to 0`$ as $`t \to \infty`$.

### Theorem

---

Let $`\bar d = \min\{d_1,...d_l\}`$, let $`\Omega_{ij} := cA_{ij} \otimes \Gamma`$ and $`\omega_{ij} := \|\Omega_{ij}\|`$, and let $`\mu:= \lambda_{\max}(\Omega^s_{22})`$. If the following are satisfied:

```math
\begin{aligned}
& \text{Uncontrolled Dissipation:} \quad - \mu >  \gamma + \omega_{21} \\ \\
& \text{Maximum Delay Threshold:} \quad \tau_0 < \frac{- \mu -\gamma - \omega_{21}}{\omega_{22}(L + \omega_{22}+ \omega_{21})} \\ \\
&\text{Control Strength:} \quad \bar d > \sigma + \gamma + \omega_{11} + \omega_{21}
\end{aligned}
```

Then network $`(1)`$ realizes exponential synchronization towards the target state $`z(t)`$.

#### Proof

---

The separation into a controlled and uncontrolled subsets of nodes imposes the following partition on the matrix $`A`$:

```math
A = \begin{bmatrix}  A_{11} & A_{12} \\ A_{21} & A_{22}  \end{bmatrix}, \quad A_{11} \in \mathbb R^{l \times l}
```

Let $`\Omega_{ij} := cA_{ij} \otimes \Gamma`$, it follows from equation $`(2)`$ that:$`^{[Ap.2]}`$

```math
\begin{aligned}
&\dot{\hat e}(t) = \hat F(t) + \Omega_{11}\hat e(t-\tau(t)) + \Omega_{12}\bar e(t-\tau(t)) -  D \hat e_i(t) -D \epsilon_i(t) \\
&\dot{\bar e}(t) = \hat F(t) + \Omega_{21}\hat e(t-\tau(t)) + \Omega_{22}\bar e(t-\tau(t)) 
\end{aligned}
```

Where $`\hat F = [\bar f_1^\top, ... \bar f_l^\top]^\top`$, $`\bar F= [\bar f_{l+1}^\top,... \bar f_N^\top]^\top`$, and $`D := \text{diag}\{d_1,..., d_l\} \otimes I_n`$.

Consider two Lyapunov functions:

```math
V_1(t) = \hat e^\top(t) \hat e(t), \quad V_2(t) = \bar e^\top(t) \bar e(t)
```

#### Proof - Part 1

---

By $`(As.1)`$ one has$`^{[Ap.3]}`$: $`\hat e^\top(t) \hat F(t) \le \gamma \|\hat e(t)\|^2`$, then differentiating $`V_1`$ gives:

```math
\begin{aligned}
&\dot V_1(t) = \dot {\hat e}^\top(t) {\hat e}(t) + \hat e^\top(t) {\dot e}(t)  \\
&\dot V_1(t) = 2\hat e^\top(t)\dot {\hat e}(t)  \\
&\dot V_1(t) = 2\hat e^\top(t) \hat F(t) + 2\hat e^\top(t)\Omega_{11}\hat e(t-\tau(t)) + 2\hat e^\top(t)\Omega_{12}\bar e(t-\tau(t)) - 2\hat e^\top(t)D\hat e(t) + 2 \hat e^\top D \epsilon(t)  \\
&\dot V_1(t) \le 2\gamma \|\hat e(t)\|^2 + 2\hat e^\top(t)\Omega_{11}\hat e(t-\tau(t)) + 2\hat e^\top(t)\Omega_{12}\bar e(t-\tau(t)) - 2\hat e^\top(t)D\hat e(t) +  2 \hat e^\top D \epsilon(t) \\
\end{aligned}
```

Note that $`\hat e^\top(t) D \hat e(t) \le \bar d \|\hat e(t)\|^2$$^{[Ap.4]}`$, $`\hat e^\top(t) D \epsilon(t) \le \bar d \|\epsilon(t)\| \| \hat e(t)\| \le \alpha e^{-\beta t} - \sigma\|\hat e(t)\|^2`$ so: 
```math
\begin{aligned}
&\dot V_1(t) \le 2(\gamma - \bar d) \|\hat e(t)\|^2 + 2\hat e^\top(t)\Omega_{11}\hat e(t-\tau(t)) + 2\hat e^\top(t)\Omega_{12}\bar e(t-\tau(t)) + 2\alpha e^{-\beta t} - 2 \sigma \|\hat e(t) \|^2  \\
\end{aligned}
```
Now by lemma 1:
```math
\begin{aligned}
&2\hat e^\top(t) \Omega_{11} \hat e(t - \tau(t)) \le \frac{\| \Omega_{11} \|^2}{\rho_1} \| \hat e(t) \|^2 + \rho_1 \| \hat e\|^2_{\tau_0} \\
&2\hat e^\top(t) \Omega_{12} \hat e(t - \tau(t)) \le \frac{\| \Omega_{12} \|^2}{\rho_2} \| \hat e(t) \|^2 + \rho_2 \| \bar e\|^2_{\tau_0} \\
\end{aligned}
```
Then we can say:
```math
\dot V_1(t) \le -a_1V(t) + b_1\|V_1\|_\tau + c_1\|V_2\|_\tau + \alpha e^{-\beta t}
```
Where:
```math
a_1 := 2(\bar d + \sigma - \gamma ), \quad b_1 := \rho_1 + \frac{\|\Omega_{11}\|^2}{\rho_1} + \frac{\|\Omega_{12}\|^2}{\rho_2}, \quad c_1 :=\rho_2, \quad \alpha' := 2 \alpha
```
#### Proof - Part 2
---
*Obtaining the Halanay Inequality for $`V_2`$*
Differentiating $`V_2(t)`$ gives:
```math
\dot V_2(t) \le 2 \gamma \| \bar e \|^2 + 2 \bar e^\top(t) \Omega_{21} \hat e(t-\tau(t)) + 2 \bar e^\top(t)\Omega_{22} \bar e (t - \tau(t)) \quad (i)
```
Youngs again gives:
```math
2 \bar e(t) \Omega_{21} \hat e(t-\tau(t)) \le \frac{\|\Omega_{21}\|^2}{\rho_3}\| \hat e\|^2_{\tau}+ \rho_3 \|\bar e(t)\|^2 \quad (ii)
```
Consider:
```math
2\bar e^\top (t)\Omega_{22}\bar e (t-\tau(t))
```
We can write it equivalently by adding and subtracting an instantanous term:
```math
\begin{aligned}
2\bar e^\top (t)\Omega_{22}\bar e (t-\tau(t)) = 2\bar e^\top (t)\Omega_{22}\bar e (t) + 2\bar e^\top (t)\Omega_{22}[\bar e (t-\tau(t)) - \bar e(t)]
\end{aligned} \quad (iii)
```
Since the matrices are real and symmetric we can use symmetrization and Rayleigh Ritz to write:$`^{[Ap.5]}`$
```math
[2\bar e^\top(t) \Omega_{22}\bar e(t) \le \lambda_{\max}(\Omega_{22} + \Omega_{22}^\top)\|\bar e \|^2] \quad (iv)
```
It can also be shown that:$`^{[Ap.6]}`$
```math
[2 \bar e^\top(t) \Omega_{22}[\bar e(t) - \bar e(t - \tau(t))]\le 2\tau_0 \|\bar e(t)\| \|\Omega_{22}\|   ((L + \|\Omega_{22}\|) \| \bar e\|_\tau + \|\Omega_{21}\|\|\hat e\|_{2\tau})] \tag v
```
Subbing $`(iv)`$ and $`(v)`$ into $`(iii)`$, letting $`\Omega_{22}^s := \frac 12 (\Omega_{22}+ \Omega_{22}^\top)`$, applying youngs with the constant $`\rho_4`$ and since $`\| \bar e(t) \| \le  \|\bar e\|_{\tau} \le  \|\bar e\|_{2\tau}`$
```math
\begin{aligned}
&2\bar e^\top (t)\Omega_{22}\bar e (t) + 2\bar e^\top (t)\Omega_{22}[\bar e (t-\tau(t)) - \bar e(t)] \\

&\le \lambda_{\max}(\Omega_{22} + \Omega_{22}^\top)\|\bar e \|^2 + 2\tau_0 \|\bar e(t)\| \|\Omega_{22}\|   ((L + \|\Omega_{22}\|) \| \bar e\|_\tau + \|\Omega_{21}\|\|\hat e\|_{2\tau}) \\

&\le 2\lambda_{\max}\Omega_{22}^s\|\bar e \|_{2\tau}^2 + 2\tau_0 \|\bar e\|_{2\tau} \|\Omega_{22}\|   ((L + \|\Omega_{22}\|) \| \bar e\|_{2\tau} + \|\Omega_{21}\|\|\hat e\|_{2\tau}) \\

&= 2\lambda_{\max}\Omega_{22}^s\|\bar e \|_{2\tau}^2 +   2\tau_0 \|\bar e\|^2_{2\tau} \|\Omega_{22}\| (L + \|\Omega_{22}\|)  +2\tau_0 \|\bar e\|_{2\tau} \|\Omega_{22}\|  \|\Omega_{21}\|\|\hat e\|_{2\tau} \\

&\le 2\lambda_{\max}\Omega_{22}^s\|\bar e \|_{2\tau}^2 +   2\tau_0 \|\bar e\|^2_{2\tau} \|\Omega_{22}\| (L + \|\Omega_{22}\|) +  \rho_4\|\bar e\|_{2\tau}^2 + \frac{\tau_0^2  \|\Omega_{22}\|^2  \|\Omega_{21}\|^2}{\rho_4}  \|\hat e\|_{2\tau}^2 \\

\end{aligned}
```
So:
```math
2\bar e^\top (t)\Omega_{22}\bar e (t-\tau(t)) \le  2\lambda_{\max}\Omega_{22}^s\|\bar e \|_{2\tau}^2 +   2\tau_0 \|\bar e\|^2_{2\tau} \|\Omega_{22}\| (L + \|\Omega_{22}\|) +  \rho_4\|\bar e\|_{2\tau}^2 + \frac{\tau_0^2  \|\Omega_{22}\|^2  \|\Omega_{21}\|^2}{\rho_4}  \|\hat e\|_{2\tau}^2 \quad (vi)
```
Subbing $`(ii)`$ and $`(vi)`$ into $`(i)`$ gives:
```math
\begin{aligned}
&\dot V_2(t)  \le 2 \gamma \| \bar e \|^2 + \frac{\|\Omega_{21}\|^2}{\rho_3}\| \hat e\|^2_{\tau}+ \rho_3 \|\bar e(t)\|^2  + 2\lambda_{\max}\Omega_{22}^s\|\bar e \|_{2\tau}^2 +   2\tau_0 \|\bar e\|^2_{2\tau} \|\Omega_{22}\| (L + \|\Omega_{22}\|) \\ &+  \rho_4\|\bar e\|_{2\tau}^2 + \frac{\tau_0^2  \|\Omega_{22}\|^2  \|\Omega_{21}\|^2}{\rho_4}  \|\hat e\|_{2\tau}^2 

\end{aligned}
```
From this we get:
```math
\dot V_2(t) \le -a_2V_2(t) + b_2\|V_2\|_{2\tau} + c_2\|V_1\|_{2\tau}, \quad t \ge \tau_0
```
Where:
```math
a_2:=-2[\lambda_{\max} (\Omega^s_{22})+ \gamma], \quad b_2:= \rho_3 + \rho_4 +2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|), \quad c_2:=\frac{\|\Omega_{21}\|^2}{\rho_3} +\frac{\tau_0^2 ||\Omega_{22}\|^2 \|\Omega_{21}\|^2}{\rho_4}
```
#### Proof - Part 3
---
*Obtaining the inequalities coefficient conditions*
Consider the two sets of coefficients:
```math
\begin{aligned}
&a_1 := 2(\bar d - \gamma - \sigma), \quad b_1 := \rho_1 + \frac{\|\Omega_{11}\|^2}{\rho_1} + \frac{\|\Omega_{12}\|^2}{\rho_2}, \quad c_1 :=\rho_2 \\ 
&a_2:=-2[\lambda_{\max} (\Omega^s_{22})+ \gamma], \quad b_2:= \rho_3 + \rho_4 +2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|), \quad c_2:=\frac{\|\Omega_{21}\|^2}{\rho_3} +\frac{\tau_0^2 ||\Omega_{22}\|^2 \|\Omega_{21}\|^2}{\rho_4}
\end{aligned}
```
We need $`h_i`$ defined in the following way to be positive:
```math
h_1 = a_1 - b_1 - c_1 > 0, \quad h_2 = a_2 - b_2 - c_1 >0
```
So therefore we need the following two equations to be posive:
```math
\begin{aligned}
&h_1 = 2(\bar d - \gamma - \sigma) - \left( \rho_1 + \frac{\|\Omega_{11}\|^2}{\rho_1} + \frac{\|\Omega_{12}\|^2}{\rho_2} \right) - \rho_2 \\ 
&h_2 = -2[\lambda_{\max} (\Omega^s_{22})+ \gamma] - (\rho_3 + \rho_4 +2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|)) - \left(\frac{\|\Omega_{21}\|^2}{\rho_3} + \frac{\tau_0^2 ||\Omega_{22}\|^2 \|\Omega_{21}\|^2}{\rho_4}\right)
\end{aligned}
```
Notice that $`\rho_i`$ are completely uncoupled so optimizing our free $`\rho_i`$ terms can be done completley with AM-GM. First $`h_1`$:
```math
\begin{aligned}
&h_1 = 2(\bar d - \gamma- \sigma) - \left( \rho_1 + \frac{\|\Omega_{11}\|^2}{\rho_1} + \frac{\|\Omega_{12}\|^2}{\rho_2} \right) - \rho_2 \\
&h_1 = 2(\bar d - \gamma- \sigma) - \left( \rho_1 + \frac{\|\Omega_{11}\|^2}{\rho_1} + \frac{\|\Omega_{12}\|^2}{\rho_2} + \rho_2\right)  
\end{aligned}
```
We want to minimize the bracketed term, so using AM-GM, optimal choices will yeild $`2\| \Omega_{11}\|`$ and $`2\| \Omega_{21} \|`$ in place of the $`\rho_i + \frac{\|\Omega_{ij}\|^2}{\rho_i}`$ type terms respectivley:
```math
h_1 = 2(\bar d - \gamma + \sigma - \| \Omega_{11}\| -\| \Omega_{12} \| )
```
So $`d_{\min} > \gamma + \sigma \| \Omega_{11}\| +\| \Omega_{21} \|`$, let $`\omega_{ij} := \|\Omega_{ij}\|`$ then this can be rewritten as:
```math
d_{\min} > \gamma + \omega_{11} + \omega_{12}
```
This is assumed, therefore we know that $`h_1 > 0`$.
Repeating the same for $`h_2`$:
```math
\begin{aligned}
&h_2 = -2[\lambda_{\max} (\Omega^s_{22})+ \gamma] - \rho_3 - \rho_4 - 2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|) - \frac{\|\Omega_{21}\|^2}{\rho_3} - \frac{\tau_0^2 ||\Omega_{22}\|^2 \|\Omega_{21}\|^2}{\rho_4} \\
&h_2 = -2[\lambda_{\max} (\Omega^s_{22})+ \gamma] - 2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|) - \left(\rho_3 + \frac{\|\Omega_{21}\|^2}{\rho_3} + \rho_4  + \frac{\tau_0^2 ||\Omega_{22}\|^2 \|\Omega_{21}\|^2}{\rho_4}\right) 
\end{aligned}
```
Applying AM-GM once again:
```math
\begin{aligned}
&h_2 = -2[\lambda_{\max} (\Omega^s_{22})+ \gamma] - 2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|) - 2\left(\|\Omega_{21}\| +  \tau_0 ||\Omega_{22}\| \|\Omega_{21}\|\right) \\

&h_2 = -2 \lambda_{\max}(\Omega^s_{22}) -2\gamma - 2\tau_0 \|\Omega_{22}\|(L + \|\Omega_{22}\|) -  2\|\Omega_{21}\| -  2\tau_0 ||\Omega_{22}\| \|\Omega_{21}\| \\

&h_2 = [-2 \lambda_{\max}(\Omega^s_{22}) -2\gamma  -  2\|\Omega_{21}\|] - \tau_0[ 2 \|\Omega_{22}\|(L + \|\Omega_{22}\|)+  2 ||\Omega_{22}\| \|\Omega_{21}\|] \\
\end{aligned}
```
Let $`N := [-2 \lambda_{\max}(\Omega^s_{22}) -2\gamma - 2\|\Omega_{21}\|]`$ and $`M:= [2 \|\Omega_{22}\|(L + \|\Omega_{22}\|)+  2 ||\Omega_{22}\| \|\Omega_{21}\|]`$ so: 
```math
\begin{aligned}
&h_2 = N - \tau_0 M > 0 
\end{aligned}
```
$`M`$ is positive and so is $`\tau_0`$, we need $`\tau_0 \lt N/M`$, which also demands that $`N>0`$, so we need:
```math
\tau_0 \lt \frac{-( \lambda_{\max}(\Omega^s_{22}) + \gamma + \|\Omega_{21}\|)}{2 \|\Omega_{22}\|(L + \|\Omega_{22}\|)+  ||\Omega_{22}\| \|\Omega_{21}\|}, \quad -\lambda_{\max}(\Omega^s_{22}) >  \gamma + \|\Omega_{21}\|
```
This can be rewritten using $`\omega_{ij} := \|\Omega_{ij}\|`$, and let $`\mu:= \lambda_{\max}(\Omega^s_{22})`$:
```math
\tau_0 \lt \frac{-\mu - \gamma - \omega_{21} }{ \omega_{22}(L + \omega_{22})+  \omega_{22} \omega_{21}}, \quad -\mu>  \gamma + \omega_{21}
```
Since this is also assumed, we have that $`h_1>0`$ and $`h_2 > 0`$, and we can conclude exponetial attractivity to the zero error state by the Halanay Lemma, which is a synchronized network.
### Appendix [1]
---
Systems Error Equation 
Let the first $`l`$ nodes be the ones that are controlled and the remaining $`N - l`$ nodes to be uncontrolled.
```math
\begin{aligned}
&\dot x(t) = f(t,x_i(t)) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau) + u_i(t), \quad 1, ... l \\
&\dot x(t) = f(t,x_i(t)) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau), \quad l+1... n
\end{aligned}\quad (1)
```
Define the synchronization error as $`e_i = x_i - z`$, with $`e = [e_1^\top, e_2^\top...e_N^\top]^\top`$, $`\hat e = [e_1^\top, e_2^\top...e_l^\top ]^\top`$ and $`\bar e = [e_{l+1}^\top... e_N^\top]^{\top}`$. 
Differentiating $`e`$ gives:
```math
\dot e = \dot x_i - \dot z
```
Subbing in $`\dot x_i`$ gives:
```math
\begin{aligned}
&\dot e(t) = f(t,x_i(t)) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau) + u_i(t) + \dot z, \quad 1, ... l \\
&\dot e(t) = f(t,x_i(t)) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau) + \dot z, \quad l+1... n
\end{aligned}
```
Subbing in $`\dot z = f(t,z(t))`$ and group it with the first term: $`\bar f_i(t) := f_i(t,x(t)) - f(t,z(t))`$:
```math
\begin{aligned}
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau) + u_i(t) , \quad 1, ... l \\
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma x_j(t-\tau) , \quad l+1... n
\end{aligned}
```
Now splitting the sum:
```math
\begin{aligned}
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau)+ c \sum_{j=1}^N a_{ij} \Gamma z(t-\tau) + u_i(t) , \quad 1, ... l \\
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau) + c \sum_{j=1}^N a_{ij} \Gamma z(t-\tau), \quad l+1... n
\end{aligned}
```
$`z(t-\tau)`$ is a constant scaling accross the sum so it can be factored out of the sum, at which point the sum is zero since $`\sum^N_{j=1} a_{ij} = 0`$ by assumption. So we are left with:
```math
\boxed{
\begin{aligned}
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau) + u_i(t) , \quad 1, ... l \\
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau), \quad l+1... n
\end{aligned}}
```
### Appendix [2]
---
Reorganizing the Error System
Recall:
```math
\begin{aligned}
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau) + u_i(t) , \quad 1, ... l \\
&\dot e(t) = \bar f_i(t) + c \sum_{j=1}^N a_{ij} \Gamma e_j(t-\tau), \quad l+1... n
\end{aligned} \quad (2)
```
TO -->
```math
\boxed{
\begin{aligned}
&\dot{\hat e}(t) = \hat F(t) + \Omega_{11}\hat e(t-\tau(t)) + \Omega_{12}\bar e(t-\tau(t)) - D\hat e(t) \\
&\dot{\bar e}(t) = \hat F(t) + \Omega_{21}\hat e(t-\tau(t)) + \Omega_{22}\bar e(t-\tau(t)) 
\end{aligned}}
```
### Appendix [3]
---
The following assumption:
```math
(x-y)^\top [f(t,x) - f(t,y)] \le \gamma \|x-y\|^2
```
Let $`x = e`$ and $`y = Z`$
Implying: 
```math
\boxed{\hat e^\top(t) \hat F(t) \le \gamma \|\hat e(t)\|^2}
```
### Appendix [4]
---
### Appendix [5]
---
Consider the following, where $`\Omega_{22}`$ is real and symmetric:
```math
2\bar e^\top(t) \Omega_{22}\bar e(t)
```
Since then matrix is real and symmetric we can use symmetrization and Reyleigh-Ritz: 
```math
2v^\top Mv =v^\top(M+M^\top)v \le \lambda_{\max}(M + M^\top)\|v\|^2
```
So:
```math
\boxed{
\begin{aligned}
&2\bar e^\top(t) \Omega_{22}\bar e(t) \le \lambda_{\max}(\Omega_{22} + \Omega_{22}^\top)\|\bar e \|^2
\end{aligned}}
```
### Appendix [6]
---
Consider the following:
```math
2 \bar e^\top(t) \Omega_{22}[\bar e(t) - \bar e(t - \tau(t))]
```
First using the fundemental theorem of calculus with the supremum over the integration bounds:
```math
\begin{aligned}
&\| \bar e(t) - \bar e(t-\tau(t)) \| = \| \int_{t-\tau(t)}^t \dot {\bar e(s)} \, ds \| \\
&\| \bar e(t) - \bar e(t-\tau(t)) \|  \le  \int_{t-\tau(t)}^t  \| \dot {\bar e(s)} \|  \, ds \\
& \le \tau_0 \sup_{s \in [t- \tau_0, t]} \| \bar {\dot e(s)} \| = \tau_0\|\bar{\dot e}\| _ \tau
\end{aligned}
```
We go onto apply this in the following way:
```math
2 \bar e^\top(t) \Omega_{22}[\bar e(t) - \bar e(t - \tau(t))] \le 2 \|\bar e(t)\| \|\Omega_{22}\| \|\bar e(t) - \bar e(t - \tau(t))\|
```
Plugging in our inequality:
```math
2 \bar e^\top(t) \Omega_{22}[\bar e(t) - \bar e(t - \tau(t))]\le 2\tau_0 \|\bar e(t)\| \|\Omega_{22}\|  \|\bar{\dot e}\| _ \tau \tag 1
```
We derived from the lipchitz condition that: SHOW WORK <-
```math
\| \dot {\bar e} \| \le (L + \|\Omega_{22}\|) \| \bar e(t)\| + \|\Omega_{21}\|\|\hat e\|_\tau
```
Taking the supremum over the lookback window on both sides:
```math
\| \dot {\bar e} \|_ \tau \le (L + \|\Omega_{22}\|) \| \bar e\|_\tau + \|\Omega_{21}\|\|\hat e\|_2\tau \quad (2)
```
Plugging this back into $`(1)`$ we have that:
```math
\boxed{2 \bar e^\top(t) \Omega_{22}[\bar e(t) - \bar e(t - \tau(t))]\le 2\tau_0 \|\bar e(t)\| \|\Omega_{22}\|(   (L + \|\Omega_{22}\|) \| \bar e\|_\tau + \|\Omega_{21}\|\|\hat e\|_2\tau)}
```
### Appendix [N]
---
Remarks
These are the revised conditions:
```math
\boxed{
\begin{aligned}
&2\bar e^\top(t) \Omega_{22}\bar e(t) \le \lambda_{\max}(\Omega_{22} + \Omega_{22}^\top)\|\bar e \|^2
\end{aligned}}
```
```math
\boxed{2 \bar e^\top(t) \Omega_{22}[\bar e(t) - \bar e(t - \tau(t))]\le 2\tau_0 \|\bar e(t)\| \|\Omega_{22}\|(   (L + \|\Omega_{22}\|) \| \bar e\|_\tau + \|\Omega_{21}\|\|\hat e\|_2\tau)}
```
?Difference 1: Grouped under a single matrix control with a single target state?
-> ?Taking the minimum on all control as a bound. Definetly weakens the result.?
An equivalent Condition set that is more convinient for those applying the theorem:
```math
\text{(C1)}\quad d_{\min} > \gamma + c\big(\|A_{11}\otimes\Gamma\| + \|A_{21}\otimes\Gamma\|\big)
```
```math
\text{(C2)}\quad \tfrac{c}{2}\lambda_{\max}\big((A_{22}+A_{22}^\top)\otimes\Gamma\big) < -\big(\gamma + \|\Omega_{21}\|\big)
```
```math
\text{(C3)}\quad \tau_0 < \tau^* := \frac{-\lambda_{\max}(\Omega_{22}^s) - \gamma - \|\Omega_{21}\|}{\|\Omega_{22}\|\big(L + \|\Omega_{22}\| + \|\Omega_{21}\|\big)} \quad (\tau^* := +\infty \text{ if } \Omega_{22} = 0)
```
---
```math
\begin{aligned}
&-c \lambda_{\max}([A_{22}\otimes \Gamma]) > \gamma \quad (Th-C.1,2,3) \\
&\tau_0 \lt \tau_{\max} := \zeta_1 / \zeta_2 \\
&\bar d > d_{\min} := \gamma + c\|\Gamma\|\|A_{11}\| + \frac{c^2\zeta_3^2 \|\Gamma\|^2}{4(\zeta_1 - \tau_0 \zeta_2) }
\end{aligned}
```

Where:
- $`\bar d = \min\{d_1,...d_l\}`$
- $`\zeta_1 = -c\lambda_{\max}([A_{22}\otimes \Gamma])`$
- $`\zeta_2 = cL\|\Gamma\| \|A_{22}\| + c^2 \|\Gamma\|^2 \|A_{22}\|^2`$
- $`\zeta_ 3 = \|A_{12}\| + \|A_{21}\| + c \tau_0 \|A_{21}\| \|A_{22}\|\|\Gamma\|`$

---
