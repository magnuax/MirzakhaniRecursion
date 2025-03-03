# MirzakhaniRecursion
 
 
Implements a recursive algorithm for calculating Weil-Peterson volumes in SymPy, using Mirzakhani's recursion relations (see [1]). 
Mostly following the procedure shown in section 4. of [2]. For the cases with no boundaries, the Dilaton equation is used instead (Theorem 22 in [2]).
Probably not the most efficient implementation, but it works and it was a fun project.

## Mirzakhani's Recursion Relation
For genus $g$ and $n$ boundary components with lengths $L_1, ..., L_n$, we denote the Weil-Peterson volume by $V_{g,n}(L_1,...,L_n)$. Using the base cases

* $V_{0,1} = 0$

* $V_{1,2} = 0$

* $V_{0,3} = 1$

* $V_{1,1}(L_1) = \frac{L_1^2}{48} + \frac{\pi^2}{12}$

The recursion relation for Weil–Petersson volumes $V_{g,n}$ are as follows, for $2g + n > 3$:
```math
\begin{aligned}
2 \frac{\partial}{\partial L_1} L_1 V_{g,n}(L) = &\int_0^\infty \int_0^\infty xy H(x + y, L_1) V_{g-1,n+1}(x, y, \hat{L}) \, dx \, dy
\\
+ \sum_{\substack{g_1 + g_2 = g \\ I \sqcup J = \{2, \ldots, n\}}} &\int_0^\infty \int_0^\infty xy H(x + y, L_1) V_{g_1,|I|+1}(x, L_I) V_{g_2,|J|+1}(y, L_J) \, dx \, dy
\\
+ \sum_{k=2}^n &\int_0^\infty x \left( H(x, L_1 + L_k) + H(x, L_1 - L_k) \right) V_{g,n-1}(x, \hat{L}_k) \, dx
\end{aligned}
```
Where $\hat{L} = (L_2, L_3, \ldots, L_n)$, $L_I = (L_{i_1}, L_{i_2}, \ldots, L_{i_m})$ for $I = \{i_1, i_2, \ldots, i_m\}$, and $\hat{L}_k = (L_2, \ldots, \hat{L}_k, \ldots, L_n)$, with the hat denoting omission. The function $H$ is defined as
```math
H(x, y) = \frac{1}{1 + \exp\left(\frac{x+y}{2}\right)} + \frac{1}{1 + \exp\left(\frac{x-y}{2}\right)}.
```

## The Dilaton equation
Since the proof of Mirzakhani's relations requires at least one boundary component, we need some other way to compute cases where $n=0$.
So, to compute the volumes $V_{g,0}$, we instead use the so-called dilaton equation [2]:
```math
\frac{\partial V_{g,n+1}}{\partial L_{n+1}}(L,2\pi i ) = 2\pi i (2g-2+n) V_{g,n}(L)	
```

## ToDo:
- Fix this readme
- Prettier code
- Make more efficient
- (tried using sp.Poly instead of general expressions - turned out a bit messy?)
- Better documentation
## References
 [1] Mirzakhani, M. Simple geodesics and Weil-Petersson volumes of moduli spaces of bordered Riemann surfaces. Invent. math. 167, 179–222 (2007). https://doi.org/10.1007/s00222-006-0013-2
 \
 [2] N. Do, Moduli spaces of hyperbolic surfaces and their WeilPetersson volumes. [arXiv:1103.4674.](https://doi.org/10.48550/arXiv.1103.4674) 

