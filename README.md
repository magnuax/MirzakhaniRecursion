# MirzakhaniRecursion
 
 
Implements a recursive algorithm for calculating Weil-Peterson volumes, using Mirzakhani's recursion relations (see [1]). Probably not the most efficient solution, but it was a fun project.


## The Recursion Relation
For genus $g$ and $n$ boundary components with lengths $L_1, ..., L_n$, we denote the Weil-Peterson volume by $V_{g,n}(L_1,...,L_n)$. Using the base cases

* $V_{0,1} = 0$

* $V_{1,2} = 0$

* $V_{0,3} = 1$

* $V_{1,1}(L_1) = \frac{L_1^2}{48} + \frac{\pi^2}{12}$

The recursion relation for Weil–Petersson volumes $V_{g,n}$ are as follows, for $2g + n > 3$:
```math
2 \frac{\partial}{\partial L_1} L_1 V_{g,n}(L) = \int_0^\infty \int_0^\infty xy H(x + y, L_1) V_{g-1,n+1}(x, y, \hat{L}) \, dx \, dy
\\
+ \sum_{\substack{g_1 + g_2 = g \\ I \sqcup J = \{2, \ldots, n\}}} \int_0^\infty \int_0^\infty xy H(x + y, L_1) V_{g_1,|I|+1}(x, L_I) V_{g_2,|J|+1}(y, L_J) \, dx \, dy
\\
+ \sum_{k=2}^n \int_0^\infty x \left( H(x, L_1 + L_k) + H(x, L_1 - L_k) \right) V_{g,n-1}(x, \hat{L}_k) \, dx
```
Where $\hat{L} = (L_2, L_3, \ldots, L_n)$, $L_I = (L_{i_1}, L_{i_2}, \ldots, L_{i_m})$ for $I = \{i_1, i_2, \ldots, i_m\}$, and $\hat{L}_k = (L_2, \ldots, \hat{L}_k, \ldots, L_n)$, with the hat denoting omission. The function $H$ is defined as
```math
H(x, y) = \frac{1}{1 + \exp\left(\frac{x+y}{2}\right)} + \frac{1}{1 + \exp\left(\frac{x-y}{2}\right)}.
```
## ToDo:
- Not working for $g>2$ due to problems with $V_{2,0}$.
- Fix this readme

## References
 [1] Mirzakhani, M. Simple geodesics and Weil-Petersson volumes of moduli spaces of bordered Riemann surfaces. Invent. math. 167, 179–222 (2007). https://doi.org/10.1007/s00222-006-0013-2
