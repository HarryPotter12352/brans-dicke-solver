# Complete Brans-Dicke theory solver

In this Python script, I provide code to compute the entirety of the Einstein tensor in the Brans-Dicke theory of gravitation where a scalar field is involved.
The parameters required for computation are

- $\begin{bmatrix}x^0 \\ x^1 \\ x^2 \\ x^3\end{bmatrix}$, the spacetime coordinates
- $g_{\mu\nu}$, the metric
- $T_{\mu\nu}$, the energy momentum
- $\phi$, the scalar field
- $V(\phi)$, the potential associated with the scalar field
- $\omega$, the dimensionless Brans-Dicke coupling parameter.

The equations we use here are given below.

The Einstein tensor is given by
$
G_{\mu\nu} = \frac{8\pi}{\phi}T_{\mu\nu} + \frac{\omega}{\phi^2}\left(\partial_{\mu}\phi\partial_{\nu}\phi - \frac{1}{2}g_{\mu\nu}\partial_{\sigma}\phi\partial^{\sigma}\phi\right) + \frac{1}{\phi}(\partial_{\mu}\partial_{\nu}\phi - g_{\mu\nu}\Box\phi) - g_{\mu\nu} \frac{V(\phi)}{2\phi} 
$
where $\Box\phi$ is the Laplace-Beltrami operator given by
$
\Box\phi = \frac{8\pi T + 2V(\phi) - \phi V'(\phi)}{3+2\omega}
$

~~Website coming soon~~