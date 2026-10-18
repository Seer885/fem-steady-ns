Spent Summer 2026 studying how to implement FEM to the steady Navier-Stokes PDE, which was done in stages to make debugging as easy as possible.

1) Applied FEM to the Stokes equations (system of PDE's)
2) Applied FEM to both the scalar, then normal convection-diffusion PDE's (non-linearity)
3) Applied FEM to the Oseen equations, aka linearized Navier-Stokes (system of PDE's)
4) Applied FEM to steady Navier-Stokes (system of PDE's with non-linearity)
5) Implemented Anderson Acceleration and SUPG stabilization to improve the convergence of the nonlinear Picard iterative solver
6) Applied the above techniques to generate streamline plot solutions to the 2D lid-driven cavity problem

Most of the above scripts used Q2/Q1 elements with uniform distributions, though the lid-driven cavity scripts used P2/P1 FEM with uniform distributions to try to emulate the results in the paper “Anderson-Accelerated Convergence of Picard Iterations for Incompressible Navier–Stokes Equations” by Sara Pollock et al. 
