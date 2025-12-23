# Known Issues and Solutions

## MethodOfLines.jl 4D Advection Boundary Conditions

**Status**: Resolved
**Package**: MethodOfLines.jl
**Date discovered**: December 2025

### Issue

When using MethodOfLines.jl for 4D advection equations with only Dirichlet (inflow) boundary conditions, the discretization throws a `BoundsError`. The error occurs during symbolic discretization when the upwind stencil tries to access indices outside the valid array bounds.

### Original Error Message

```
BoundsError: attempt to access 5x5x5x5 Array{SymbolicUtils.BasicSymbolic{Real}, 4} at index [CartesianIndex{4}[CartesianIndex(0, 2, 2, 2), CartesianIndex(1, 2, 2, 2), ...]]
```

### Root Cause

The upwind scheme requires information from neighboring cells for boundary handling. With only inflow (x=0) boundary conditions specified, the scheme cannot properly handle the outflow (x=1) boundaries and attempts to access ghost points outside the array.

### Solution

Add Neumann (zero gradient) boundary conditions at the outflow boundaries in addition to the Dirichlet conditions at inflow:

```julia
@parameters t x1 x2 x3 x4
@variables rho(..)
Dx1 = Differential(x1)
Dx2 = Differential(x2)
Dx3 = Differential(x3)
Dx4 = Differential(x4)

bcs = [
    rho(0.0, x1, x2, x3, x4) ~ 0.0,  # Initial condition
    # Dirichlet (zero value) at x=0 boundaries (inflow)
    rho(t, 0.0, x2, x3, x4) ~ 0.0,
    rho(t, x1, 0.0, x3, x4) ~ 0.0,
    rho(t, x1, x2, 0.0, x4) ~ 0.0,
    rho(t, x1, x2, x3, 0.0) ~ 0.0,
    # Neumann (zero gradient) at x=1 boundaries (outflow)
    Dx1(rho(t, 1.0, x2, x3, x4)) ~ 0.0,
    Dx2(rho(t, x1, 1.0, x3, x4)) ~ 0.0,
    Dx3(rho(t, x1, x2, 1.0, x4)) ~ 0.0,
    Dx4(rho(t, x1, x2, x3, 1.0)) ~ 0.0,
]
```

### Reference

From the [MethodOfLines.jl documentation](https://docs.sciml.ai/MethodOfLines/dev/):

> Note that the WENO Scheme is often unstable in more than 1 spatial dimension due to difficulties with boundary handling, this can be avoided if you supply 2 or more bcs per boundary in the dimension along which an advection term acts.

While this note specifically mentions WENO Scheme, the same principle applies to UpwindScheme in higher dimensions: providing complete boundary conditions (both Dirichlet and Neumann) at each boundary ensures the discretization has sufficient information to construct the stencils.

### Performance Note

4D PDE solving with MethodOfLines is computationally intensive. With resolution=5 (5^4=625 grid points), the 4D advection test takes approximately 16 minutes. For faster iteration during development, use smaller resolutions (resolution=4) or coarser time steps.
