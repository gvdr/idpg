# Deprecated PDE Code (January 2026)

This directory contains archived PDE evolution code that was replaced during the January 2026 refactoring.

## Files

### PDEEvolution_explicit_euler.jl
Original hand-coded PDE evolution with explicit Euler time-stepping.

**Why archived:**
- Explicit Euler has poor stability (requires very small dt for stiff problems)
- No adaptive time-stepping or error control
- Replaced with OrdinaryDiffEq.jl integration

### PDESciML_methodoflines.jl
MethodOfLines.jl + ModelingToolkit symbolic PDE solver.

**Why archived:**
- MethodOfLines.jl has "Low maturity" status (per SciML classification)
- 96 open issues, known scalability problems
- Broke between Julia 1.10 and 1.11
- Heavy dependencies (ModelingToolkit, DomainSets) for limited benefit
- Our B^d_+ domain requires masking anyway

## Reverting

If you need to restore this code:

```julia
# Copy back to src/modules/
cp archive/deprecated_pde/PDEEvolution_explicit_euler.jl src/modules/PDEEvolution.jl
cp archive/deprecated_pde/PDESciML_methodoflines.jl src/modules/PDESciML.jl

# Re-add to Project.toml:
# MethodOfLines = "94925ecb-adb7-4558-8ed8-f975c56a0bf4"
# ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
# DomainSets = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"

# Re-add include and exports in IDPG.jl
```

See `docs/pde_refactoring_2025.md` for full rationale.
