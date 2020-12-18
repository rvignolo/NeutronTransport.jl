```@meta
CurrentModule = NeutronTransport
```

# NeutronTransport

```@index
```

```@autodocs
Modules = [NeutronTransport]
```

# Motivation

We have selected the Method of Characteristics (MOC) as the first implemented lattice formulation since it overcomes the main limitation for the Collision Probability (CP) method: it does not produce full square matrices of order equal to the number of regions in the domain times the energy groups, so it is possible to solve problems with many more regions.

In this context, an efficient ray tracing algorithm for unstructured meshes (see [RayTracing.jl](https://github.com/rvignolo/RayTracing.jl)) and a non-linear power method were developed. These two things allow computing the scalar flux per region and energy group by collecting all mean angular fluxes over track segments and region sources.

# Introduction

Lattice computations, where a representative lattice is solved by means of an adequate transport formulation, are the first step of deterministic reactor analysis and calculations. The Collision Probability Method (CP) and the Method of Characteristics (MOC) are among these formulations. Both of them compute the scalar flux over well-discretized energy and region meshes that is later used to obtain homogenized cross sections for core-level calculations.

At first, representative lattices used to be unique 2D square *pin cells* geometries. Nowadays, with the increment of computing power, lattice domains became arrangements of *pin cells* that may represent one or many fuel elements. Collision Probability methods can handle these complex geometries but they become impractical as the numer of regions increases. On the other hand, over the last years, the Method of Characteristics gain particular interest because it can handle fuel arrangements or even complete cores, overcoming CP limitations.

## Method of Characteristics

The Method of Characteristics obtains the neutron flux distribution by solving the characteristic form of the transport equation over tracks that emulate neutron trajectories across the problem domain.

## Nomenclature


## Formulation

The steady-state integro-differential neutron transport equation is extensively discussed in the bibliography. In its multigroup form, the continuous dependence on energy is replaced by a discrete scheme, yielding:

```math
\boldsymbol{\hat{\Omega}} \cdot \nabla \psi_g(\vec{x}, \boldsymbol{\hat{\Omega}})
 + \Sigma^t_g(\vec{x}) \cdot \psi_g(\vec{x}, \boldsymbol{\hat{\Omega}}) = q_g(\vec{x}, \boldsymbol{\hat{\Omega}}),
```
