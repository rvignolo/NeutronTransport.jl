```@meta
CurrentModule = NeutronTransport
```

# NeutronTransport

```@index
```

```@autodocs
Modules = [NeutronTransport]
```

# Introduction

Lattice computations, where a representative lattice is solved by means of an adequate transport formulation, are the first step of deterministic reactor analysis and calculations. The Collision Probability Method (CP) and the Method of Characteristics (MOC) are among these formulations. Both of them compute the scalar flux over well-discretized energy and region meshes that is later used to obtain homogenized cross sections for core-level calculations.

At first, representative lattices used to be unique 2D square *pin cells* geometries. Nowadays, with the increment of computing power, lattice domains became arrangements of *pin cells* that may represent one or many fuel elements. Collision Probability methods can handle these complex geometries but they become impractical as the numer of regions increases. On the other hand, over the last years, the Method of Characteristics gain particular interest because it can handle fuel arrangements or even complete cores, overcoming CP limitations.

## Method of Characteristics

The Method of Characteristics (MOC) was the selected formulation since it overcomes the main limitation for the Collision Probability (CP) method: it does not produce full square matrices of order equal to the number of regions in the domain times the energy groups, so it is possible to solve problems with many more regions.

In this context, an efficient ray tracing algorithm for unstructured meshes (see [RayTracing.jl](https://github.com/rvignolo/RayTracing.jl)) was developed in first place. Now, this package allows computing the scalar flux per region and energy group by collecting all mean angular fluxes over track segments and region sources.

The Method of Characteristics obtains the neutron flux distribution by solving the characteristic form of the transport equation over tracks that emulate neutron trajectories across the problem domain.

## Nomenclature

| .                              |       .                     |  .              |          .             |
|--------------------------------|-----------------------------|-----------------|------------------------|
| **Symbols:**                   |                             |                 |                        |
| ``\vec{x}``                    | Position vector.            | ``\ell``        | Track segment length.  |
| ``\boldsymbol{\hat{\Omega}}``  | Direction of motion versor. | ``k``           | Multiplication factor. |
| ``\varphi``                    | Azimuthal angle.            | ``Q``           | Scalar source.         |
| ``\theta``                     | Polar angle.                | ``q``           | Angular source.        |
| ``\psi``                       | Angular flux.                |                 |                        |
| ``\phi``                       | Scalar flux.                 | **Subindexes:** |                        |
| ``\Sigma^t``                   | Total cross-section.        | ``g``           | Energy group.          |
| ``\Sigma^a``                   | Absorption cross-section.   | ``i``           | Flat source region.    |
| ``\Sigma^f``                   | Fission cross-section.      | ``a``           | Azimuthal angle.       |
| ``\Sigma^s``                   | Scattering cross-section.   | ``p``           | Polar angle.           |
| ``\nu``                        | Neutrons per fission.        | ``m``           | Solid angle.           |
| ``\chi``                       | Fission spectrum.           | ``k``           | Track segment.         |


## Formulation

The steady-state integro-differential neutron transport equation is extensively discussed in the bibliography. In its multigroup form, the continuous dependence on energy is replaced by a discrete scheme, yielding:

```math
\boldsymbol{\hat{\Omega}} \cdot \nabla \psi_g(\vec{x}, \boldsymbol{\hat{\Omega}})
 + \Sigma^t_g(\vec{x}) \cdot \psi_g(\vec{x}, \boldsymbol{\hat{\Omega}}) = q_g(\vec{x}, \boldsymbol{\hat{\Omega}}).
```

This equation defines a linear first-order partial differential equation (PDE) if we assume that the flux is known inside the source term. Then, it can be casted to a family of linear first-order ordinary differential equations (ODEs) by means of the method of characteristics. In this context, we would want to transform the linear first-order PDE into an ODE along the appropriate curve $\vec{x}(s)$; i.e. something of the form:
