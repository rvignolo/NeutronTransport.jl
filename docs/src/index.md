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

The Method of Characteristics obtains the neutron flux distribution by solving the characteristic form of the transport equation over tracks that emulate neutron trajectories across the problem domain.

It was the selected formulation since it overcomes the main limitation for the Collision Probability (CP) method: it does not produce full square matrices of order equal to the number of regions in the domain times the energy groups, so it is possible to solve problems with many more regions.

In this context, an efficient ray tracing algorithm for unstructured meshes (see [RayTracing.jl](https://github.com/rvignolo/RayTracing.jl)) was developed in first place. Now, this package allows computing the scalar flux per region and energy group by collecting all mean angular fluxes over track segments and region sources.

## Nomenclature

The following is the nomenclature used in the documentation:

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

This equation defines a linear first-order partial differential equation (PDE) if we assume that the flux is known inside the source term ``q_g(\vec{x}, \boldsymbol{\hat{\Omega}})``. Then, it can be casted to a family of linear first-order ordinary differential equations (ODEs) by means of the method of characteristics. In this context, we would want to transform the linear first-order PDE into an ODE along the appropriate curve ``\vec{x}(s)``; i.e. something of the form:

```math
\frac{d}{ds} \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}}) = F(\vec{x}(s), \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}})).
```

Using the chain rule, we get:

```math
\frac{d}{ds} \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}}) = \frac{\partial \psi}{\partial x} \frac{dx}{ds} +
                     \frac{\partial \psi}{\partial y} \frac{dy}{ds} +
                     \frac{\partial \psi}{\partial z} \frac{dz}{ds} =
                     \frac{d \vec{x}}{ds} \cdot \nabla \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}}),
```

which, inserted in equation the previous equation, yields:

```math
\frac{d \vec{x}}{ds} \cdot \nabla \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}}) = F(\vec{x}(s), \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}})).
```

By setting ``d \vec{x} / ds = \boldsymbol{\hat{\Omega}}``, we get:

```math
\boldsymbol{\hat{\Omega}} \cdot \nabla \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}}) = F(\vec{x}(s), \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}})),
```

which is equivalent to the original PDE if:

```math
F(\vec{x}(s), \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}})) = q(\vec{x}(s), \boldsymbol{\hat{\Omega}}) - \Sigma^t(\vec{x}(s)) \cdot \psi(\vec{x}(s), \boldsymbol{\hat{\Omega}}).
```

Now, solving equation ``\vec{x}(s)`` yields to the parametric equations of the characteristic curves:

```math
 \vec{x}(s) = \vec{x}_0 + s \boldsymbol{\hat{\Omega}},
```

for arbitrary ``\vec{x}_0``. It can be conclude that, along the characteristic curves ``\vec{x}(s)``, the original PDE becomes the following ODE:

```math
\frac{d}{ds} \psi_g(s, \boldsymbol{\hat{\Omega}}) + \Sigma^t_g(s) \cdot \psi_g(s, \boldsymbol{\hat{\Omega}}) = q_g(s, \boldsymbol{\hat{\Omega}}),
```

where the abuse of notation ``\vec{x}_0 + s \boldsymbol{\hat{\Omega}} \rightarrow s`` has been used.

At the beginning of the section, we stated that the source term is known. However, in reality, the source term in the transport equation is unknown, since it is given by:

```math
q_g(s, \boldsymbol{\hat{\Omega}}) =
 \sum_{g^\prime=1}^G \int_{4\pi} \Sigma^s_{g^\prime \rightarrow g}(s, \boldsymbol{\hat{\Omega}^\prime} \rightarrow \boldsymbol{\hat{\Omega}}) \cdot \psi_{g^\prime}(s, \boldsymbol{\hat{\Omega}^\prime}) \, d\boldsymbol{\hat{\Omega}^\prime} + \frac{\chi_g(s)}{4\pi k_{\text{eff}}} \sum_{g^\prime=1}^G \int_{4\pi} \nu\Sigma^f_{g^\prime}(s) \cdot \psi_{g^\prime}(s, \boldsymbol{\hat{\Omega}^\prime}) \, d\boldsymbol{\hat{\Omega}^\prime}.
```

for eigenvalue problems (without external sources). Additionally, when considering isotropic scattering, the previous expression reduces to:

```math
 q_g(s, \boldsymbol{\hat{\Omega}}) =
 \frac{1}{4\pi} \sum_{g^\prime=1}^G \Sigma^s_{g^\prime \rightarrow g}(s) \cdot \phi_{g^\prime}(s) + \frac{\chi_g(s)}{4\pi k_{\text{eff}}} \sum_{g^\prime=1}^G \nu\Sigma^f_{g^\prime}(s) \cdot \phi_{g^\prime}(s)  = q_g(s).
```

Consequently, even though the ODE may be integrated to find a solution, an iterative scheme over the source term must be used to converge to the true solution. Commonly, each of these iterations is known as a *transport sweep*.