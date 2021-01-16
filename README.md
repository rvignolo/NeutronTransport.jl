# NeutronTransport

| **Documentation** |
|:------------ |
| [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rvignolo.github.io/NeutronTransport.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rvignolo.github.io/NeutronTransport.jl/dev/) |
|**Build Status** |
| [![Build Status](https://github.com/rvignolo/NeutronTransport.jl/workflows/CI/badge.svg)](https://github.com/rvignolo/NeutronTransport.jl/actions) |

## Description

**NeutronTransport** is a reactor physics program that solves the steady-state multigroup [neutron transport equation](https://en.wikipedia.org/wiki/Neutron_transport#Neutron_transport_equation) by means of the Method of Characteristics approximation over unstructured grids. It relies on [RayTracing.jl](https://github.com/rvignolo/RayTracing.jl) for the tracking procedure.

## Examples

These are some popular examples solved with NeutronTransport:

| ![](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/pincell-g1.png)  | ![](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/bwr-g2.png) | ![](demo/c5g7-g7.png) |
|:-------------:|:-------------:|:-------------:|
| [*Pincell*](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/pincell.jl) | [4 by 4 BWR lattice (2 Gd pins)](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/bwr.jl) | [C5G7 Benchmark](https://github.com/rvignolo/NeutronTransport.jl/blob/main/demo/c5g7.jl) |
