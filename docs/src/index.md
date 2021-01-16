# Introduction

**NeutronTransport** is a high-performance library designed to achieve fast and advanced reactor physics calculations.

Currently, there isn't any reactor physics program written in Julia. **NeutronTransport** main goal is to open that gate so people from the field get to know this amazing programming language. On the other hand, since it is free software (as in free speech), it may be developed collaboratively by volunteer engineers/computer programmers who do not want to be, as [StammlerAbbate](@cite) used to say, users that *will never understand the programs they are to use and, as computer slaves, consider them as black boxes, blindly trusting their results*.

# Getting Started

The package can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the `Pkg` REPL mode and run:

```julia
pkg> add NeutronTransport
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("NeutronTransport")
```
