# Benchmarks

Run from the repository root:

```julia
julia --project=benchmark -e 'import Pkg; Pkg.develop([Pkg.PackageSpec(path="../RayTracing.jl"), Pkg.PackageSpec(path=".")]); Pkg.instantiate(); include("benchmark/benchmarks.jl")'
```

For a quick smoke run:

```julia
julia --project=benchmark benchmark/benchmarks.jl --quick
```

The suite uses small synthetic fixtures plus selected meshes from `demo/` so benchmark
coverage can grow without turning the demo scripts into test harnesses.
