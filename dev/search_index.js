var documenterSearchIndex = {"docs":
[{"location":"references.html#References","page":"References","title":"References","text":"","category":"section"},{"location":"references.html","page":"References","title":"References","text":"","category":"page"},{"location":"moc.html#Method-of-Characteristics","page":"Method of Characteristics","title":"Method of Characteristics","text":"","category":"section"},{"location":"moc.html#Introduction","page":"Method of Characteristics","title":"Introduction","text":"","category":"section"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"Lattice computations, where a representative lattice is solved by means of an adequate transport formulation, are the first step of deterministic reactor calculations and analysis. The Collision Probability Method (CP) and the Method of Characteristics (MOC) are among these formulations. Both of them compute the scalar flux over well-discretized energy and region meshes that is later used to obtain homogenized cross sections for core-level calculations.","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"At first, representative lattices used to be unique 2D square pin cells geometries. Nowadays, with the increment of computing power, lattice domains became arrangements of pin cells that may represent one or many fuel elements. Collision Probability methods can handle these complex geometries but they become impractical as the numer of regions increases. On the other hand, over the last years, the Method of Characteristics gain particular interest because it can handle fuel arrangements or even complete cores, overcoming CP limitations.","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"The Method of Characteristics obtains the neutron flux distribution by solving the characteristic form of the transport equation over tracks that emulate neutron trajectories across the problem domain.","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"It was the selected formulation since it overcomes the main limitation for the Collision Probability (CP) method: it does not produce full square matrices of order equal to the number of regions in the domain times the energy groups, so it is possible to solve problems with many more regions.","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"In this context, an efficient ray tracing algorithm for unstructured meshes (see RayTracing.jl) was developed in first place. Now, this package allows computing the scalar flux per region and energy group by collecting all mean angular fluxes over track segments and region sources.","category":"page"},{"location":"moc.html#Nomenclature","page":"Method of Characteristics","title":"Nomenclature","text":"","category":"section"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"The following is the nomenclature used in the documentation:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":". . . .\nSymbols:   \nvecx Position vector. ell Track segment length.\nboldsymbolhatOmega Direction of motion versor. k Multiplication factor.\nvarphi Azimuthal angle. Q Scalar source.\ntheta Polar angle. q Angular source.\npsi Angular flux.  \nphi Scalar flux. Subindexes: \nSigma^t Total cross-section. g Energy group.\nSigma^a Absorption cross-section. i Flat source region.\nSigma^f Fission cross-section. a Azimuthal angle.\nSigma^s Scattering cross-section. p Polar angle.\nnu Neutrons per fission. m Solid angle.\nchi Fission spectrum. k Track segment.","category":"page"},{"location":"moc.html#Formulation","page":"Method of Characteristics","title":"Formulation","text":"","category":"section"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"The steady-state integro-differential neutron transport equation is extensively discussed in the bibliography. In its multigroup form, the continuous dependence on energy is replaced by a discrete scheme, yielding:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"boldsymbolhatOmega cdot nabla psi_g(vecx boldsymbolhatOmega)\n + Sigma^t_g(vecx) cdot psi_g(vecx boldsymbolhatOmega) = q_g(vecx boldsymbolhatOmega)","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"This equation defines a linear first-order partial differential equation (PDE) if we assume that the flux is known inside the source term q_g(vecx boldsymbolhatOmega). Then, it can be casted to a family of linear first-order ordinary differential equations (ODEs) by means of the method of characteristics. In this context, we would want to transform the linear first-order PDE into an ODE along the appropriate curve vecx(s); i.e. something of the form:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"fracdds psi(vecx(s) boldsymbolhatOmega) = F(vecx(s) psi(vecx(s) boldsymbolhatOmega))","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"Using the chain rule, we get:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"fracdds psi(vecx(s) boldsymbolhatOmega) = fracpartial psipartial x fracdxds +\n                     fracpartial psipartial y fracdyds +\n                     fracpartial psipartial z fracdzds =\n                     fracd vecxds cdot nabla psi(vecx(s) boldsymbolhatOmega)","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"which, inserted in the previous equation, yields:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"fracd vecxds cdot nabla psi(vecx(s) boldsymbolhatOmega) = F(vecx(s) psi(vecx(s) boldsymbolhatOmega))","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"By setting d vecx  ds = boldsymbolhatOmega, we get:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"boldsymbolhatOmega cdot nabla psi(vecx(s) boldsymbolhatOmega) = F(vecx(s) psi(vecx(s) boldsymbolhatOmega))","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"which is equivalent to the original PDE if:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"F(vecx(s) psi(vecx(s) boldsymbolhatOmega)) = q(vecx(s) boldsymbolhatOmega) - Sigma^t(vecx(s)) cdot psi(vecx(s) boldsymbolhatOmega)","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"Now, solving equation vecx(s) yields to the parametric equations of the characteristic curves:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":" vecx(s) = vecx_0 + s boldsymbolhatOmega","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"for arbitrary vecx_0. It can be conclude that, along the characteristic curves vecx(s), the original PDE becomes the following ODE:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"fracdds psi_g(s boldsymbolhatOmega) + Sigma^t_g(s) cdot psi_g(s boldsymbolhatOmega) = q_g(s boldsymbolhatOmega)","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"where the abuse of notation vecx_0 + s boldsymbolhatOmega rightarrow s has been used.","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"At the beginning of the section, we stated that the source term is known. However, in reality, the source term in the transport equation is unknown, since it is given by:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"q_g(s boldsymbolhatOmega) =\n sum_g^prime=1^G int_4pi Sigma^s_g^prime rightarrow g(s boldsymbolhatOmega^prime rightarrow boldsymbolhatOmega) cdot psi_g^prime(s boldsymbolhatOmega^prime)  dboldsymbolhatOmega^prime + fracchi_g(s)4pi k_texteff sum_g^prime=1^G int_4pi nuSigma^f_g^prime(s) cdot psi_g^prime(s boldsymbolhatOmega^prime)  dboldsymbolhatOmega^prime","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"for eigenvalue problems (without external sources). Additionally, when considering isotropic scattering, the previous expression reduces to:","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":" q_g(s boldsymbolhatOmega) =\n frac14pi sum_g^prime=1^G Sigma^s_g^prime rightarrow g(s) cdot phi_g^prime(s) + fracchi_g(s)4pi k_texteff sum_g^prime=1^G nuSigma^f_g^prime(s) cdot phi_g^prime(s)  = q_g(s)","category":"page"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"Consequently, even though the ODE may be integrated to find a solution, an iterative scheme over the source term must be used to converge to the true solution. Commonly, each of these iterations is known as a transport sweep.","category":"page"},{"location":"moc.html#Discretization-scheme","page":"Method of Characteristics","title":"Discretization scheme","text":"","category":"section"},{"location":"moc.html","page":"Method of Characteristics","title":"Method of Characteristics","text":"In order to find a numerical solution to the transport equation we must transform the continuous operators and dependencies into discrete versions. Up to this point we have only applied a discretization scheme to the energy dependence. We now apply a discretization scheme to the other continuous dependencies: solid angle boldsymbolhatOmega and spatial variables.","category":"page"},{"location":"index.html#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"NeutronTransport is a high-performance library designed to achieve fast and advanced reactor physics calculations.","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Currently, there isn't any reactor physics program written in Julia. NeutronTransport main goal is to open that gate so people from the field get to know this amazing programming language. On the other hand, since it is free software (as in free speech), it may be developed collaboratively by volunteer engineers/computer programmers who do not want to be, as R. J. G. Stamm'ler, M. J. Abbate (1983) used to say, users that will never understand the programs they are to use and, as computer slaves, consider them as black boxes, blindly trusting their results.","category":"page"},{"location":"index.html#Getting-Started","page":"Introduction","title":"Getting Started","text":"","category":"section"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"The package can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"pkg> add NeutronTransport","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"Or, equivalently, via the Pkg API:","category":"page"},{"location":"index.html","page":"Introduction","title":"Introduction","text":"julia> import Pkg; Pkg.add(\"NeutronTransport\")","category":"page"}]
}
