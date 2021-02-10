# Orchestrated Trios: Compiling for Efficient Communication in Quantum Programs with 3-Qubit Gates

This is the source code repository for the paper
[Orchestrated Trios: Compiling for Efficient Communication in Quantum Programs with 3-Qubit Gates][paper]
to be published in the proceedings of ASPLOS '21, the 26th International Conference on Architectural Support for Programming Languages and Operating Systems, April 2021.

The benchmarks used can be found [here](https://github.com/jmbaker94/quantumcircuitbenchmarks).

[paper]: https://caseyduckering.com/#orchestrated-trios


### Dependencies

Python dependencies:
- Cirq==0.8, Qiskit==0.14.0
- networkx
- [quantumcircuitbenchmarks@7d5452d](https://github.com/jmbaker94/quantumcircuitbenchmarks/tree/7d5452d3080cca41c4525c8134eb7fad773e41cb)
- numpy, matplotlib
- latextools

Julia 1.4/1.5 dependencies:
- LightGraphs, MetaGraphs, GraphPlot
- OrderedCollections, LinearAlgebra, Statistics
- PyCall
- PyPlot
- [cduck/Tweaks.jl@c9dc3fd](https://github.com/cduck/Tweaks.jl/tree/c9dc3fdb866d5a3a05c23e9bfff1839c5565d6d1)
