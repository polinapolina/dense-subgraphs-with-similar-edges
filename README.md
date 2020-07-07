## dense-subgraphs-with-similar-edges
Implementations for the experiments included in the paper:

Polina Rozenshtein, Giulia Preti, Aristides Gionis, and Yannis Velegrakis. "Mining Dense Subgraphs with Similar Edges." ECML PKDD 2020

### Compile C code solving Parametric MinCut 
* `make -f makefile` compiles C library with parametric MinCut for DenSim
* `make -f makefile_baseline` compiles C library with parametric MinCut for the baselines

This C implementation of Parametric MinCut is based on the implementation from Chandran, Bala G., and Dorit S. Hochbaum. "A computational study of the pseudoflow and push-relabel algorithms for the maximum flow problem." Operations research 57.2 (2009): 358-376.)

