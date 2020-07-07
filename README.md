## dense-subgraphs-with-similar-edges
Implementations for the experiments included in the paper:

Polina Rozenshtein, Giulia Preti, Aristides Gionis, and Yannis Velegrakis. "Mining Dense Subgraphs with Similar Edges." ECML PKDD 2020

### Compile C code solving Parametric MinCut 
* `make -f makefile` compiles C library with parametric MinCut for DenSim
* `make -f makefile_baseline` compiles C library with parametric MinCut for the baselines

This C implementation of Parametric MinCut is based on the implementation from Chandran, Bala G., and Dorit S. Hochbaum. "A computational study of the pseudoflow and push-relabel algorithms for the maximum flow problem." Operations research 57.2 (2009): 358-376.)

### main code in Python 3.6
* list of dependencies requirements.txt
* installation: pip3 install -r requirements.txt

### running DenSim
* run ``construct_metagraph.py`` to create a metagraph for MinCut input
* run ``main.py`` on the constructed metagraph
* use ``-h`` for arguments

### running baselines
* run ``construct_baseline.py`` to create metagraphs
* run ``main_baseline.py`` on the constructed metagraph
* use ``-h`` for arguments


