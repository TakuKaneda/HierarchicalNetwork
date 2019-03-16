# HierarchicalNetwork

Solve the multi-stage stochastic linear programming for the hierarchical distribution network.

**REQUIRE** `Julia v1.0` and `Gurobi` solver with `JuMP` package

<!-- # Table of contents -->
1. [Data Generation](#data_generation)
2. [Method](#method)
  1. [Perfect Foresight](#perfectforesight)
  2. [Scenario-Based Robust MPC with Nested Dantzig-Wolfe](#sbrmpcDW)

## Data Generation <a name="data_generation"></a>
For the data generation (e.g. network topology, demand profile, lattice of PV etc..), please check the notebook [`test/data_generation.ipynb`](test/data_generation.ipynb). The generated data is saved in `data/Case1`.

## Method <a name="method"></a>
Here we present the description of methods for managing uncertainty.

### Perfect Foresight <a name="perfectforesight"></a>
`run_parallel_perfect_foresight.jl` implements parallelly the perfect foresight policy on the network to decide the interface flow values (`p_in` and `p_out`) which will be used for SDDP implementation. The result is saved as `JLD2` format.

Implementation is done for one branch of the medium voltage network.
I would recommend to set the number of samples, `NSamples`, greater than 500 in order to get 'plausible' data.

### Scenario-Based Robust MPC with Nested Dantzig-Wolfe <a name="sbrmpcDW"></a>
[`test/run_dantzig_wolfe.jl`](`test/run_dantzig_wolfe.jl`) implements the SBR-MPC with the nested DW to solve the problem. Again, we consider on one branch of the medium voltage network.
