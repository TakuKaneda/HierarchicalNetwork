# HierarchicalNetwork

Solve the multi-stage stochastic linear programming for the hierarchical distribution network.

**REQUIRE** `Julia v1.0` and `Gurobi` solver with `JuMP` package

`run_parallel_perfect_foresight.jl` implements parallelly the perfect foresight policy on the network to decide the interface flow values (`p_in` and `p_out`) which will be used for SDDP implementation. The result is saved as `JLD2` format.


Implementation is done for one branch of the medium voltage network.
I would recommend to set the number of samples, `NSamples`, greater than 500 in order to get 'plausible' data.

For the data generation, please check the notebook [`test/data_generation.ipynb`](test/data_generation.ipynb). The generated data is saved in `data/Case1`.
