NSamples = 500;
NProcessors = 4;
using Distributed
addprocs(NProcessors)  # add processors

@everywhere include("src/load_data.jl");
# here we consider only one branch!
@everywhere PGenerationMax = generator_capacity_rate .* feeder_capacity/NBranches;

# load model
@everywhere include("test/model_perfect_foresight.jl")

# create model
@everywhere PFProblem = PerfectForesightCreation()
Output = pmap(ImplementPerfectForesight,1:NSamples)
p_in_average, p_out_average = ComputePinPoutAverageParallel(Output)
## save as jld2 file
save("data/Case1/sddp_interface_data.jld2",
    "p_out_data",p_out_average,
    "p_in_data",p_in_average);
