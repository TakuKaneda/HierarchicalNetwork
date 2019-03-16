NSamples = 100

include("../src/load_data.jl");
# here we consider only one branch!
PGenerationMax = generator_capacity_rate .* feeder_capacity/NBranches;

# load model
include("model_perfect_foresight.jl")
## model creation
RealScenario  = SamplePath(TransProb, NSamples);
@time PFProblem = PerfectForesightCreation()
## algorithm result
AssignTime = zeros(NSamples);
SolveTime = zeros(NSamples);
TotalCost = zeros(NSamples);
p_out_Store = Array{Float64}(undef,  NSubnetworks[2], H, NSamples);

## algorithm implementation
for i = 1:NSamples
    # fix balance_RHS to the Net demand values
    start_t = time();
    for sn=1:NSubnetworks[2], n = LVBatteryNodes, t = 1:H
        fix(PFProblem.balance_RHS[sn,n,t],  LVDemandDataDict[n][t] - (LVPVNodeDict[n] * NormPV[t][RealScenario[t,i]]));
    end
    AssignTime[i] = time()-start_t;
    println("solve sample $i")
    start_t = time();
    optimize!(PFProblem.m)
    SolveTime[i] =  time()-start_t;
    p_out_Store[:,:,i] = value.(PFProblem.p_out);
    TotalCost[i] = objective_value(PFProblem.m);
end
## save the cost result
save("trial_solution/N"* string(NSamples) *"_PF_result.jld2","PFTotalCost",TotalCost,
    "NSamples",NSamples,"RealScenario",RealScenario);

## Compute the average p_in and p_out
p_in_average, p_out_average = ComputePinPoutAverage(RealScenario,p_out_Store,NSamples)
## save as jld2 file
# save("data/Case1/sddp_interface_data.jld2",
#     "p_out_data",p_out_average,
#     "p_in_data",p_in_average);
