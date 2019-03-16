# data load
include("../src/load_data.jl");

##################### Algorithm Hyperparameter #####################
NScenarios = 5;  # number of scenarios for SBR MPC with nested DW
K = 20; # max num of iteration for the Algorithm
FutureStepSize = 0; # size of looking a head
KnowCurrentOutcome = false; # wheather knowing the current outcome or not
####################################################################
# model load
include("model_dantzig_wolfe.jl");

# here we consider only one branch!
PGenerationMax = generator_capacity_rate .* feeder_capacity/NBranches; # scale the generator capacity

##################### Experiment Parameter #########################
NSamples = 100; # numeber of samples to test the Alogrithms
RealScenario = SamplePath(TransProb,NSamples); # Generate real scenarios
####################################################################

## Implementation
SBRTotalCost = zeros(NSamples,H);
SBRTotalStorage = zeros(NSamples,H);
SBROnlineTime = zeros(NSamples,H);
for i = 1:NSamples
    SBRTotalCost[i,:], SBRTotalStorage[i,:], SBROnlineTime[i,:] = ImplementMPC(SBRMPCSubproblem,RealScenario[:,i],FutureStepSize,KnowCurrentOutcome,i)
end
## data save in JLD2 format
save("trial_solution/N"* string(NSamples) * "_FS"*string(FutureStepSize) *"_nestedSBRDW_result.jld2","SBRTotalCost",SBRTotalCost,
    "SBRTotalStorage",SBRTotalStorage,"SBROnlineTime",SBROnlineTime,
    "NSamples",NSamples,"RealScenario",RealScenario,"FutureStepSize",FutureStepSize,
    "KnowCurrentOutcome",KnowCurrentOutcome);
####################################################################
####################################################################
# implement perfect foresight
include("model_perfect_foresight.jl")

PFTotalCost = zeros(NSamples);

## algorithm implementation
for i = 1:NSamples
    # fix balance_RHS to the Net demand values
    start_t = time();
    for sn=1:NSubnetworks[2], n = LVBatteryNodes, t = 1:H
        fix(PFProblem.balance_RHS[sn,n,t],  LVDemandDataDict[n][t] - (LVPVNodeDict[n] * NormPV[t][RealScenario[t,i]]));
    end
    # AssignTime[i] = time()-start_t;
    println("solve sample $i")
    start_t = time();
    optimize!(PFProblem.m)
    SolveTime[i] =  time()-start_t;
    println("   solve time: ", SolveTime[i])
    PFTotalCost[i] = objective_value(PFProblem.m);
end
@time PFProblem = PerfectForesightCreation()
save("trial_solution/N"* string(NSamples) *"_PF_result.jld2","PFTotalCost",PFTotalCost,
    "NSamples",NSamples,"RealScenario",RealScenario);
