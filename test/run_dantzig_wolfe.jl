# data load
include("../src/load_data.jl");

##################### Algorithm Hyperparameter #####################
NScenarios = 5;  # number of scenarios for SBR MPC with nested DW
K = 10; # max num of iteration for the Algorithm
FutureStepSize = 1; # size of looking a head
KnowCurrentOutcome = true; # wheather knowing the current outcome or not
####################################################################
# model load
include("model_dantzig_wolfe.jl");

# here we consider only one branch!
PGenerationMax = generator_capacity_rate .* feeder_capacity/NBranches; # scale the generator capacity

##################### Experiment Parameter #########################
NSamples = 50; # numeber of samples to test the Alogrithms
RealScenario = SamplePath(TransProb,NSamples); # Generate real scenarios
####################################################################

## Implementation
SBRTotalCost = zeros(NSamples,H);
SBRTotalStorage = zeros(NSamples,H);
SBROnlineTime = zeros(NSamples,H);
for i = 1:NSamples
    SBRTotalCost[i,:], SBRTotalStorage[i,:], SBROnlineTime[i,:] = ImplementMPC(SBRMPCSubproblem,RealScenario[:,i],FutureStepSize,KnowCurrentOutcome,i)
end
## assign the early stage costs for 2:NSamples
for i = 2:NSamples
     TimeStages = 1:findall(NOutcomes.!=1)[1]-1;
     SBRTotalCost[i,TimeStages] = SBRTotalCost[1,TimeStages] ;
end
## data save in JLD2 format
save("trial/N"* string(NSamples) * "_FS"*string(FutureStepSize)*"_"*string(KnowCurrentOutcome)*"_nestedSBRDW_result.jld2","SBRTotalCost",SBRTotalCost,
    "SBRTotalStorage",SBRTotalStorage,"SBROnlineTime",SBROnlineTime,
    "NSamples",NSamples,"RealScenario",RealScenario,"FutureStepSize",FutureStepSize,
    "KnowCurrentOutcome",KnowCurrentOutcome);
