# data load
include("../src/load_data.jl");

##################### Algorithm Hyperparameter #####################
NScenarios = 4;  # number of scenarios for SBR MPC with nested DW
K = 10; # max num of iteration for the Algorithm
FutureStepSize = 0; # size of looking a head
KnowCurrentOutcome = true; # wheather knowing the current outcome CAUTION: bug in false
####################################################################
##################### Experiment Parameter #########################
NSamples = 50; # numeber of samples to test the Alogrithms
RealScenario = SamplePath(TransProb,NSamples); # Generate real scenarios
####################################################################

# model load
include("model_dantzig_wolfe.jl");

# here we consider only one branch!
PGenerationMax = generator_capacity_rate .* feeder_capacity/NBranches; # scale the generator capacity

## Implementation
SBRSolutions = Array{Any}(undef,NSamples);
for i = 1:NSamples
    SBRSolutions[i] = ImplementMPC(SBRMPCSubproblem,RealScenario[:,i],FutureStepSize,KnowCurrentOutcome,i)
end
## assign the early-stage values for 2:NSamples
if NSamples > 1
    for i = 2:NSamples
         TimeStages = 1:findall(NOutcomes.!=1)[1] - FutureStepSize - 1;
         SBRSolutions[i].TotalCost[TimeStages] = SBRSolutions[1].TotalCost[TimeStages] ;
         SBRSolutions[i].TotalProduction[TimeStages] = SBRSolutions[1].TotalProduction[TimeStages] ;
         SBRSolutions[i].TotalStorage[TimeStages] = SBRSolutions[1].TotalStorage[TimeStages] ;
         SBRSolutions[i].TotalLoadShedding[TimeStages] = SBRSolutions[1].TotalLoadShedding[TimeStages] ;
         SBRSolutions[i].OnlineTime[TimeStages] = SBRSolutions[1].OnlineTime[TimeStages] ;
         SBRSolutions[i].IterationCount[TimeStages] = SBRSolutions[1].IterationCount[TimeStages] ;
    end
end
## data save in JLD2 format
# save("trial/N"* string(NSamples) * "_FS"*string(FutureStepSize)*"_"*string(KnowCurrentOutcome)*"_SBRnestedDW_result.jld2",
#     "SBRSolutions",SBRSolutions,
#     "NSamples",NSamples,"RealScenario",RealScenario,"FutureStepSize",FutureStepSize,
#     "KnowCurrentOutcome",KnowCurrentOutcome);
