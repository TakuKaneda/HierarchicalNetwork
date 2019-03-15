using Pkg, JuMP, StatsBase, DataFrames, CSV, DelimitedFiles, FileIO, JLD2
# using LightGraphs, GraphPlot  # for plotting graph

include("load_function.jl");

## load network data
network_data = load("data/Case1/network_data.jld2");

LVLoadNodeDict = network_data["LVLoadNodeDict"];
LVChildrenNodes = network_data["LVChildrenNodes"];
LVDemandDataDict = network_data["LVDemandDataDict"];
LVPVNodeDict = network_data["LVPVNodeDict"];
MVChildrenNodes = network_data["MVChildrenNodes"];
MVLoadNodeDict = network_data["MVLoadNodeDict"];
feeder_capacity = network_data["feeder_capacity"];
NBranches = network_data["NBranches"];
# clear memory
network_data = nothing;
## load lattice data
lattice_data = load("data/Case1/lattice_data.jld2");

TransProb = lattice_data["TransProb"];
NormPV = lattice_data["NormPV"];
NOutcomes = lattice_data["NOutcomes"];
# clear memory
lattice_data = nothing;

## load generator data
generator_data = load("data/Case1/generator_data.jld2")
generator_capacity_rate = generator_data["generator_capacity_rate"]
MargCost = generator_data["MargCost"];
generator_data = nothing;

## load battery data (Tesla Powerwall 2)
battery_data = load("data/Case1/battery_data.jld2")
BatteryCapacity = battery_data["BatteryCapacity"];
BatteryChargeEfficiency = battery_data["BatteryChargeEfficiency"];
BatteryDischargeEfficiency = battery_data["BatteryDischargeEfficiency"];
BatteryChargeRate = battery_data["BatteryChargeRate"];
InitialStorage = battery_data["InitialStorage"];
battery_data = nothing;

## node & line data
NLayers = 2;
MVNNodes = length(MVChildrenNodes);
MVNLines = MVNNodes * 2 - 1;  # CAUTION! this needs to be generalized
LVNNodes = length(LVChildrenNodes);
LVNLines = LVNNodes; # there is no lower layer and the root is connect to MV

LVDemandNodes = collect(keys(LVDemandDataDict));
LVPVNodes = collect(keys(LVPVNodeDict));
LVBatteryNodes = collect(keys(LVPVNodeDict));
LVNonPVNodes = setdiff(LVDemandNodes,LVBatteryNodes);
LVNonBatteryNodes = setdiff(LVDemandNodes, LVBatteryNodes) # node wihtout PV and Battery
LVNonDemandNodes = setdiff(1:LVNNodes, collect(keys(LVDemandDataDict)))
NGenerators = length(generator_capacity_rate);

NSubnetworks = [1, MVNNodes-1];
NLayerNodes = [MVNNodes, LVNNodes];
NLayerLines = [MVNLines, LVNLines];
## Other parameters
PGenerationMax = generator_capacity_rate .* feeder_capacity
PGenerationMin = zeros(NGenerators);
VOLL = 5.0;
H = size(NormPV,1);
T = 1:H;
SLimit = 10000;
