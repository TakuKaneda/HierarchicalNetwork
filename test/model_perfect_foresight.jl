using JuMP, Clp,GLPK,Ipopt, Statistics, Cbc
# using Gurobi

struct MyProblem
    m
    pflow
    pgeneration
    storage
    batterycharge
    batterydischarge
    loadshedding
    productionshedding
    # p_in
    p_out
    balance_RHS

    BatteryDynamics_stage1
    BatteryDynamics
    FlowMax
    FlowMin
    StorageMax
    BatteryChargeMax
    BatteryDischargeMax
    # Pin_Flow_equality
    # Pin_Pout_equality
    Pout_Flow_equality
    InterfaceFlow_equality
    RootFlowZero
    MVBalance_rootnode
    MVBalance
    LVBalance_PVBattery
    LVBalance_NonPV
    LVBalance_NonDemand
    GenerationMax
    GenerationMin
end
##

mutable struct Solutions # change to immutable so as to be converted to SharedArray
    "struct that stores solutions over a sample"
    pflow#::Array{Float64,2}
    pgeneration#::Array{Float64,2}
    storage#::Array{Float64,2}
    batterycharge#::Array{Float64,2}
    batterydischarge#::Array{Float64,2}
    loadshedding#::Array{Float64,2}
    productionshedding#::Array{Float64,2}
    p_in#::Array{Float64,1}
    p_out#::Array{Float64,1}
    StageCost#::Array{Float64,1}
    ObjectiveValue

    # constructor
    # maybe there is a better way to assign default values
    Solutions() = new(
        Any,
        Any,
        Any,
        Any,
        Any,
        Any,
        Any,
        Any,
        Any,
        Any,
        Any
    )
end
##
function PerfectForesightCreation()
    "
    create a struct of perfect foresight policy problem
    "
    ## Build model
    # m = Model(with_optimizer(Gurobi.Optimizer, LogToConsole=0));
    # m = Model(with_optimizer(Clp.Optimizer, LogLevel=1));
    # m = Model(with_optimizer(GLPK.Optimizer, msg_lev=0));
    # m = Model(with_optimizer(Ipopt.Optimizer, print_level=0));
    m = Model(with_optimizer(Cbc.Optimizer, loglevel=0));
    ## Variables
    @variable(m, pflow[l = 1:NLayers,sn =1:NSubnetworks[l],i=1:NLayerLines[l],1:H]);
    @variable(m, pgeneration[1:NGenerators,1:H]);
    @variable(m, storage[sn =1:NSubnetworks[2],n in LVBatteryNodes,1:H] >= 0);
    @variable(m, batterycharge[sn =1:NSubnetworks[2],n in LVBatteryNodes,1:H] >= 0);
    @variable(m, batterydischarge[sn =1:NSubnetworks[2],n in LVBatteryNodes,1:H] >= 0);
    @variable(m, loadshedding[l = 1:NLayers,sn =1:NSubnetworks[l],n=1:NLayerNodes[l],1:H] >= 0);
    @variable(m, productionshedding[l = 1:NLayers,sn =1:NSubnetworks[l],n=1:NLayerNodes[l],1:H] >= 0);
    # @variable(m, p_in[l=1:NLayers,sn=1:NSubnetworks[l],1:H]);
    @variable(m, p_out[sn=1:NSubnetworks[2],1:H]);
    @variable(m, balance_RHS[sn =1:NSubnetworks[2], n in LVBatteryNodes,1:H]);

    ## Objective - minimize cost of generation and load shedding
    # we need to devide by 4 since 15-min model
    @objective(m, Min,
        1/4 * (sum(MargCost[g]*pgeneration[g,t] for g=1:NGenerators, t=1:H)
        + VOLL * sum(loadshedding))
    );

    ## Constraints
    # dynamics
    @constraint(m, BatteryDynamics_stage1[sn =1:NSubnetworks[2],n in LVBatteryNodes],
         (storage[sn,n,1] - InitialStorage
         - BatteryChargeEfficiency * batterycharge[sn,n,1] / 4  # 15-min time step
         + batterydischarge[sn,n,1]/BatteryDischargeEfficiency / 4 # 15-min time step
         == 0)
    );
    @constraint(m, BatteryDynamics[sn =1:NSubnetworks[2],n in LVBatteryNodes, t=2:H],
        (storage[sn,n,t] - storage[sn,n,t-1]
         - BatteryChargeEfficiency * batterycharge[sn,n,t] / 4
         + batterydischarge[sn,n,t]/BatteryDischargeEfficiency / 4
        == 0)
    );

    # Flow Limits
    @constraint(m, FlowMax[l = 1:NLayers,sn =1:NSubnetworks[l],i=1:NLayerLines[l],t = 1:H],
        (pflow[l,sn,i,t] <= SLimit[l])
    );
    @constraint(m, FlowMin[l = 1:NLayers,sn =1:NSubnetworks[l],i=1:NLayerLines[l],t = 1:H],
        ( - pflow[l,sn,i,t] <= SLimit[l])
    );

    # Storage Capacity
    @constraint(m, StorageMax[sn =1:NSubnetworks[2],n in LVBatteryNodes, t=1:H],
        (storage[sn,n,t] <= BatteryCapacity)
    );

    # Charging Capacity
    @constraint(m, BatteryChargeMax[sn =1:NSubnetworks[2],n in LVBatteryNodes, t = 1:H],
        (batterycharge[sn,n,t] <= BatteryChargeRate)
    );

    # Discharging Capacity
    @constraint(m, BatteryDischargeMax[sn =1:NSubnetworks[2],n in LVBatteryNodes, t = 1:H],
        (batterydischarge[sn,n,t] <= BatteryChargeRate)
    );

    # p_in & pflow equality
    # @constraint(m, Pin_Flow_equality[sn=1:NSubnetworks[2],t = 1:H],
    #     (p_in[2,sn,t] - pflow[1,1,sn+NLayerNodes[1],t] == 0)
    # );
    # # p_in & p_out equality
    # @constraint(m, Pin_Pout_equality[sn=1:NSubnetworks[2],t = 1:H],
    #     (p_in[2,sn,t] - p_out[sn,t] == 0)
    # );
    @constraint(m, Pout_Flow_equality[sn=1:NSubnetworks[2],t = 1:H],
        (p_out[sn,t] - pflow[1,1,sn+NLayerNodes[1],t] == 0)
    );
    # pflow of layer 1 and layer 2 equality
    @constraint(m, InterfaceFlow_equality[sn=1:NSubnetworks[2],t = 1:H],
        (pflow[1,1,sn+NLayerNodes[1],t] - pflow[2,sn,1,t] == 0)
    );
    # pflow of the root is zero
    @constraint(m,RootFlowZero[t=1:H],pflow[1,1,1,t]==0);
    # MV Balancing - root node
    @constraint(m, MVBalance_rootnode[t = 1:H],
        (sum(pgeneration[g,t] for g = 1:NGenerators)
        + loadshedding[1,1,1,t] - productionshedding[1,1,1,t]
        + sum(pflow[1,1,m,t] for m in MVChildrenNodes[1])
        - pflow[1,1,1,t]
        == 0.0 # no demand
        )
    );
    # MV Balancing - usual nodes (transformers)
    @constraint(m, MVBalance[n = 2:NLayerNodes[1],t = 1:H],
        (loadshedding[1,1,n,t] - productionshedding[1,1,n,t]
        + sum(pflow[1,1,m,t] for m in MVChildrenNodes[n])
        - pflow[1,1,n,t]
        == 0.0 # no demand
        )
    );
    # LV Balancing - PV & Battery houses
    @constraint(m, LVBalance_PVBattery[sn=1:NSubnetworks[2], n = LVBatteryNodes, t = 1:H],
        (batterydischarge[sn,n,t] - batterycharge[sn,n,t]
        + loadshedding[2,sn,n,t] - productionshedding[2,sn,n,t]
        + sum(pflow[2,sn,m,t] for m in LVChildrenNodes[n])
        - pflow[2,sn,n,t]
        ==  balance_RHS[sn,n,t]
        #= LVDemandDataDict[n][t] # demand data
        - (LVPVNodeDict[n] * NormPV[t][RealPath[t]]) =# # pv output scaled by capacity
        )
    );

    # LV Balancing - no PV and Battery houses
    @constraint(m, LVBalance_NonPV[sn=1:NSubnetworks[2], n = LVNonPVNodes, t = 1:H],
        (loadshedding[2,sn,n,t] - productionshedding[2,sn,n,t]
        + sum(pflow[2,sn,m,t] for m in LVChildrenNodes[n])
        - pflow[2,sn,n,t]
        == LVDemandDataDict[n][t] # only demand data
        )
    );
    # LV Balancing - non-house nodes
    @constraint(m, LVBalance_NonDemand[sn=1:NSubnetworks[2], n = LVNonDemandNodes, t = 1:H],
        (loadshedding[2,sn,n,t] - productionshedding[2,sn,n,t]
        + sum(pflow[2,sn,m,t] for m in LVChildrenNodes[n])
        - pflow[2,sn,n,t]
        == 0.0 # only demand data
        )
    );

    # Generation Limits
    @constraint(m, GenerationMax[g = 1:NGenerators, t=1:H],
        (pgeneration[g,t] <= PGenerationMax[g])
    );
    @constraint(m, GenerationMin[g = 1:NGenerators, t=1:H],
        ( - pgeneration[g,t] <= - PGenerationMin[g])
    );
    return MyProblem(m,pflow,pgeneration,storage,batterycharge,batterydischarge,
    loadshedding,productionshedding,#=p_in,=#p_out,balance_RHS,
    BatteryDynamics_stage1,BatteryDynamics,FlowMax,
    FlowMin ,StorageMax,BatteryChargeMax,BatteryDischargeMax #=,Pin_Flow_equality,
    Pin_Pout_equality=#,Pout_Flow_equality,InterfaceFlow_equality,RootFlowZero,MVBalance_rootnode,
    MVBalance,LVBalance_PVBattery,LVBalance_NonPV,LVBalance_NonDemand,GenerationMax,GenerationMin)
end
##
function ComputePinPoutAverage(SampleScenario,p_out_Store,NSamples)
    p_out_sum = [zeros(NSubnetworks[2],NOutcomes[t]) for t = 1:H];
    p_out_average = [zeros(NSubnetworks[2],NOutcomes[t]) for t = 1:H];
    p_out_count = [zeros(NOutcomes[t]) for t = 1:H];
    p_in_sum = [zeros(NOutcomes[t]) for t = 1:H];
    p_in_average = [zeros(NOutcomes[t]) for t = 1:H];
    p_in_count = [zeros(NOutcomes[t]) for t = 1:H];
    for i = 1:NSamples, t = 1:H
        p_out_sum[t][:,SampleScenario[t,i]] += p_out_Store[:,t,i]
        p_out_count[t][SampleScenario[t,i]] += 1
        # I use median instead of mean for p_in
        p_in_sum[t][SampleScenario[t,i]] += median(p_out_Store[:,t,i])
        p_in_count[t][SampleScenario[t,i]] += 1
    end
    for t = 1:H, k = 1:NOutcomes[t]
        p_out_average[t][:,k] = p_out_sum[t][:,k]/p_out_count[t][k];
        p_in_average[t][k] = p_in_sum[t][k]/p_in_count[t][k];
    end
    return p_in_average, p_out_average
end
##
function ComputePinPoutAverageParallel(Output)
    NSamples = length(Output)
    p_out_Store = Array{Float64}(undef,  NSubnetworks[2], H, NSamples);
    SampleScenario  = ones(Int64,H,NSamples);
    for i = 1:NSamples
        p_out_Store[:,:,i] = Output[i][1];
        SampleScenario[:,i] = Output[i][2];
    end
    p_out_sum = [zeros(NSubnetworks[2],NOutcomes[t]) for t = 1:H];
    p_out_average = [zeros(NSubnetworks[2],NOutcomes[t]) for t = 1:H];
    p_out_count = [zeros(NOutcomes[t]) for t = 1:H];
    p_in_sum = [zeros(NOutcomes[t]) for t = 1:H];
    p_in_average = [zeros(NOutcomes[t]) for t = 1:H];
    p_in_count = [zeros(NOutcomes[t]) for t = 1:H];
    for i = 1:NSamples, t = 1:H
        p_out_sum[t][:,SampleScenario[t,i]] += p_out_Store[:,t,i]
        p_out_count[t][SampleScenario[t,i]] += 1
        # I use median instead of mean for p_in
        p_in_sum[t][SampleScenario[t,i]] += median(p_out_Store[:,t,i])
        p_in_count[t][SampleScenario[t,i]] += 1
    end
    for t = 1:H, k = 1:NOutcomes[t]
        p_out_average[t][:,k] = p_out_sum[t][:,k]/p_out_count[t][k];
        p_in_average[t][k] = p_in_sum[t][k]/p_in_count[t][k];
    end
    return p_in_average, p_out_average
end

##
function ImplementPerfectForesight(idx)
    OnePath = SamplePath(TransProb);
    for sn=1:NSubnetworks[2], n = LVBatteryNodes, t = 1:H
        fix(PFProblem.balance_RHS[sn,n,t],  LVDemandDataDict[n][t] - (LVPVNodeDict[n] * NormPV[t][OnePath[t]]));
    end
    optimize!(PFProblem.m)
    p_out_Store = value.(PFProblem.p_out);
    return p_out_Store, OnePath
end
##
function PerfectForesight(RealPath, solutions, idx=NaN)
    "
    implement perfect foresight policy
    Args:
        - RealPath: Real path of the 'current' sample
        - solutions: an instance of the struct Solutions
    "
    ## Build model
    m = Model(with_optimizer(Gurobi.Optimizer, LogToConsole=1));

    ## Variables
    @variable(m, pflow[l = 1:NLayers,sn =1:NSubnetworks[l],i=1:NLayerLines[l],1:H]);
    @variable(m, pgeneration[1:NGenerators,1:H]);
    @variable(m, storage[sn =1:NSubnetworks[2],n in LVBatteryNodes,1:H] >= 0);
    @variable(m, batterycharge[sn =1:NSubnetworks[2],n in LVBatteryNodes,1:H] >= 0);
    @variable(m, batterydischarge[sn =1:NSubnetworks[2],n in LVBatteryNodes,1:H] >= 0);
    @variable(m, loadshedding[l = 1:NLayers,sn =1:NSubnetworks[l],n=1:NLayerNodes[l],1:H] >= 0);
    @variable(m, productionshedding[l = 1:NLayers,sn =1:NSubnetworks[l],n=1:NLayerNodes[l],1:H] >= 0);
    @variable(m, p_in[l=1:NLayers,sn=1:NSubnetworks[l],1:H]);
    @variable(m, p_out[sn=1:NSubnetworks[2],1:H]);

    ## Objective - minimize cost of generation and load shedding
    # we need to devide by 4 since 15-min model
    @objective(m, Min,
        1/4 * (sum(MargCost[g]*pgeneration[g,t] for g=1:NGenerators, t=1:H)
        + VOLL * sum(loadshedding))
    );

    ## Constraints
    # dynamics
    @constraint(m, BatteryDynamics_stage1[sn =1:NSubnetworks[2],n in LVBatteryNodes],
         (storage[sn,n,1] - InitialStorage
         - BatteryChargeEfficiency * batterycharge[sn,n,1] / 4  # 15-min time step
         + batterydischarge[sn,n,1]/BatteryDischargeEfficiency / 4 # 15-min time step
         == 0)
    );
    @constraint(m, BatteryDynamics[sn =1:NSubnetworks[2],n in LVBatteryNodes, t=2:H],
        (storage[sn,n,t] - storage[sn,n,t-1]
         - BatteryChargeEfficiency * batterycharge[sn,n,t] / 4
         + batterydischarge[sn,n,t]/BatteryDischargeEfficiency / 4
        == 0)
    );

    # Flow Limits
    @constraint(m, FlowMax[l = 1:NLayers,sn =1:NSubnetworks[l],i=1:NLayerLines[l],t = 1:H],
        (pflow[l,sn,i,t] <= SLimit[l])
    );
    @constraint(m, FlowMin[l = 1:NLayers,sn =1:NSubnetworks[l],i=1:NLayerLines[l],t = 1:H],
        ( - pflow[l,sn,i,t] <= SLimit[l])
    );

    # Storage Capacity
    @constraint(m, StorageMax[sn =1:NSubnetworks[2],n in LVBatteryNodes, t=1:H],
        (storage[sn,n,t] <= BatteryCapacity)
    );

    # Charging Capacity
    @constraint(m, BatteryChargeMax[sn =1:NSubnetworks[2],n in LVBatteryNodes, t = 1:H],
        (batterycharge[sn,n,t] <= BatteryChargeRate)
    );

    # Discharging Capacity
    @constraint(m, BatteryDischargeMax[sn =1:NSubnetworks[2],n in LVBatteryNodes, t = 1:H],
        (batterydischarge[sn,n,t] <= BatteryChargeRate)
    );

    # p_in & pflow equality
    @constraint(m, Pin_Flow_equality[sn=1:NSubnetworks[2],t = 1:H],
        (p_in[2,sn,t] - pflow[1,1,sn+NLayerNodes[1],t] == 0)
    );
    # p_in & p_out equality
    @constraint(m, Pin_Pout_equality[sn=1:NSubnetworks[2],t = 1:H],
        (p_in[2,sn,t] - p_out[sn,t] == 0)
    );
    # pflow of layer 1 and layer 2 equality
    @constraint(m, InterfaceFlow_equality[sn=1:NSubnetworks[2],t = 1:H],
        (pflow[1,1,sn+NLayerNodes[1],t] - pflow[2,sn,1,t] == 0)
    );
    # pflow of the root is zero
    @constraint(m,RootFlowZero[t=1:H],pflow[1,1,1,t]==0);
    # MV Balancing - root node
    @constraint(m, MVBalance_rootnode[t = 1:H],
        (sum(pgeneration[g,t] for g = 1:NGenerators)
        + loadshedding[1,1,1,t] - productionshedding[1,1,1,t]
        + sum(pflow[1,1,m,t] for m in MVChildrenNodes[1])
        - pflow[1,1,1,t]
        == 0.0 # no demand
        )
    );
    # MV Balancing - usual nodes (transformers)
    @constraint(m, MVBalance[n = 2:NLayerNodes[1],t = 1:H],
        (loadshedding[1,1,n,t] - productionshedding[1,1,n,t]
        + sum(pflow[1,1,m,t] for m in MVChildrenNodes[n])
        - pflow[1,1,n,t]
        == 0.0 # no demand
        )
    );
    # LV Balancing - PV & Battery houses
    @constraint(m, LVBalance_PVBattery[sn=1:NSubnetworks[2], n = LVBatteryNodes, t = 1:H],
        (batterydischarge[sn,n,t] - batterycharge[sn,n,t]
        + loadshedding[2,sn,n,t] - productionshedding[2,sn,n,t]
        + sum(pflow[2,sn,m,t] for m in LVChildrenNodes[n])
        - pflow[2,sn,n,t]
        == LVDemandDataDict[n][t] # demand data
        - (LVPVNodeDict[n] * NormPV[t][RealPath[t]]) # pv output scaled by capacity
        )
    );

    # LV Balancing - no PV and Battery houses
    @constraint(m, LVBalance_NonPV[sn=1:NSubnetworks[2], n = LVNonPVNodes, t = 1:H],
        (loadshedding[2,sn,n,t] - productionshedding[2,sn,n,t]
        + sum(pflow[2,sn,m,t] for m in LVChildrenNodes[n])
        - pflow[2,sn,n,t]
        == LVDemandDataDict[n][t] # only demand data
        )
    );
    # LV Balancing - non-house nodes
    @constraint(m, LVBalance_NonDemand[sn=1:NSubnetworks[2], n = LVNonDemandNodes, t = 1:H],
        (loadshedding[2,sn,n,t] - productionshedding[2,sn,n,t]
        + sum(pflow[2,sn,m,t] for m in LVChildrenNodes[n])
        - pflow[2,sn,n,t]
        == 0.0 # only demand data
        )
    );

    # Generation Limits
    @constraint(m, GenerationMax[g = 1:NGenerators, t=1:H],
        (pgeneration[g,t] <= PGenerationMax[g])
    );
    @constraint(m, GenerationMin[g = 1:NGenerators, t=1:H],
        ( - pgeneration[g,t] <= - PGenerationMin[g])
    );
    ## Solve
    # @time status = solve(m);
    optimize!(m);
    # println("Objective value: ", getobjectivevalue(m))

    ## Store Results
    solutions.pflow = value.(pflow);
    solutions.pgeneration = value.(pgeneration);
    solutions.storage = value.(storage);
    solutions.batterycharge = value.(batterycharge);
    solutions.batterydischarge = value.(batterydischarge);
    solutions.loadshedding = value.(loadshedding);
    solutions.productionshedding = value.(productionshedding);
    solutions.p_in = value.(p_in);
    solutions.p_out = value.(p_out);
    solutions.ObjectiveValue = objective_value(m);
    solutions.StageCost = [1/4*(sum(MargCost[g]*value.(pgeneration[g,t]) for g = 1:NGenerators)
            +VOLL*sum(sum(sum(value.(loadshedding[l,sn,n,t]) for n=1:NLayerNodes[l]) for sn =1:NSubnetworks[l]) for  l = 1:NLayers)) for t=1:H]
    if isnan(idx)
        println("Total cost is $(sum(solutions.StageCost[:]))")
    else
        println("Sample $idx, Total cost is $(sum(solutions.StageCost[:]))")
    end
    return #solutions.StageCost[:]
end
