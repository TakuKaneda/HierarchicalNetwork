## model for nested Dantzig Wolfe algorithm
## there are generators at the root node
using JuMP, Clp
# using Gurobi

ImplementedSolution = Array{Any}(undef, NLayers,NSubnetworks[NLayers],H)  # solutions actually implemented
TemporalSolutionSet = Array{Any}(undef, NLayers,NSubnetworks[NLayers],K+1)  # solution for the Algorithm

mutable struct MPCSolutionStruct
    "
    struct that stores a solution of a subproblem
    "
    # vars
    pflow
    pgeneration
    loadshedding
    productionshedding
    storage
    batterycharge
    batterydischarge
    p_in
    p_out
    lambda
    costupperbound
    SubnetworkCost # cost of the current layer
    ObjectiveValue # objective value
    CostForUpperLayer # candidate cost for the upper layer
    # dual
    Coupling
    LambdaSum
    # constructor
    MPCSolutionStruct() = new(
        # vars
        Any, # pflow
        Any, # pgeneration
        Any, # loadshedding
        Any, # productionshedding
        Any, # storage
        Any, # batterycharge
        Any, # batterydischarge
        Any, # p_in
        Any, # p_out
        Any, # lambda
        Any, # costupperbound
        Any, # SubnetworkCost
        Any, # ObjectiveValue
        Any, # CostForUpperLayer
        # dual
        Any, # Coupling
        Any # LambdaSum
        )
end
##
mutable struct SBRMPCProblemStruct
    "
    struct that stores a solution of a subproblem
    "
    # model
    m
    # vars
    pflow
    pgeneration
    loadshedding
    productionshedding
    storage
    batterycharge
    batterydischarge
    p_in
    p_out
    lambda
    costupperbound
    balance_RHS
    # constraints
    CostBound
    FixLoadShedding_Demand
    FixLoadShedding_NonDemand
    FixPin2zero
    Balance_Root
    Balance_Battery
    Balance_NonBattery
    BatteryDynamics_Stage1
    BatteryDynamics_CurrentStage
    BatteryDynamics
    FlowMax
    FlowMin
    GenerationMax
    GenerationMin
    ProductionSheddingMax
    StorageMax
    BatteryChargeMax
    BatteryDischargeMax
    Pin_Flow_equality
    Pout_Flow_equality
    Equal_pgeneration
    Equal_loadshedding
    Equal_productionshedding
    Equal_storage
    Equal_batterycharge
    Equal_batterydischarge
    Equal_pflow
    Coupling
    LambdaSum
    FixInterfaceFlow
    FixLambda
end
##
mutable struct MPCResultStruct
    "struct to store the result of one-day MPC implementation"
    TotalCost
    TotalProduction
    TotalStorage
    TotalLoadShedding
    OnlineTime
    ImplementedSolution
end
## definition of subproblems for Dantzig Wolfe
function SBRMPCSubproblem(NetworkID, MPCScenarios, TimeChoice, Iteration, dual_pi = zeros(H), LastIteration=false, Last_p_in = Any, Last_lambda = Any)
    """
    solve a subproblem at (l,sn) with inputs of the upper layer
    args    - NetworkID: ID of Subproblem [layerID, subnetID]
            - MPCScenarios[t,s]: outcome at stage t of scenario s
            - TimeChoice
            - Iteration: iteration upto now
            - dual_pi: dual variable of coupling constraint from the upper layer master
            - LastIteration: wheather the iteration is the last one or not
            - Last_p_in: if we are at the last iteration, Last_p_in will be the suggestion from the upper layer
            - Last_lambda: Last_lambda will be the suggestion from the upper layer
    return  - sol: the struct Solution_Multi which has the solution
    """
    CLayer = NetworkID[1]; # current layer
    CSubnet = NetworkID[2]; # current subnetwork

    # only consider the unique samples
    MPCScenarios = unique(MPCScenarios,dims=2);
    # number of scenarios
    MPCScenarioSize = size(MPCScenarios,2);

    ## model
    # m = Model(with_optimizer(Gurobi.Optimizer, LogToConsole=0));
    m = Model(with_optimizer(Clp.Optimizer, LogLevel=0));
    ## size of step: stage optimization
    TimeEnd = TimeChoice + size(MPCScenarios,1) - 1;

    ## Variables
    @variable(m, loadshedding[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize] >= 0)
    @variable(m, productionshedding[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize] >= 0)
    @variable(m, storage[n in LayerBatteryNodes[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize] >= 0)
    @variable(m, batterycharge[n in LayerBatteryNodes[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize] >= 0)
    @variable(m, batterydischarge[n in LayerBatteryNodes[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize] >= 0)
    @variable(m, pflow[1:NLayerLines[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize])
    if CLayer == 1
        @variable(m, p_out[1: NSubnetworks[CLayer+1], TimeChoice:TimeEnd])  # independent on scenarios
        @variable(m, lambda[1:Iteration] >=0)
        @variable(m, pgeneration[1:NGenerators, TimeChoice:TimeEnd, 1:MPCScenarioSize]) # only the root node of the root layer has generators
    elseif CLayer == NLayers
        @variable(m, p_in[TimeChoice:TimeEnd])  # independent on scenarios
    else
        @variable(m, p_out[1: NSubnetworks[CLayer+1], TimeChoice:TimeEnd]) # independent on scenarios
        @variable(m, p_in[TimeChoice:TimeEnd]) # independent on scenarios
        @variable(m, lambda[1:Iteration] >=0)
    end
    @variable(m, costupperbound)
    @variable(m, balance_RHS[n in 1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:MPCScenarioSize]);

    # find the children sub-networks under the current sub-network
    if CLayer < NLayers
        ChildrenNetworks = [NSubnetworks[CLayer+1]*(CSubnet-1)+j for j = 1:NSubnetworks[CLayer+1]] # set of children networks
    end
    ## objective
    if CLayer == 1 # root layer: not the reduced cost
        @objective(m, Min,
            (costupperbound + sum(lambda[j] * sum(TemporalSolutionSet[CLayer+1,sn,j].CostForUpperLayer for sn in ChildrenNetworks)　for j = 1:Iteration))
        );
    elseif CLayer == NLayers # bottom layer: no lower layers
        @objective(m, Min,
            (costupperbound
             - sum((dual_pi[t]) * (-1) * p_in[t] for t = TimeChoice:TimeEnd))
        );
    else # inner layers: with reduced cost & lower layers
        @objective(m, Min,
            (costupperbound - sum((dual_pi[t]) * (-1) * p_in[t] for t = TimeChoice:TimeEnd) +
            sum(lambda[j] * sum(TemporalSolutionSet[CLayer+1,sn,j].CostForUpperLayer for sn in ChildrenNetworks) for j = 1:Iteration))
        );
    end

    ## constraints
    # For the cost upper bound
    VOLP = 1e-4; # cost for productionshedding (but not include in the actual cost)
    if CLayer == 1
        @constraint(m, CostBound[s = 1:MPCScenarioSize],
            (1/4*(sum(MargCost[g] * pgeneration[g,t,s] for g = 1:NGenerators, t = TimeChoice:TimeEnd)
            + sum(VOLL * loadshedding[n,t,s] for n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd))
            + sum(VOLP * productionshedding[n,t,s] for n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd)
            <= costupperbound)
        )
    else
        @constraint(m, CostBound[s = 1:MPCScenarioSize],
            (1/4*sum(VOLL * loadshedding[n,t,s] for n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd)
            + sum(VOLP * productionshedding[n,t,s] for n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd)
             <= costupperbound)
        )
    end
    # For the first iteration (initial solution)
    # Set loadshedding to be the demend (because it is always feasible)
    # Also Set p_in to be zero (always feasible)
    if Iteration == 1
        @constraint(m, FixLoadShedding_Demand[n in LayerDemandNodes[CLayer], t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            loadshedding[n,t,s] == DemandDataDict[CLayer][n][t]
        )
        @constraint(m, FixLoadShedding_NonDemand[n in LayerNonDemandNodes[CLayer], t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            loadshedding[n,t,s] == 0.0
        )
        # @constraint(m, FixProductionShedding_PV[n in LayerPVNodes[CLayer], t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
        #     productionshedding[n,t,s] == PVNodeDict[CLayer][n] * NormPV[t][MPCScenarios[t-TimeChoice+1,s]]
        # )
        if CLayer > 1
            @constraint(m, FixPin2zero[t = TimeChoice:TimeEnd], p_in[t] == 0.0)
        end
    end

    # Balance
    if CLayer > 1
        # with battery (and PV)
        @constraint(m, Balance_Battery[n in LayerBatteryNodes[CLayer],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (loadshedding[n,t,s] - productionshedding[n,t,s] - pflow[n,t,s]
             + sum(pflow[i,t,s] for i in LayerChildrenNodes[CLayer][n])
             + batterydischarge[n,t,s] - batterycharge[n,t,s]
             == balance_RHS[n,t,s])
        )
        # without battery
        @constraint(m, Balance_NonBattery[n in LayerNonBatteryNodes[CLayer],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (loadshedding[n,t,s] - productionshedding[n,t,s] - pflow[n,t,s]
             + sum(pflow[i,t,s] for i in LayerChildrenNodes[CLayer][n])
             == balance_RHS[n,t,s])
        )
    else
        # root node
        @constraint(m, Balance_Root[t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (sum(pgeneration[g,t,s] for g=1:NGenerators)
            + loadshedding[1,t,s] - productionshedding[1,t,s] - pflow[1,t,s]
             + sum(pflow[i,t,s] for i in LayerChildrenNodes[CLayer][1])
             == 0.0)
        )
        # with demand and battery (and PV)
        @constraint(m, Balance_Battery[n in setdiff(LayerBatteryNodes[CLayer],1),t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (loadshedding[n,t,s] - productionshedding[n,t,s] - pflow[n,t,s]
             + sum(pflow[i,t,s] for i in LayerChildrenNodes[CLayer][n])
             + batterydischarge[n,t,s] - batterycharge[n,t,s]
             == balance_RHS[n,t,s])
        )
        # without demand
        @constraint(m, Balance_NonBattery[n in setdiff(LayerNonBatteryNodes[CLayer],1), t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (loadshedding[n,t,s] - productionshedding[n,t,s] - pflow[n,t,s]
             + sum(pflow[i,t,s] for i in LayerChildrenNodes[CLayer][n])
             == balance_RHS[n,t,s])
        )
    end
    # balance_RHS assignment
    @constraint(m, balance_RHS_NonDemand[n in LayerNonDemandNodes[CLayer], t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
        balance_RHS[n,t,s] == 0.0
    )
    @constraint(m, balance_RHS_NonBatteryDemand[n in setdiff(LayerDemandNodes[CLayer],LayerPVNodes[CLayer]), t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
        balance_RHS[n,t,s] == DemandDataDict[CLayer][n][t]
    )
    @constraint(m, balance_RHS_PVDemand[n in LayerPVNodes[CLayer], t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
        balance_RHS[n,t,s] == DemandDataDict[CLayer][n][t] - (PVNodeDict[CLayer][n] * NormPV[t][MPCScenarios[t-TimeChoice+1,s]])
    )
    # Dynamics
    if TimeChoice == 1  # At stage 1
        @constraint(m, BatteryDynamics_Stage1[n in LayerBatteryNodes[CLayer], s=1:MPCScenarioSize],
            (storage[n,1,s] - BatteryChargeEfficiency * batterycharge[n,1,s] / 4 # 15min
            + batterydischarge[n,1,s]/BatteryDischargeEfficiency / 4  # 15min
            == InitialStorage)
        );
        if TimeEnd - TimeChoice > 0
            @constraint(m, BatteryDynamics[n in LayerBatteryNodes[CLayer], t=TimeChoice+1:TimeEnd, s=1:MPCScenarioSize],
                (storage[n,t,s] - BatteryChargeEfficiency * batterycharge[n,t,s] / 4 # 15min
                + batterydischarge[n,t,s]/BatteryDischargeEfficiency / 4 # 15min
                == storage[n,t-1,s])
            );
        end
    elseif TimeChoice < H
        @constraint(m, BatteryDynamics_CurrentStage[n in LayerBatteryNodes[CLayer], s=1:MPCScenarioSize],
            (storage[n,TimeChoice,s] - BatteryChargeEfficiency * batterycharge[n,TimeChoice,s]  / 4 # 15min
            + batterydischarge[n,TimeChoice,s]/BatteryDischargeEfficiency  / 4 # 15min
            == ImplementedSolution[CLayer,CSubnet,TimeChoice-1].storage[n] )
        );
        if TimeEnd - TimeChoice > 0
            @constraint(m, BatteryDynamics[n in LayerBatteryNodes[CLayer], t=TimeChoice+1:TimeEnd, s=1:MPCScenarioSize],
                (storage[n,t,s] - BatteryChargeEfficiency * batterycharge[n,t,s]  / 4 # 15min
                + batterydischarge[n,t,s]/BatteryDischargeEfficiency / 4 # 15min
                == storage[n,t-1,s])
            );
        end
    else # TimeChoice == H
        @constraint(m, BatteryDynamics_CurrentStage[n in LayerBatteryNodes[CLayer], s=1:MPCScenarioSize],
            (storage[n,TimeChoice,s] - BatteryChargeEfficiency * batterycharge[n,TimeChoice,s]  / 4 # 15min
            + batterydischarge[n,TimeChoice,s]/BatteryDischargeEfficiency / 4 # 15min
            == ImplementedSolution[CLayer,CSubnet,TimeChoice-1].storage[n])
        );
    end
    # Flow Limits
    @constraint(m, FlowMax[i = 1:NLayerLines[CLayer],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize], SLimit[CLayer] >= pflow[i,t,s])
    @constraint(m, FlowMin[i = 1:NLayerLines[CLayer],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize], pflow[i,t,s] >= -SLimit[CLayer])
    if CLayer == 1
        # Generation Limits
        @constraint(m, GenerationMax[g = 1:NGenerators,t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],(PGenerationMax[g] >= pgeneration[g,t,s]))
        @constraint(m, GenerationMin[g = 1:NGenerators,t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],(pgeneration[g,t,s] >= PGenerationMin[g]))
        # Production Shedding Limits
        @constraint(m, ProductionSheddingMax[t = TimeChoice:TimeEnd, s=1:MPCScenarioSize], sum(pgeneration[g,t,s] for g=1:NGenerators) >= productionshedding[1,t,s])
    end
    # Load Shedding Limits
    # @constraint(m, LoadSheddingMax[n = 1:NLayerNodes[CLayer],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize], D[CLayer,n,t,MPCScenarios[t-TimeChoice+1,s]] >= loadshedding[n,t,s])
    # Storage Limits
    @constraint(m, StorageMax[n in LayerBatteryNodes[CLayer],t=TimeChoice:TimeEnd, s=1:MPCScenarioSize],storage[n,t,s] <= BatteryCapacity);
    # Charging Limits
    @constraint(m, BatteryChargeMax[n in LayerBatteryNodes[CLayer],t=TimeChoice:TimeEnd, s=1:MPCScenarioSize],batterycharge[n,t,s] <= BatteryChargeRate);
    # Discharging Limits
    @constraint(m, BatteryDischargeMax[n in LayerBatteryNodes[CLayer],t=TimeChoice:TimeEnd, s=1:MPCScenarioSize],batterydischarge[n,t,s] <= BatteryChargeRate);
    # Equality of p_in and pflow & Equality of p_out and pflow
    if CLayer == 1
        @constraint(m, Pout_Flow_equality[i = 1:NSubnetworks[CLayer+1],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (p_out[i,t] + pflow[NLayerNodes[CLayer] + i,t,s] == 0.0)
        )
        @constraint(m, RootFlowZero[t=TimeChoice:TimeEnd, s=1:MPCScenarioSize],pflow[1,t,s]==0);
    elseif CLayer == NLayers
        @constraint(m, Pin_Flow_equality[t = TimeChoice:TimeEnd, s=1:MPCScenarioSize], p_in[t] + pflow[1,t,s] == 0.0)
    else
        @constraint(m, Pin_Flow_equality[t = TimeChoice:TimeEnd, s=1:MPCScenarioSize], p_in[t] + pflow[1,t,s] == 0.0)
        @constraint(m, Pout_Flow_equality[i = 1:NSubnetworks[CLayer+1],t = TimeChoice:TimeEnd, s=1:MPCScenarioSize],
            (p_out[i,t] + pflow[NLayerNodes[CLayer] + i,t,s] == 0.0)
        )
    end
    # Equality constraints for the variable
    if MPCScenarioSize > 1
        if CLayer == 1
            @constraint(m, Equal_pgeneration[g=1:NGenerators, s=1:MPCScenarioSize],pgeneration[g,TimeChoice,s] == sum(pgeneration[g,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
        end
        @constraint(m, Equal_loadshedding[n=1:NLayerNodes[CLayer], s=1:MPCScenarioSize],loadshedding[n,TimeChoice,s] == sum(loadshedding[n,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
        # @constraint(m, Equal_productionshedding[n=1:NLayerNodes[CLayer], s=1:MPCScenarioSize],productionshedding[n,TimeChoice,s] == sum(productionshedding[n,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
        @constraint(m, Equal_storage[n in LayerBatteryNodes[CLayer], s=1:MPCScenarioSize],storage[n,TimeChoice,s] == sum(storage[n,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
        @constraint(m, Equal_batterycharge[n in LayerBatteryNodes[CLayer], s=1:MPCScenarioSize],batterycharge[n,TimeChoice,s] == sum(batterycharge[n,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
        @constraint(m, Equal_batterydischarge[n in LayerBatteryNodes[CLayer], s=1:MPCScenarioSize],batterydischarge[n,TimeChoice,s] == sum(batterydischarge[n,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
        @constraint(m, Equal_pflow[i=1:NLayerLines[CLayer],s=1:MPCScenarioSize],pflow[i,TimeChoice,s] == sum(pflow[i,TimeChoice,k] for k=1:MPCScenarioSize)/MPCScenarioSize)
    end
    # Coupling Constraint: equality of p_out and p_in
    if CLayer < NLayers
        @constraint(m, Coupling[i = 1:NSubnetworks[CLayer+1], t = TimeChoice:TimeEnd],
            (p_out[i,t] + (-1)*(sum(lambda[j] * TemporalSolutionSet[CLayer+1,i,j].p_in[t] for j = 1:Iteration))
             == 0.0)
        )
        @constraint(m, LambdaSum, sum(lambda) == 1)
    end

    # At the last iteration
    # Set p_in and lambda to be the suggestion Last_p_in and Last_lambda
    if LastIteration && CLayer > 1
        @constraint(m, FixInterfaceFlow[t = TimeChoice:TimeEnd], p_in[t] == Last_p_in[t-TimeChoice+1])
        if Last_lambda != Any
            @constraint(m, FixLambda, lambda .== Last_lambda )
        end
    end
    # if TimeChoice > 1 #&& LastIteration
    #     # print(m);
    # end

    # If you want to store the problem struct
    # if TimeChoice == 1
    #     if CLayer == 1
    #     elseif CLayer == NLayers
    #     else
    #     end
    # elseif TimeChoice == H
    #     if CLayer == 1
    #     elseif CLayer == NLayers
    #     else
    #     end
    # else # 1<TimeChoice<H
    #     if CLayer == 1
    #     elseif CLayer == NLayers
    #     else
    #     end
    # end

    # Solve the problem
    optimize!(m);

    # store solution
    sol = MPCSolutionStruct()
    # compute the max prod shedding to extract from the objective value
    max_productionshedding = maximum([sum(value.(productionshedding[n,t,s]) for n in 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd)  for  s in 1:MPCScenarioSize])
    sol.ObjectiveValue = objective_value(m) - VOLP * max_productionshedding
    sol.SubnetworkCost = 1/4 * sum(VOLL * value.(loadshedding[n,TimeChoice,1]) for n = 1:NLayerNodes[CLayer])
    if CLayer == 1
        sol.SubnetworkCost += 1/4 * sum(MargCost[g] * value.(pgeneration[g,TimeChoice,1]) for g = 1:NGenerators)
    end
    sol.CostForUpperLayer = value.(costupperbound)  - VOLP * max_productionshedding
    if CLayer < NLayers
        sol.CostForUpperLayer += sum(value.(lambda)[j] * TemporalSolutionSet[CLayer+1,sn,j].CostForUpperLayer for sn in ChildrenNetworks, j in 1:Iteration)
    end
    # println("Objective: ",objective_value(m));
    # println("SubnetworkCost: ",sol.SubnetworkCost)

    # vars
    sol.pflow = value.(pflow)[:,TimeChoice,1]
    sol.loadshedding = value.(loadshedding)[:,TimeChoice,1]
    sol.productionshedding = value.(productionshedding)[:,TimeChoice,1]
    sol.storage = value.(storage)[:,TimeChoice,1]
    sol.batterycharge = value.(batterycharge)[:,TimeChoice,1]
    sol.batterydischarge = value.(batterydischarge)[:,TimeChoice,1]
    sol.costupperbound = value.(costupperbound)
    if CLayer == 1
        sol.pgeneration = value.(pgeneration)[:,TimeChoice,1]
        sol.p_out = value.(p_out)
        sol.lambda = value.(lambda)
        # dual
        sol.Coupling = dual.(Coupling)
        sol.LambdaSum = dual.(LambdaSum)
    elseif CLayer == NLayers
        sol.p_in = value.(p_in)
    else
        sol.p_in = value.(p_in)
        sol.p_out = value.(p_out)
        sol.lambda = value.(lambda)
        # dual
        sol.Coupling = dual.(Coupling)
        sol.LambdaSum = dual.(LambdaSum)
    end
    # return the solution
    return sol
end

## Implementation function
function ImplementMPC(MPCType, RealScenario, FutureStepSize = 0, KnowCurrentOutcome = false, sampleID = 1)
    "
    Implement a MPC of MPCType on a RealScenario
    "
    ReducedCost = [[] for t in T]
    UpperBound = [[] for t in T]
    LowerBound = [[] for t in T]
    TotalCost = zeros(H);
    TotalStorage = zeros(H);
    TotalPS = zeros(H);
    TotalLS = zeros(H);
    TotalProduction = zeros(H);
    OnlineTime = zeros(H);
    TimeStages = T;
    if sampleID > 1
        TimeStages =findall(NOutcomes.!=1)[1]:H;
    end
    for t in TimeStages
        # Phase 1: implement nested DW with MPC scenarios
        println("Sample $sampleID Stage $t : Phase 1 - Implement nested DW")
        if !KnowCurrentOutcome
            println("Stage $t : unknown current outcome ")
        else
            println("Stage $t : outcome ", RealScenario[t])
        end
        s_time = time()
        # scenario generation
        MPCScenarios = MPCScenarioGeneration(t,FutureStepSize,KnowCurrentOutcome,RealScenario,NScenarios);
        println("Generated mpc scenarios: ", MPCScenarios)
        # start algorithm
        iter = 0
        for k = 1:K+1
            if k == 1
                println("Stage $t : initial iteration")
                # set the initial solution of Subproblems
                # solve the subproblems from the bottom
                for l = NLayers:-1:1
                    for sn = 1:NSubnetworks[l]
                        # println("Stage $t, PROBLEM ($l,$sn) Outcome ",RealScenario[l,sn,t]," : initial iteration")
                        # println(" MPC Scenario: ",MPCScenarios[l,sn,:,:]')
                        TemporalSolutionSet[l,sn,k] = MPCType([l,sn],MPCScenarios,t,k)
                        # println("=======================================\n")
                    end
                end
            else
                iter += 1;
                println("Stage $t : $iter-th iteration")
                # solve from the bottom layer
                for l = NLayers:-1:2
                    for sn = 1:NSubnetworks[l]
                        up_n = ceil(Int64, sn/NSubnetworks[l]) # node id of the ancestor
                        up_i = Int64((sn-1) % NSubnetworks[l] + 1) # branch id of the ancesotr which connects to the current node n
                        Ancestor = (l-1,up_n)
                        interface_price = TemporalSolutionSet[l-1,up_n,k-1].Coupling[up_i,:]
                        # println("Stage $t, PROBLEM ($l,$sn) Outcome ",RealScenario[t]," : iteration $iter, price: ", collect(interface_price),", Ancestor $Ancestor-$up_i")
                        # println(" MPC Scenario: ",MPCScenarios[l,sn,:,:]')
                        TemporalSolutionSet[l,sn,k] = MPCType([l,sn],MPCScenarios,t,k, interface_price)
                        # println("=======================================\n")
                    end
                end
                # for the root node (master)
                # println("Stage $t, PROBLEM (1,1) Outcome ",RealScenario[1,1,t]," : iteration $iter")
                # println(" MPC Scenario: ",MPCScenarios[1,1,:,:]')
                TemporalSolutionSet[1,1,k] = MPCType([1,1],MPCScenarios,t,k);
                # println("=======================================\n")
                # compute the reduced cost for the master
                rc = sum(TemporalSolutionSet[2,i,k].ObjectiveValue for i=1:NSubnetworks[2]) - TemporalSolutionSet[1,1,k].LambdaSum
                ReducedCost[t] = push!(ReducedCost[t],rc)
                # upper and lower bound
                UpperBound[t] = push!(UpperBound[t], TemporalSolutionSet[1,1,k].ObjectiveValue)
                LowerBound[t] = push!(LowerBound[t], (TemporalSolutionSet[1,1,k].ObjectiveValue
                            + minimum([TemporalSolutionSet[2,m,k].ObjectiveValue for m = 1:NSubnetworks[2]])
                            - TemporalSolutionSet[1,1,k].LambdaSum))
                # check the reduced cost:
                # println("Upper Bound: ",round(UpperBound[t][iter],digits=2) , " Lower Bound: ",round(LowerBound[t][iter],digits=2) , " Reduced Cost: ",round(ReducedCost[t][iter],digits=2))
                if ReducedCost[t][iter] >= -1e-10 && iter > 1
                    break
                end
            end
        end
        println("Stage $t : Algorithm finished!!!!!")
        # println("   Number of Iteratoin: $iter")


        ## Phase 2: excecute the decision.
        # we need to resolve the problem with fixed interface flow (p_in/p_out) and fixed lambda
        println("-------------------------------")
        println("Sample $sampleID Stage $t : Phase 2 - excecute the decision")
        if !KnowCurrentOutcome
            # we observe the current outcome NOW
            println("Stage $t : current outcome observed: ", RealScenario[t])
        end
        println("Stage $t : start solve with fixed interface flow and lambda")
        LastIteration = iter + 1
        bestlambda = zeros(LastIteration,NLayers-1,NSubnetworks[NLayers])
        final_p_in = zeros(NLayers,NSubnetworks[NLayers], size(MPCScenarios,1))
        bestlambda[:,1,1] = TemporalSolutionSet[1,1,LastIteration].lambda

        # Layer 1
        # println("\nObjecive of Master (Whole Network): ",TemporalSolutionSet[1,1,LastIteration].ObjectiveValue)
        # println("   Cost at (1,1): ",round(TemporalSolutionSet[1,1,LastIteration].SubnetworkCost,digits=5))
        ImplementedSolution[1,1,t] = TemporalSolutionSet[1,1,LastIteration]
        TotalCost[t] = ImplementedSolution[1,1,t].SubnetworkCost
        TotalStorage[t] = sum(ImplementedSolution[1,1,t].storage)
        TotalPS[t] = sum(ImplementedSolution[1,1,t].productionshedding)
        TotalLS[t] = sum(ImplementedSolution[1,1,t].loadshedding)
        TotalProduction[t] = sum(ImplementedSolution[1,1,t].pgeneration)
        # Layer 2 to NLayers
        for l = 2:NLayers
            for sn = 1:NSubnetworks[l]
                up_n = ceil(Int64, sn/NSubnetworks[l]) # node id of the ancestor
                up_i = Int64((sn-1) % NSubnetworks[l] + 1) # branch id of the ancesotr which connects to the current node n
                # interface_price = TemporalSolutionSet[l-1,up_n,k-1].Coupling[up_i,:]
                # up_n = ceil(Int64, sn/NBranches) # node id of the ancestor
                # up_i = Int64((sn-1) % NBranches + 1) # branch id of the ancesotr which connects to the current node n
                interface_price = TemporalSolutionSet[l-1,up_n,LastIteration].Coupling[up_i,:]
                # [final_p_in[l,sn,u-t+1] = sum(bestlambda[k,l-1,up_n] .* TemporalSolutionSet[l,sn,k].p_in[u] for k=1:LastIteration) for u=t:T]
                ## I did want to compute the p_in from the best lambda and temporal p_in but
                ## it didnt work so changed it to the corresponding p_out
                final_p_in[l,sn,:] = ImplementedSolution[l-1,up_n,t].p_out[up_i,:]
                # println("---final p_in at ($l,$sn): ",final_p_in[l,sn,:], " with price: ",interface_price)
                if l < NLayers
                   for k=1:LastIteration
                       bestlambda[k,l,sn] = sum(TemporalSolutionSet[l,sn,q].lambda[k] .* bestlambda[q,l-1,up_n] for q = k:LastIteration)
                   end
                   ImplementedSolution[l,sn,t] = MPCType([l,sn],RealScenario[t:t+size(MPCScenarios,1)-1],t,LastIteration,interface_price,
                                        true, final_p_in[l,sn,:], bestlambda[:,l,sn])
                   # ImplementedSolution[l,sn,t] = MPCType([l,sn],MPCScenarios,t,LastIteration,interface_price,
                   #                      true, final_p_in[l,sn,:], bestlambda[:,l,sn])
                else
                    ImplementedSolution[l,sn,t] = MPCType([l,sn],RealScenario[t:t+size(MPCScenarios,1)-1],t,LastIteration,interface_price,true,final_p_in[l,sn,:])
                   # ImplementedSolution[l,sn,t] = MPCType([l,sn],MPCScenarios,t,LastIteration,interface_price,true,final_p_in[l,sn,:])
                end
                # println("    Cost at ($l,$sn) : ", round(ImplementedSolution[l,sn,t].SubnetworkCost,digits=5))
                # println("    Storage at ($l,$sn) : ", round(sum(ImplementedSolution[l,sn,t].storage)))
                TotalStorage[t] += sum(ImplementedSolution[l,sn,t].storage)
                TotalCost[t] += ImplementedSolution[l,sn,t].SubnetworkCost
                TotalPS[t] += sum(ImplementedSolution[l,sn,t].productionshedding)
                TotalLS[t] += sum(ImplementedSolution[l,sn,t].loadshedding)
           end
        end
        OnlineTime[t] = time() - s_time
        println("Total cost over the network (\$): ",round(TotalCost[t],digits=5))
        println("Total production over the network (kW): ",round(TotalProduction[t],digits=5))
        println("Total storage over the network (kWh): ",round(TotalStorage[t],digits=5))
        println("Total load shedding over the network (kW): ",round(TotalLS[t],digits=5))
        println("Total production shedding over the network (kW): ",round(TotalPS[t],digits=5))
        println("Computation time (s): ",round(OnlineTime[t],digits=2))
        println("===============================")
    end
    # println("Total cost: ",round(sum(TotalCost)),", at each stage: ", TotalCost)
    # println("===============================\n\n")
    return MPCResultStruct(TotalCost,TotalProduction,TotalStorage,TotalLS,OnlineTime,ImplementedSolution) #TotalCost, TotalProduction, TotalStorage, OnlineTime
end

# function SAMPCSubproblem(NetworkID, MPCScenarios, TimeChoice, Iteration, dual_pi = zeros(T), LastIteration=false, Last_p_in = Any, Last_lambda = Any)
#     """
#     solve a subproblem at (l,sn) with inputs of the upper layer
#     args    - NetworkID: ID of Subproblem [layerID, subnetID]
#             - MPCScenarios[t,s]: outcome at stage t of scenario s
#             - TimeChoice
#             - Iteration: iteration upto now
#             - dual_pi: dual variable of coupling constraint from the upper layer master
#             - LastIteration: wheather the iteration is the last one or not
#             - Last_p_in: if we are at the last iteration, Last_p_in will be the suggestion from the upper layer
#             - Last_lambda: Last_lambda will be the suggestion from the upper layer
#     return  - sol: the struct Solution_Multi which has the solution
#     """
#     CLayer = NetworkID[1]; # current layer
#     CSubnet = NetworkID[2]; # current subnetwork
#
#     ## model
#     m = Model(with_optimizer(Gurobi.Optimizer, LogToConsole=0));
#
#     ## Variables
#     @variable(m, loadshedding[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:NScenarios] >= 0)
#     @variable(m, productionshedding[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:NScenarios] >= 0)
#     @variable(m, storage[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:NScenarios] >= 0)
#     @variable(m, batterycharge[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:NScenarios] >= 0)
#     @variable(m, batterydischarge[1:NLayerNodes[CLayer], TimeChoice:TimeEnd, 1:NScenarios] >= 0)
#     @variable(m, pflow[1:NLayerLines[CLayer], TimeChoice:TimeEnd, 1:NScenarios])
#     if CLayer == 1
#         @variable(m, pgeneration[1:NGenerators, TimeChoice:TimeEnd, 1:NScenarios])
#         @variable(m, p_out[1:NBranches, TimeChoice:TimeEnd])  # independent on scenarios
#         @variable(m, lambda[1:Iteration] >=0)
#     elseif CLayer == NLayers
#         @variable(m, p_in[TimeChoice:TimeEnd])  # independent on scenarios
#     else
#         @variable(m, p_out[1:NBranches, TimeChoice:TimeEnd]) # independent on scenarios
#         @variable(m, p_in[TimeChoice:TimeEnd]) # independent on scenarios
#         @variable(m, lambda[1:Iteration] >=0)
#     end
#
#     ## objective
#     if CLayer < NLayers
#         ChildrenNetworks = [NBranches*(CSubnet-1)+j for j = 1:NBranches] # set of children networks
#         Child_s = ChildrenNetworks[1] - 1 # for enumerate child network ID
#     end
#     if CLayer == 1 # root layer: not the reduced cost
#         @objective(m, Min,
#             (1/NScenarios * sum(sum( MargCost[g,t] * pgeneration[g,t,s] for g=1:NGenerators)
#              + sum(VOLL[CLayer,n,t] * loadshedding[n,t,s] for n = 1:NLayerNodes[CLayer]) for t = TimeChoice:TimeEnd, s = 1:NScenarios)
#              + sum(lambda[j] * sum(TemporalSolutionSet[CLayer+1,sn,j].CostForUpperLayer for sn in ChildrenNetworks) for j = 1:Iteration))
#         );
#     elseif CLayer == NLayers # bottom layer: no lower layers
#         @objective(m, Min,
#             (1/NScenarios * sum(#= MargCost[CLayer,n,t] * pgeneration[n,t,s] +=# VOLL[CLayer,n,t] * loadshedding[n,t,s] for n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd, s = 1:NScenarios)
#              - sum((dual_pi[t]) * (-1) * p_in[t] for t = TimeChoice:TimeEnd))
#         );
#     else # inner layers: with reduced cost & lower layers
#         @objective(m, Min,
#             (1/NScenarios * sum(#= MargCost[CLayer,n,t] * pgeneration[n,t,s] +=# VOLL[CLayer,n,t] * loadshedding[n,t,s] for n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd, s = 1:NScenarios)
#              - sum((dual_pi[t]) * (-1) * p_in[t] for t = TimeChoice:TimeEnd) +
#             sum(lambda[j] * sum(TemporalSolutionSet[CLayer+1,sn,j].CostForUpperLayer for sn in ChildrenNetworks) for j = 1:Iteration))
#         );
#     end
#
#     ## constraints
#     # Set loadshedding to be the demend (because it is always feasible)
#     if Iteration == 1
#         @constraint(m, FixLoadShedding[n = 1:NLayerNodes[CLayer], t = TimeChoice:TimeEnd, s=1:NScenarios], loadshedding[n,t,s] == D[CLayer,n,t,MPCScenarios[t-TimeChoice+1,s]])
#         if CLayer > 1
#             @constraint(m, FixPin2zero[t = TimeChoice:TimeEnd], p_in[t] == 0.0)
#         end
#     end
#
#     # Balance
#     if CLayer > 1
#         @constraint(m, Balance[n = 1:NLayerNodes[CLayer],t = TimeChoice:TimeEnd, s=1:NScenarios],
#             (#=pgeneration[n,t,s] + =#loadshedding[n,t,s] - productionshedding[n,t,s] - pflow[n,t,s]
#              + sum(pflow[i,t,s] for i in ChildrenLines[n])
#              + batterydischarge[n,t,s] - batterycharge[n,t,s]
#              == D[CLayer,n,t,MPCScenarios[t-TimeChoice+1,s]])
#         )
#     else
#         @constraint(m, Balance_rootnode[t = TimeChoice:TimeEnd, s=1:NScenarios],
#             (sum(pgeneration[g,t,s] for g=1:NGenerators) + loadshedding[1,t,s] - productionshedding[1,t,s] - pflow[1,t,s]
#              + sum(pflow[i,t,s] for i in ChildrenLines[1])
#              + batterydischarge[1,t,s] - batterycharge[1,t,s]
#              == D[CLayer,1,t,MPCScenarios[t-TimeChoice+1,s]])
#         )
#         @constraint(m, Balance[n = 2:NNodes,t = TimeChoice:TimeEnd, s=1:NScenarios],
#             (#=pgeneration[n,t,s] + =#loadshedding[n,t,s] - productionshedding[n,t,s] - pflow[n,t,s]
#              + sum(pflow[i,t,s] for i in ChildrenLines[n])
#              + batterydischarge[n,t,s] - batterycharge[n,t,s]
#              == D[CLayer,n,t,MPCScenarios[t-TimeChoice+1,s]])
#         )
#     end
#     # Dynamics
#     if TimeChoice == 1
#         @constraint(m, BatteryDynamics_Stage1[n = 1:NLayerNodes[CLayer], s=1:NScenarios],
#             (storage[n,1,s] - BatteryChargeEfficiency[CLayer,n] * batterycharge[n,1,s]
#             + batterydischarge[n,1,s]/BatteryDischargeEfficiency[CLayer,n]
#             == InitialStorage[CLayer,n])
#         );
#         @constraint(m, BatteryDynamics[n = 1:NLayerNodes[CLayer], t=2:T, s=1:NScenarios],
#             (storage[n,t,s] - BatteryChargeEfficiency[CLayer,n] * batterycharge[n,t,s]
#             + batterydischarge[n,t,s]/BatteryDischargeEfficiency[CLayer,n]
#             == storage[n,t-1,s])
#         );
#     elseif TimeChoice < T
#         @constraint(m, BatteryDynamics_CurrentStage[n = 1:NLayerNodes[CLayer], s=1:NScenarios],
#             (storage[n,TimeChoice,s] - BatteryChargeEfficiency[CLayer,n] * batterycharge[n,TimeChoice,s]
#             + batterydischarge[n,TimeChoice,s]/BatteryDischargeEfficiency[CLayer,n]
#             == ImplementedSolution[CLayer,CSubnet,TimeChoice-1].storage[n]#=[n,TimeChoice-1]=#)
#         );
#         @constraint(m, BatteryDynamics[n = 1:NLayerNodes[CLayer], t=TimeChoice+1:T, s=1:NScenarios],
#             (storage[n,t,s] - BatteryChargeEfficiency[CLayer,n] * batterycharge[n,t,s]
#             + batterydischarge[n,t,s]/BatteryDischargeEfficiency[CLayer,n]
#             == storage[n,t-1,s])
#         );
#     else # TimeChoice == T
#         @constraint(m, BatteryDynamics_CurrentStage[n = 1:NLayerNodes[CLayer], s=1:NScenarios],
#             (storage[n,TimeChoice,s] - BatteryChargeEfficiency[CLayer,n] * batterycharge[n,TimeChoice,s]
#             + batterydischarge[n,TimeChoice,s]/BatteryDischargeEfficiency[CLayer,n]
#             == ImplementedSolution[CLayer,CSubnet,TimeChoice-1].storage[n]#=[n,TimeChoice-1]=#)
#         );
#     end
#     # Flow Limits
#     @constraint(m, FlowMax[i = 1:NLayerLines[CLayer],t = TimeChoice:TimeEnd, s=1:NScenarios], LineCap[CLayer,i,t] >= pflow[i,t,s])
#     @constraint(m, FlowMin[i = 1:NLayerLines[CLayer],t = TimeChoice:TimeEnd, s=1:NScenarios], pflow[i,t,s] >= -LineCap[CLayer,i,t])
#     if CLayer == 1
#         # Generation Limits
#         @constraint(m, GenerationMax[g = 1:NGenerators,t = TimeChoice:TimeEnd, s=1:NScenarios],(PMax[g,t] >= pgeneration[g,t,s]))
#         @constraint(m, GenerationMin[g = 1:NGenerators,t = TimeChoice:TimeEnd, s=1:NScenarios],(pgeneration[g,t,s] >= PMin[g,t]))
#         # Production Shedding Limits
#         @constraint(m, ProductionSheddingMax[t = TimeChoice:TimeEnd, s=1:NScenarios], sum(pgeneration[g,t,s] for g=1:NGenerators) >= productionshedding[1,t,s])
#     end
#     # Load Shedding Limits
#     @constraint(m, LoadSheddingMax[n = 1:NLayerNodes[CLayer],t = TimeChoice:TimeEnd, s=1:NScenarios], D[CLayer,n,t,MPCScenarios[t-TimeChoice+1,s]] >= loadshedding[n,t,s])
#     # Storage Limits
#     @constraint(m, StorageMax[n = 1:NLayerNodes[CLayer],t=TimeChoice:TimeEnd, s=1:NScenarios],storage[n,t,s] <= BatteryCapacity[CLayer,n]);
#     # Charging Limits
#     @constraint(m, BatteryChargeMax[n = 1:NLayerNodes[CLayer],t=TimeChoice:TimeEnd, s=1:NScenarios],batterycharge[n,t,s] <= BatteryChargeRate[CLayer,n]);
#     # Discharging Limits
#     @constraint(m, BatteryDischargeMax[n = 1:NLayerNodes[CLayer],t=TimeChoice:TimeEnd, s=1:NScenarios],batterydischarge[n,t,s] <= BatteryChargeRate[CLayer,n]);
#     # Equality of p_in and pflow & Equality of p_out and pflow
#     if CLayer == 1
#         @constraint(m, Pout_Flow_equality[i = 1:NBranches,t = TimeChoice:TimeEnd, s=1:NScenarios],
#             (p_out[i,t] + pflow[NNodes + i,t,s] == 0.0)
#         )
#     elseif CLayer == NLayers
#         @constraint(m, Pin_Flow_equality[t = TimeChoice:TimeEnd, s=1:NScenarios], p_in[t] + pflow[1,t,s] == 0.0)
#     else
#         @constraint(m, Pin_Flow_equality[t = TimeChoice:TimeEnd, s=1:NScenarios], p_in[t] + pflow[1,t,s] == 0.0)
#         @constraint(m, Pout_Flow_equality[i = 1:NBranches,t = TimeChoice:TimeEnd, s=1:NScenarios],
#             (p_out[i,t] + pflow[NNodes + i,t,s] == 0.0)
#         )
#     end
#     # Equality constraints for the variable
#     if CLayer == 1
#         @constraint(m, Equal_pgeneration[g=1:NGenerators, s=1:NScenarios],pgeneration[g,TimeChoice,s] == sum(pgeneration[g,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#     end
#     @constraint(m, Equal_loadshedding[n=1:NLayerNodes[CLayer], s=1:NScenarios],loadshedding[n,TimeChoice,s] == sum(loadshedding[n,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#     @constraint(m, Equal_productionshedding[n=1:NLayerNodes[CLayer], s=1:NScenarios],productionshedding[n,TimeChoice,s] == sum(productionshedding[n,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#     @constraint(m, Equal_storage[n=1:NLayerNodes[CLayer], s=1:NScenarios],storage[n,TimeChoice,s] == sum(storage[n,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#     @constraint(m, Equal_batterycharge[n=1:NLayerNodes[CLayer], s=1:NScenarios],batterycharge[n,TimeChoice,s] == sum(batterycharge[n,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#     @constraint(m, Equal_batterydischarge[n=1:NLayerNodes[CLayer], s=1:NScenarios],batterydischarge[n,TimeChoice,s] == sum(batterydischarge[n,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#     @constraint(m, Equal_pflow[i=1:NLayerLines[CLayer],s=1:NScenarios],pflow[i,TimeChoice,s] == sum(pflow[i,TimeChoice,k] for k=1:NScenarios)/NScenarios)
#
#     # Coupling Constraint: equality of p_out and p_in
#     if CLayer < NLayers
#         @constraint(m, Coupling[i = 1:NBranches, t = TimeChoice:TimeEnd],
#             (p_out[i,t] + (-1)*(sum(lambda[j] * TemporalSolutionSet[CLayer+1,i,j].p_in[t] for j = 1:Iteration))
#              == 0.0)
#         )
#         @constraint(m, LambdaSum, sum(lambda) == 1)
#     end
#
#     # At the last iteration
#     # Set p_in and lambda to be the suggestion Last_p_in and Last_lambda
#     if LastIteration && CLayer > 1
#         @constraint(m, FixInterfaceFlow[t = TimeChoice:TimeEnd], p_in[t] == Last_p_in[t-TimeChoice+1])
#         if Last_lambda != Any
#             @constraint(m, FixLambda, lambda .== Last_lambda )
#         end
#     end
#     ## Solve the problem
#     optimize!(m);
#
#     # store solution
#     sol = MPCSolutionStruct()
#     sol.ObjectiveValue = objective_value(m)
#     sol.SubnetworkCost = sum(VOLL[CLayer,n,TimeChoice] * value.(loadshedding[n,TimeChoice,1]) for n = 1:NLayerNodes[CLayer])
#     if CLayer == 1
#         sol.SubnetworkCost += sum( MargCost[g,TimeChoice] * value.(pgeneration[g,TimeChoice,1]) for g=1:NGenerators)
#     end
#     sol.CostForUpperLayer = objective_value(m)
#     if CLayer > 1
#         sol.CostForUpperLayer += sum((dual_pi[t]) * (-1) * value.(p_in)[t] for t = TimeChoice:TimeEnd)
#     end
#     # println("Objective: ",objective_value(m));
#     # println("SubnetworkCost: ",sol.SubnetworkCost)
#     # vars
#     sol.pflow = value.(pflow)[:,TimeChoice,1]
#     sol.loadshedding = value.(loadshedding)[:,TimeChoice,1]
#     sol.productionshedding = value.(productionshedding)[:,TimeChoice,1]
#     sol.storage = value.(storage)[:,TimeChoice,1]
#     sol.batterycharge = value.(batterycharge)[:,TimeChoice,1]
#     sol.batterydischarge = value.(batterydischarge)[:,TimeChoice,1]
#     sol.costupperbound = "no"
#     if CLayer == 1
#         sol.pgeneration = value.(pgeneration)[:,TimeChoice,1]
#         sol.p_out = value.(p_out)
#         sol.lambda = value.(lambda)
#         # dual
#         sol.Coupling = dual.(Coupling)
#         sol.LambdaSum = dual.(LambdaSum)
#     elseif CLayer == NLayers
#         sol.p_in = value.(p_in)
#     else
#         sol.p_in = value.(p_in)
#         sol.p_out = value.(p_out)
#         sol.lambda = value.(lambda)
#         # dual
#         sol.Coupling = dual.(Coupling)
#         sol.LambdaSum = dual.(LambdaSum)
#     end
#     # return the solution
#     return sol
# end
