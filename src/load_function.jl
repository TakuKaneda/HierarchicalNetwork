using Pkg, JuMP, StatsBase, DataFrames, CSV, DelimitedFiles, JSON, FileIO, JLD2
# using LightGraphs, GraphPlot  # for plotting graph

function ReadSubnetworkData(path_to_nodes, path_to_edges, has_lower_layer=false)
    "
    Read data of a subnetwork from Nodes.txt and Edegs.txt
        return  ChildrenNodes:  -Array{Array{Int64,1},1}
                                ChildrenNodes[n] has the node ID of its children
                                e.g. ChildrenNodes[1] = [2, 3, 4]
                                If the subnetwork has lower layers,
                                the corresponding interface nodes also have
                                line ID beggining with NNodes+1,...
                LoadNodeDict:   -Dict{Int64,Float64}
                                Dictionary of node ID to the nominal load
                                Possibly empty for medium voltage network
                g:              -SimpleDiGraph{Int64}
                                Graph struct of LightGraphs which contains the
                                info of the subnetwork (not include the line for
                                lower layers

    "
    nodes_raw = readdlm(path_to_nodes,',', Float64)
    edges_raw = readdlm(path_to_edges,',', Float64)

    ## Create children dictionary with orignal id
    NNodes = size(nodes_raw,1);
    child_dict =  Dict()#Dict{Int64,Int64}()
    for e = 1:size(edges_raw)[1]
        c_node = edges_raw[e,1];
        p_node = edges_raw[e,2];
        if haskey(child_dict,p_node)
            child_dict[p_node] = push!(child_dict[p_node],c_node)
        else
            child_dict[p_node] = [c_node]
        end
    end
    ## Create dictionary of the orignal id -> new node id (from 1 to NNodes)
    root_id = findall(nodes_raw[:,end] .== 1.0); # find() in v0.6
    node_dict = Dict();
    n_id = root_id;
    nodes_queue = Any[nodes_raw[n_id,1][1]]
    n_count = 1
    while n_count <= NNodes
        if !haskey(node_dict,nodes_queue[1])
            node_dict[nodes_queue[1]] = n_count; # assign new node id
            n_count+=1
        end
        if haskey(child_dict, nodes_queue[1])
            c_set = child_dict[nodes_queue[1]] # push set of child nodes to the queue
            for c = 1:length(c_set)
                push!(nodes_queue,c_set[c]);
            end
        end
        deleteat!(nodes_queue,1)
    end
    ## Create dict for nodes with load
    LoadNodeDict = Dict{Int64,Float64}();
    for i in findall(nodes_raw[:,4].==1.0)
        LoadNodeDict[node_dict[nodes_raw[i,1]]] = nodes_raw[i,5]
    end
    ## Creating children nodes with new index (1-NNodes)
    ChildrenNodes = [Int64[] for i = 1:length(node_dict)]
    for k in keys(node_dict)
        if haskey(child_dict,k)
            children = child_dict[k]
            for c in children
                ChildrenNodes[node_dict[k]] = push!(ChildrenNodes[node_dict[k]],node_dict[c])
            end
        end
    end
    ## Create digraph
    g = SimpleDiGraph(NNodes)
    for n = 1:NNodes
        children = ChildrenNodes[n]
        for c in children
            add_edge!(g,c,n)
        end
    end
    ## if there are lower layers -> add edges for them
    if has_lower_layer
        # Adding the id of lines for low voltage grids (NNodes+1,...)
        for c in range(2,stop=length(node_dict)) # except the root node
            ChildrenNodes[c] = push!(ChildrenNodes[c], c+length(node_dict)-1)
        end
    end
    return ChildrenNodes, LoadNodeDict, g
end

function ConvertGeneralTransProb2Array(path_to_json)
    "
    convert a SINGLE TransProb data into Array
    TransProb[t] = Array{Float64,2}:    Transition probability matrix from
                                        outcome k of stage t to outcome j of
                                        stage t+1
      e.g. TransProb[2][3,4]: TransProb from outcome 3 of stage 2 to outcome 4 of stage 3
    "
    tpstr = String(read(path_to_json)); # read .json
    tpjson = JSON.parse(tpstr)  ;

    TransProb = Array{Array}(undef,length(tpjson));
    # first stage: vector
    if length(tpjson[string(2)]) == 1
        TransProb[1] = [convert(Float64,tpjson[string(2)])];
    else
        TransProb[1] = convert(Array{Float64,1},tpjson[string(2)])';
    end

    # from stage 2: possibely matrix
    for t=2:length(tpjson)
        if length(tpjson[string(t)]) == 1
            if length(tpjson[string(t+1)]) == 1
                TransProb[t] = [convert(Float64,tpjson[string(t+1)])];
            else
                TransProb[t] = convert(Array{Float64,1},tpjson[string(t+1)])';
            end
        else
            if length(tpjson[string(t+1)][1]) == 1
                for k = 1:length(tpjson[string(t+1)])
                    if k == 1
                        TransProb[t] = [convert(Float64,tpjson[string(t+1)][k])];
                    else
                        TransProb[t] = vcat(TransProb[t], [convert(Float64,tpjson[string(t+1)][k])]);
                    end
                end
            else
                for k=1:length(tpjson[string(t+1)])
                    if k == 1
                        TransProb[t] = convert(Array{Float64,1},tpjson[string(t+1)][k])'
                    else
                        TransProb[t] = vcat(TransProb[t],convert(Array{Float64,1},tpjson[string(t+1)][k])')
                    end
                end
            end
        end
    end
    return TransProb
end

function ConvertGeneralLatticeNodeValue2Array(path_to_json)
    "
    convert a SINGLE Node value (PV output) data into Array
    NodeValue[t] = Array{Float64,1}:  PV output at stage t
      e.g. NodeValue[36][3]: PV output at outcome 3 of stage 36
    "
    nodestr = String(read(path_to_json)); # read .json
    nodejson = JSON.parse(nodestr)  ;

    NodeValue = Array{Array}(undef,length(nodejson));
    for t = 1:length(nodejson)
        if length(nodejson[string(t)]) == 1
            NodeValue[t] = [convert(Float64,nodejson[string(t)])];
        else
            NodeValue[t] = convert(Array{Float64,1},nodejson[string(t)]);
        end
    end
    return NodeValue
end

function FlatMarginalCostCurve(x, capacity, price_cap = 140.0/1000)
    "marginal cost function
    args:   x           power in kWh
            capacity
            price_cap
    "
    if x > capacity
        return(price_cap)
    elseif x == capacity
        return(price_cap)
    elseif 0.99*capacity < x <= 1*capacity
        return(price_cap)
    elseif 0.98*capacity < x <= 0.99*capacity
        return(95.0/1000)
    elseif 0.95*capacity < x <= 0.98*capacity
        return(80.0/1000)
    elseif 0.9*capacity < x <= 0.95*capacity
        return(70.0/1000)
    elseif 0.85*capacity < x <= 0.9*capacity
        return(65.0/1000)
    elseif 0.75*capacity < x <= 0.85*capacity
        return(55.0/1000)
    elseif 0.6*capacity < x <= 0.75*capacity
        return(50.0/1000)
    elseif 0.4*capacity < x <= 0.6*capacity
        return(40.0/1000)
    elseif 0.2*capacity < x <= 0.4*capacity
        return(35.0/1000)
    elseif 0 <= x <= 0.2*capacity
        return(f1(x,capacity,price_cap))
    else
        return(0.0/1000)
    end
end

function f1(x, capacity = 3.0, price_cap = 140.0/1000)
    # step function
    if 0 < x <= 4/30*capacity
        return(10.0/1000)
    elseif 4/30*capacity < x <= 5/30 * capacity
        return(20.0/1000)
    elseif  5/30 * capacity < x <=  6/30 * capacity
        return(25.0/1000)
    else
        return(0.0/1000)
    end
end

function SamplePath(TransProb, NSamples=1)
    "
    Generate sample paths using TransProb.
    return the resulting sample path of outcome at every stage
    sample_path[t,i]: outcome of layer l at stage t of sample i
    NOTE: All layer has the same outcome!
    "
    n_stage = length(TransProb) + 1
    n_outcome =  maximum([size(TransProb[t],1) for t in 1:length(TransProb)]);
    sample_path = ones(Int64,n_stage, NSamples)
    for i = 1:NSamples
        for t=1:n_stage-1
            sample_path[t+1,i] = sample(1:length(TransProb[t][sample_path[t,i],:]),
                                StatsBase.Weights(TransProb[t][sample_path[t,i],:]))
        end
    end
    if NSamples == 1
        return sample_path[:,1]
    else
        return sample_path
    end
end

function MPCScenarioGeneration(t,FutureStepSize,KnowCurrentOutcome,RealScenario,NScenarios)
    "generate scenarios for MPC implementation
    args    - t: current time step
            - FutureStepSize: step size of looking a head
            - KnowCurrentOutcome: true or false, if true the outcome of the
                                current stage will be the real outcome.
            - RealScenario: real outcome path
            - NScenarios: num of scenario to generate
    output  - MPCScenarios: Array{Int64,2}: (FutureStepSize+1,NSamples) outcomes
    "
    if t + FutureStepSize <= H
        MPCScenarios = zeros(Int64, FutureStepSize+1, NScenarios)
    else
        MPCScenarios = zeros(Int64, H - t + 1, NScenarios)
    end
    # if we DO NOT KNOW the current outcome
    if !KnowCurrentOutcome
        if t > 1
            [MPCScenarios[1,i] = sample(1:length(TransProb[t-1][RealScenario[t-1],:]),
                StatsBase.Weights(TransProb[t-1][RealScenario[t-1],:])) for i =1:NScenarios];
        else
            MPCScenarios[1,:] .= 1;
        end
    else # if we KNOW the current outcome
        MPCScenarios[1,:] .= RealScenario[t]
    end

    # future outcomes by sampling with transition probability
    if FutureStepSize > 0
        for i = 1:NScenarios, tt = 2:size(MPCScenarios,1)
            MPCScenarios[tt,i] = sample(1:length(TransProb[t+tt-2][MPCScenarios[tt-1,i],:]),
                    StatsBase.Weights(TransProb[t+tt-2][MPCScenarios[tt-1,i],:]))
        end
    end
    return MPCScenarios
end
