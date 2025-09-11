include("./graph.jl")


struct MultiLevelGraph #<: AbstractGraph
    num_levels::Int
    levels::Vector{String}
    coarse_to_fine_graphs::Vector{Vector{SubGraph}} # Or lists of edges ? level->quotient node->expanded graph
    fine_neighbors::Vector{Vector{Dict{Int,Any}}}   # for each level, for each node, have fine nbrs
    mixed_nbr_weights::Dict{Set{Tuple{Vararg{String}}},Any}
    id_to_partitions::Vector{Vector{Tuple{Vararg{String}}}}
    partition_to_ids::Vector{Dict{Tuple{Vararg{String}},Int}}
    graphs_by_level::Vector{BaseGraph}             # the graphs at each level
end

const Graph = MultiLevelGraph

""""""
function edge_weight(
    graph::MultiLevelGraph,
    edge_key#::Set{Tuple{Vararg{String}}}
)
    n1, n2 = collect(edge_key)
    if length(n1) == length(n2)
        level = length(n1)
        id_1 = graph.partition_to_ids[level][n1]
        id_2 = graph.partition_to_ids[level][n2]
        g = graph.graphs_by_level[level].simple_graph
        edge_weight = g.weights[id_1, id_2]
    else
        edge_weight = graph.mixed_nbr_weights[Set([n1, n2])]
    end
    return edge_weight
end

""""""
function get_fine_neighbors(graph::MultiLevelGraph, node::Tuple, nbr::Tuple)
    node_level = length(node)
    nbr_level = length(nbr)
    @assert nbr_level < graph.num_levels

    levels = graph.levels
    node_id = graph.partition_to_ids[node_level][node]
    partition_to_ids = graph.partition_to_ids

    if node_level > nbr_level
        g = graph.graphs_by_level[node_level]
        all_nbrs_at_node_level = neighbors(g.simple_graph, node_id)
        fine_nbrs = Vector{Int}(undef, 0)
        for node_nbr in all_nbrs_at_node_level
            nbr_atr = g.node_attributes[node_nbr]
            nbr_name = Tuple([nbr_atr[k] for k in levels[1:node_level]])
            if nbr_name[1:nbr_level] == nbr
                id = partition_to_ids[nbr_level+1][nbr_name[1:nbr_level+1]]
                push!(fine_nbrs, id)
            end
        end
        fine_nbrs = collect(Set(fine_nbrs))
    else # node_level <= nbr_level
        nbr_id = partition_to_ids[node_level][nbr[1:node_level]]
        fine_nbrs = graph.fine_neighbors[node_level][node_id][nbr_id]
        for ii = 1:(nbr_level - node_level)
            nbr_id = partition_to_ids[node_level+ii][nbr[1:node_level+ii]]
            fine_nbrs = fine_nbrs[nbr_id]
        end
        fine_nbrs = collect(keys(fine_nbrs))
    end
    return fine_nbrs
end

""""""
function build_quotient_ids!(
    id_to_partitions::Vector{Vector{Tuple{Vararg{String}}}},
    partition_to_ids::Vector{Dict{Tuple{Vararg{String}},Int}},
    base_graph::BaseGraph, 
    levels::Vector{String};
    cur_level::Int=-1
)
    num_levels = length(levels)
    if cur_level == -1
        cur_level = num_levels
    end

    partition_to_ids[cur_level] = Dict{Tuple{Vararg{String}}, Int}()
    cur_quotient_id = 1
    for (id, attribute) in enumerate(base_graph.node_attributes)
        partition_array = Vector{String}(undef, cur_level)
        levels_for_tuple = levels[1:cur_level]
        for (level_index, level_name) in enumerate(levels_for_tuple)
            partition_array[level_index] = attribute[level_name]
        end
        partition = Tuple(partition_array)
        if haskey(partition_to_ids[cur_level], partition)
            continue
        end
        partition_to_ids[cur_level][partition] = cur_quotient_id
        cur_quotient_id += 1
    end

    num_q_nodes = length(partition_to_ids[cur_level])
    id_to_partitions[cur_level] = Vector{Tuple{Vararg{String}}}(undef,
                                                                 num_q_nodes)
    for (q_node, q_id) in partition_to_ids[cur_level]
        id_to_partitions[cur_level][q_id] = q_node
    end

    if cur_level > 1
        build_quotient_ids!(id_to_partitions, partition_to_ids, base_graph,
                            levels, cur_level=cur_level-1)
    end
end

""""""
function build_oriented_coarse_nbrs!(
    coarse_node_attributes::Vector{Dict{String,Any}},
    partition::Vector{String},
    base_graph::BaseGraph, levels::Vector{String},
    coarse_id_to_partition::Vector{Tuple{Vararg{String}}},
    coarse_partition_to_id::Dict{Tuple{Vararg{String}},Int}
)
    oriented_nbrs_col = base_graph.oriented_nbrs_col
    @assert oriented_nbrs_col != nothing
    fine_node_attributes = base_graph.node_attributes
    coarse_nm_to_brd_fine_nd = Dict{Tuple{Vararg{String}}, Int}()
    for (node_id, node_atr) in enumerate(fine_node_attributes) #can be streamlined
        q_node = Tuple([node_atr[k] for k in partition])
        brd_node = false
        if 0 in node_atr[oriented_nbrs_col]
            brd_node = true
        else
            for nbr in neighbors(base_graph.simple_graph, node_id)
                nbr_atr = fine_node_attributes[nbr]
                q_nbr = Tuple([nbr_atr[k] for k in partition])
                if q_nbr != q_node
                    brd_node = true
                    break
                end
            end
        end
        if brd_node
            coarse_nm_to_brd_fine_nd[q_node] = node_id
        end
    end
    oriented_brd_fine_nds = Dict{Tuple{Vararg{String}}, Vector{Int}}()
    return
    for q_node in keys(coarse_nm_to_brd_fine_nd) #[("GASTON",),("MCDOWELL",),("IREDELL",)]
        # if q_node != ("CARTERET",)#("GASTON",)
        #     continue
        # end
        if q_node != ("NEW HANOVER",)
            continue
        end
        coarse_oriented_nbrs = Vector{Int}(undef, 0)
        brd_fine_node = coarse_nm_to_brd_fine_nd[q_node]
        oriented_nbrs = fine_node_attributes[brd_fine_node][oriented_nbrs_col]
        oriented_nbrs_int_huh, int_to_ext = build_oriented_nbrs_int_huh(
            q_node,
            oriented_nbrs,
            fine_node_attributes,
            partition
            )
        cur_node = brd_fine_node
        start_trace = nothing
        #start_trace_node = oriented_nbrs[start_trace]
        # println()
        # for i = 1:10
        #     println("start run")
        # end
        # println(cur_node)
        #println(fine_node_attributes[brd_fine_node]["pct21"], " start_trace ", fine_node_attributes[oriented_nbrs[start_trace]]["pct21"])
        cur_trace = int_to_ext[1]
        breakcounter = 0
        #println(cur_node, oriented_nbrs, oriented_nbrs_int_huh, int_to_ext)
        #println("start ", start_trace, " ", cur_node)
        start_node = 0
        while true
            fnd_int = false
            nbr = oriented_nbrs[cur_trace]
            while !fnd_int
                cur_trace = mod(cur_trace, length(oriented_nbrs)) + 1
                #println("at ",fine_node_attributes[oriented_nbrs[cur_trace]]["pct21"])
                #println("inner loop ", cur_trace, " ", oriented_nbrs_int_huh[cur_trace],
                #       " ", oriented_nbrs[cur_trace], " cur_node ", fine_node_attributes[cur_node]["pct21"]),
                       #" ", coarse_oriented_nbrs)
                if oriented_nbrs_int_huh[cur_trace]
                    fnd_int = true
                    #println("stop")
                else
                    nbr = oriented_nbrs[cur_trace]
                    if nbr == 0
                        nbr_node = 0
                    else
                        nbr_node = Tuple([fine_node_attributes[nbr][k] for k in partition])
                        nbr_id = coarse_partition_to_id[nbr_node]
                        push!(coarse_oriented_nbrs, nbr_id)
                        #println("exterior in county ",nbr_node)
                    end
                end
            end
            prev_node = cur_node
            cur_node = oriented_nbrs[cur_trace]
            entry_node = oriented_nbrs[mod(cur_trace-2,length(oriented_nbrs))+1] #last exterior node
            println("entry_node ", entry_node, " old oriented_nbrs ", oriented_nbrs)
            #println("final exterior node ",fine_node_attributes[entry_node]["pct21"])
            if breakcounter == 0
                start_node = oriented_nbrs[cur_trace]
                coarse_oriented_nbrs = []
            end
            oriented_nbrs = fine_node_attributes[cur_node][oriented_nbrs_col]
            oriented_nbrs_int_huh = build_oriented_nbrs_int_huh(
                q_node,
                oriented_nbrs,
                fine_node_attributes,
                partition
                )[1]
            potential_trace_start = findall(jj->jj==prev_node,oriented_nbrs)
            trace_start = 0
            for jj in potential_trace_start
                if !oriented_nbrs_int_huh[mod(jj, length(oriented_nbrs)) + 1]
                    trace_start = jj
                    break
                end
            end
            println("current int ", trace_start)
            println("prev_node = cur_node ", prev_node, ", ", cur_node)
            println("new oriented_nbrs ", oriented_nbrs)
            println("circshift(oriented_nbrs,-trace_start)", circshift(oriented_nbrs,-trace_start))
            println()
            cur_trace = mod(minimum(findall(jj->jj==entry_node,circshift(oriented_nbrs,-trace_start)))+trace_start-2,length(oriented_nbrs))+1
            # println("nbrs ",[fine_node_attributes[v]["pct21"] for v in oriented_nbrs])
            # println("enter at ",fine_node_attributes[oriented_nbrs[cur_trace]]["pct21"])
            if (cur_node == start_node && cur_trace == start_trace) || breakcounter >100
                break
            end
            if breakcounter == 0
                start_trace = cur_trace
            end
            breakcounter += 1
        end
        q_id = coarse_partition_to_id[q_node]
        repetition_free_coarse_oriented_nbrs = Vector{Int}(undef, 0)
        if coarse_oriented_nbrs[1] != coarse_oriented_nbrs[length(coarse_oriented_nbrs)]
            push!(repetition_free_coarse_oriented_nbrs, coarse_oriented_nbrs[1])
        end
        for i = 2:length(coarse_oriented_nbrs)
            if coarse_oriented_nbrs[i] != coarse_oriented_nbrs[i-1]
                push!(repetition_free_coarse_oriented_nbrs, coarse_oriented_nbrs[i])
            end
        end
        coarse_node_attributes[q_id][oriented_nbrs_col] = repetition_free_coarse_oriented_nbrs
        println(q_node,coarse_partition_to_id[q_node],repetition_free_coarse_oriented_nbrs)
    end
end

function build_oriented_nbrs_int_huh( #create b for a single node
    q_node::Tuple{String},
    oriented_nbrs::Vector{Any},
    fine_node_attributes::Vector{Dict{String,Any}},
    partition::Vector{String}
)
    oriented_nbrs_int_huh = Vector{Bool}(undef, length(oriented_nbrs))
    int_to_ext = Vector{Int}(undef,0)
    for (ii, o_nbr) in enumerate(oriented_nbrs)
        if o_nbr < 1
            oriented_nbrs_int_huh[ii] = false
        else
            o_nbr_name = Tuple([fine_node_attributes[o_nbr][k]
                                for k in partition])
            oriented_nbrs_int_huh[ii] = (o_nbr_name == q_node)
        end
        if (ii > 1 && oriented_nbrs_int_huh[ii-1] &&
            !oriented_nbrs_int_huh[ii])
            push!(int_to_ext, ii-1)
        end
    end
    if (oriented_nbrs_int_huh[length(oriented_nbrs)] &&
        !oriented_nbrs_int_huh[1])
        push!(int_to_ext, length(oriented_nbrs))
    end
    return oriented_nbrs_int_huh, int_to_ext
end

""""""
function build_quotient_nodes(
    base_graph::BaseGraph, levels::Vector{String},
    partition::Vector{String},
    id_to_partition::Vector{Tuple{Vararg{String}}},
    partition_to_id::Dict{Tuple{Vararg{String}},Int}
)
    num_nodes = length(id_to_partition)
    node_attributes = Vector{Dict{String,Any}}(undef, num_nodes)

    for q_id = 1:num_nodes
        node_attributes[q_id] = Dict{String,Any}()
    end

    oriented_nbrs_col = base_graph.oriented_nbrs_col

    for node in base_graph.node_attributes
        q_node = Tuple([node[k] for k in partition])
        q_id = partition_to_id[q_node]

        for part in partition
            node_attributes[q_id][part] = node[part]
        end

        add_keys = Set{String}()
        for attribute in keys(node)
            if attribute ∈ levels || attribute == oriented_nbrs_col
                continue
            end
            push!(add_keys, attribute)
        end

        if length(node_attributes[q_id]) == length(partition)
            for attribute in add_keys
                if typeof(node[attribute]) <: String
                    node_attributes[q_id][attribute] = Set([node[attribute]])
                else
                    node_attributes[q_id][attribute] = deepcopy(node[attribute])
                end
            end
        else
            for attribute in add_keys
                if typeof(node[attribute]) <: Dict
                    for (key,val) in node[attribute]
                        if haskey(node_attributes[q_id][attribute], key)
                            node_attributes[q_id][attribute][key] += val
                        else
                            node_attributes[q_id][attribute][key] = val
                        end
                    end
                elseif typeof(node_attributes[q_id][attribute]) <: Set
                    if !(typeof(node[attribute]) <: Set)
                        push!(node_attributes[q_id][attribute], node[attribute])
                    else
                        union!(node_attributes[q_id][attribute], 
                               node[attribute])
                    end
                else
                    node_attributes[q_id][attribute] += node[attribute]
                end
            end
        end
    end

    if oriented_nbrs_col != nothing
        #return node_attributes, partition, base_graph,levels, id_to_partition, partition_to_id
        build_oriented_coarse_nbrs!(node_attributes, partition, base_graph,
                                    levels, id_to_partition, partition_to_id)
    end


    return node_attributes
end

""""""
function build_quotient_edges(
    base_graph::BaseGraph, levels::Vector{String},
    partition::Vector{String},
    id_to_partition::Vector{Tuple{Vararg{String}}},
    partition_to_id::Dict{Tuple{Vararg{String}},Int}
)::Dict{Set{Int},Dict{String,Any}}
    edge_attributes = Dict{Set{Int},Dict{String,Any}}()

    base_node_attributes = base_graph.node_attributes
    for (edge, attributes) in base_graph.edge_attributes
        n_ind1, n_ind2 = edge
        q_node1 = Tuple([base_node_attributes[n_ind1][k] for k in partition])
        q_node2 = Tuple([base_node_attributes[n_ind2][k] for k in partition])
        q_ind1, q_ind2 = partition_to_id[q_node1], partition_to_id[q_node2]
        if q_ind1 == q_ind2
            continue
        end

        q_edge = Set([q_ind1, q_ind2])
        if haskey(edge_attributes, q_edge)
            for (attribute, value) in attributes
                if typeof(value) <: Real
                    edge_attributes[q_edge][attribute] += value
                elseif typeof(value) <: String
                    @assert edge_attributes[q_edge][attribute] == value
                else
                    @assert false
                end
            end
        else
            edge_attributes[q_edge] = Dict{String,Any}()
            for (attribute, value) in attributes
                edge_attributes[q_edge][attribute] = value
            end
        end
    end
    return edge_attributes
end

""""""
function build_quotient_graph(
    base_graph::BaseGraph, levels::Vector{String},
    partition::Vector{String},
    id_to_part::Vector{Tuple{Vararg{String}}},
    part_to_id::Dict{Tuple{Vararg{String}},Int}
)#::BaseGraph
    node_attributes = build_quotient_nodes(base_graph, levels, partition,
                                           id_to_part, part_to_id)
    edge_attributes = build_quotient_edges(base_graph, levels, partition,
                                           id_to_part, part_to_id)
    # @show edge_attributes
    num_nodes = length(id_to_part)
    num_edges = length(edge_attributes)
    total_pop = base_graph.total_pop

    simple_graph = SimpleWeightedGraph(length(id_to_part))
    edge_weights = base_graph.edge_weights

    for edge in keys(edge_attributes)
        n1, n2 = edge
        weight = edge_attributes[Set([n1, n2])][edge_weights]
        add_edge!(simple_graph, n1, n2)
        simple_graph.weights[n1, n2] = weight
        simple_graph.weights[n2, n1] = weight
    end

    return BaseGraph(
        num_nodes,
        num_edges,
        total_pop,
        base_graph.pop_col,
        base_graph.bpop_col,
        base_graph.vap_col,
        base_graph.bvap_col,
        edge_weights,
        base_graph.area_col,
        base_graph.node_border_col,
        base_graph.edge_perimeter_col,
        base_graph.oriented_nbrs_col,
        base_graph.mcd_col,
        simple_graph,
        node_attributes,
        edge_attributes,
    )
end


""""""
function build_qnode_coarse_to_fine_graphs(
    coarse_partition_to_ids::Dict{Tuple{Vararg{String}},Int},
    fine_partition_to_ids::Dict{Tuple{Vararg{String}},Int},
    base_graph::BaseGraph
)::Vector{SubGraph}
    num_coarse_nodes = length(coarse_partition_to_ids)
    fine_node_assignments = Vector{Set{Int}}([Set{Int}()
                                     for ii=1:num_coarse_nodes])

    for (dn_part, dn_id) in fine_partition_to_ids
        un_part = dn_part[1:length(dn_part)-1]
        un_id = coarse_partition_to_ids[un_part]
        push!(fine_node_assignments[un_id], dn_id)
    end

    subGraphs = Vector{SubGraph}(undef, num_coarse_nodes)

    for coarse_node_id = 1:num_coarse_nodes
        sub_nodes = collect(fine_node_assignments[coarse_node_id])
        sg, vmap = induced_subgraph(base_graph.simple_graph, sub_nodes)
        subGraphs[coarse_node_id] = SubGraph(sg, vmap)
    end
    return subGraphs
end

""""""
function build_fine_neighbor(
    node_name::Tuple{Vararg{String}},
    subgraph::SubGraph,
    graphs_by_level_view::AbstractVector{BaseGraph},
    coarse_to_fine_graphs_view::AbstractVector{Vector{SubGraph}},
    id_to_partitions_view::AbstractVector{Vector{Tuple{Vararg{String}}}}
)
    remaining_levels = length(graphs_by_level_view)
    @assert remaining_levels == length(coarse_to_fine_graphs_view)
    @assert remaining_levels == length(id_to_partitions_view)

    fine_neighbors = Dict{Int, Any}()
    coarse_level = length(node_name)

    for sub_nbr_ind = 1:nv(subgraph.graph)
        nbr_ind = subgraph.vmap[sub_nbr_ind]
        is_nbr = false
        simple_graph = graphs_by_level_view[1].simple_graph
        for nbr_nbr_ind in neighbors(simple_graph, nbr_ind)
            nbr_nbr_name = id_to_partitions_view[1][nbr_nbr_ind][1:coarse_level]
            if node_name == nbr_nbr_name
                is_nbr = true
                break
            end
        end
        if is_nbr
            if remaining_levels == 1
                fine_neighbors[nbr_ind] = Dict{Int, Any}()
            else
                ssbgrph = coarse_to_fine_graphs_view[1][nbr_ind]
                gblv = view(graphs_by_level_view, 2:remaining_levels)
                c2fgv = view(coarse_to_fine_graphs_view, 2:remaining_levels)
                id2pv = view(id_to_partitions_view, 2:remaining_levels)
                fn = build_fine_neighbor(node_name, ssbgrph, gblv, c2fgv, id2pv)
                fine_neighbors[nbr_ind] = fn
            end
        end
    end
    return fine_neighbors
end

""""""
function build_fine_neighbors(graphs_by_level_view,
                              coarse_to_fine_graphs_view,
                              id_to_partitions_view)
    graph = graphs_by_level_view[1]
    simple_graph = graph.simple_graph
    num_nodes = graph.num_nodes
    remaining_levels = length(graphs_by_level_view)

    fine_neighbors = Vector{Dict{Int, Any}}(undef, num_nodes)
    #println("num_nodes ", num_nodes)

    for node_ind = 1:num_nodes
        fine_neighbors[node_ind] = Dict{Int, Any}()
        node_partition = id_to_partitions_view[1][node_ind]
        for nbr_ind in neighbors(simple_graph, node_ind)
            #
            nbr_partition = id_to_partitions_view[1][nbr_ind]
            #println("here", id_to_partitions_view[1][node_ind], " -- ",
            #      id_to_partitions_view[1][nbr_ind], " :: ",
            #      node_ind, ", ", nbr_ind)
            #
            nbr_sub_graph = coarse_to_fine_graphs_view[1][nbr_ind]
            gblv = view(graphs_by_level_view, 2:remaining_levels)
            c2fgv = view(coarse_to_fine_graphs_view, 2:remaining_levels)
            id2pv = view(id_to_partitions_view, 2:remaining_levels)
            fn = build_fine_neighbor(node_partition, nbr_sub_graph, gblv,
                                     c2fgv, id2pv)
            fine_neighbors[node_ind][nbr_ind] = fn
        end
    end
    return fine_neighbors
end

""""""
function MultiLevelGraph(
    base_graph::BaseGraph,
    levels::Union{String, Vector{String}}
)#::MultiLevelGraph
    if typeof(levels) <: String
        levels = [levels]
    end
    num_levels = length(levels)
    id_to_partitions = Vector{Vector{Tuple{Vararg{String}}}}(undef,
                                                             num_levels)
    partition_to_ids = Vector{Dict{Tuple{Vararg{String}},Int}}(undef,
                                                               num_levels)
    build_quotient_ids!(id_to_partitions, partition_to_ids, base_graph, levels)
    graphs_by_level = Vector{BaseGraph}(undef, num_levels)
    graphs_by_level[num_levels] = base_graph
    for ii = num_levels-1:-1:1
        graphs_by_level[ii] = build_quotient_graph(graphs_by_level[ii+1],
                                                   levels,
                                                   levels[1:ii],
                                                   id_to_partitions[ii],
                                                   partition_to_ids[ii])
    end

    coarse_to_fine_graphs = Vector{Vector{SubGraph}}(undef, num_levels)
    for ii = num_levels-1:-1:1
        coarse_to_fine_graphs[ii] = build_qnode_coarse_to_fine_graphs(
                                                         partition_to_ids[ii],
                                                         partition_to_ids[ii+1],
                                                         graphs_by_level[ii+1])
    end

    fine_neighbors = Vector{Vector{Dict{Int,Any}}}(undef, num_levels-1)
    g2l_view = view(graphs_by_level, 1:num_levels)
    c2fgrphs_view = view(coarse_to_fine_graphs, 1:num_levels)
    id2part_view = view(id_to_partitions, 1:num_levels)
    for ii = 1:num_levels-1
        fine_neighbors[ii] = build_fine_neighbors(g2l_view, c2fgrphs_view,
                                                  id2part_view)
        remaining_levels = length(g2l_view)
        g2l_view = view(g2l_view, 2:remaining_levels)
        c2fgrphs_view = view(c2fgrphs_view, 2:remaining_levels)
        id2part_view = view(id2part_view, 2:remaining_levels)
    end

    #println("TODO: change to count mixed weights rather than connections only")
    # edge_weight_col
    mixed_nbr_weights = Dict{Set{Tuple{Vararg{String}}},Any}()
    for ii = 1:num_levels-1 #level of first node, if not census block
        for u = 1:graphs_by_level[ii].num_nodes #first node index
            simple_graph = graphs_by_level[ii].simple_graph
            for v in neighbors(simple_graph, u)
                build_connection_ew!(u,ii,v,ii,
                                        mixed_nbr_weights,
                                        fine_neighbors[ii][u][v],
                                        num_levels,
                                        graphs_by_level,
                                        id_to_partitions,
                                        base_graph.edge_weights)
            end
        end
    end

    if ((base_graph.num_nodes!=length(id_to_partitions[end])) ||
        (base_graph.num_nodes!=length(keys(partition_to_ids[end]))))
        println("Base graph not recovered level label.")
        println(" Possibly due to non-unique label")
        @assert false
    end

    return MultiLevelGraph(
        num_levels,
        levels,
        coarse_to_fine_graphs,
        fine_neighbors,
        mixed_nbr_weights,
        id_to_partitions,
        partition_to_ids,
        graphs_by_level
    )
end
"""
"""
function MultiLevelGraph(
    filepath::AbstractString,
    pop_col::AbstractString,
    levels::Union{String, Vector{String}};
    inc_node_data::Set{String}=Set(),
    edge_weights::String="connections",
    bpop_col=nothing,
    vap_col=nothing,
    bvap_col=nothing,
    area_col=nothing,
    node_border_col=nothing,
    edge_perimeter_col=nothing,
    oriented_nbrs_col=nothing,
    mcd_col=nothing,
    adjacency::String="rook"
)

    base_graph = BaseGraph(filepath, pop_col,
                    inc_node_data=inc_node_data,
                    edge_weights=edge_weights,
                    bpop_col=bpop_col,
                    vap_col=vap_col,
                    bvap_col=bvap_col,
                    area_col=area_col,
                    node_border_col=node_border_col,
                    edge_perimeter_col=edge_perimeter_col,
                    oriented_nbrs_col=oriented_nbrs_col,
                    mcd_col=mcd_col,
                    adjacency=adjacency)
    return MultiLevelGraph(base_graph, levels)
end
"""
compute total connection from edge_weights
"""
function build_connection_ew!(coarse_node_id,
    coarse_level,
    fine_node_id,
    fine_level,
    mixed_nbr_connection,
    finer_nbrs_dict, #fine neighbors of coarse node in and including fine node
    num_levels,
    gbl,#graphs by level
    idtp,
    edge_weights)
    coarse_name = idtp[coarse_level][coarse_node_id]
    if length(finer_nbrs_dict) == 0
        nbr_list = filter(i->(idtp[fine_level][i][coarse_level]==coarse_name[coarse_level])
                        ,neighbors(gbl[num_levels].simple_graph, fine_node_id))
        mixed_nbr_connection[Set([coarse_name,idtp[fine_level][fine_node_id]])]= sum(gbl[num_levels].edge_attributes[Set([fine_node_id,v])][edge_weights] for v in nbr_list)
    else
        nbr_length = 0.0
        for finer_node_id in keys(finer_nbrs_dict)
            build_connection_ew!(coarse_node_id,
                coarse_level,
                finer_node_id,
                fine_level+1,
                mixed_nbr_connection,
                finer_nbrs_dict[finer_node_id], #fine neighbors of u in and including v
                num_levels,
                gbl,
                idtp,
                edge_weights)
            nbr_length += mixed_nbr_connection[Set([coarse_name,idtp[fine_level+1][finer_node_id]])]
        end
        mixed_nbr_connection[Set([coarse_name,idtp[fine_level][fine_node_id]])] = nbr_length
    end
    return mixed_nbr_connection
end
"""
compute total connection length rather than count
"""
function build_connection_length!(coarse_node_id,
    coarse_level,
    fine_node_id,
    fine_level,
    mixed_nbr_connection,
    finer_nbrs_dict, #fine neighbors of coarse node in and including fine node
    num_levels,
    gbl,#graphs by level
    idtp)
    coarse_name = idtp[coarse_level][coarse_node_id]
    if length(finer_nbrs_dict) == 0
        nbr_list = filter(i->(idtp[fine_level][i][coarse_level]==coarse_name[coarse_level])
                        ,neighbors(gbl[num_levels].simple_graph, fine_node_id))
        mixed_nbr_connection[Set([coarse_name,idtp[fine_level][fine_node_id]])]= sum(gbl[num_levels].edge_attributes[Set([fine_node_id,v])]["length"] for v in nbr_list)
    else
        nbr_length = 0.0
        for finer_node_id in keys(finer_nbrs_dict)
            build_connection_length!(coarse_node_id,
                coarse_level,
                finer_node_id,
                fine_level+1,
                mixed_nbr_connection,
                finer_nbrs_dict[finer_node_id], #fine neighbors of u in and including v
                num_levels,
                gbl,
                idtp)
            nbr_length += mixed_nbr_connection[Set([coarse_name,idtp[fine_level+1][finer_node_id]])]
        end
        mixed_nbr_connection[Set([coarse_name,idtp[fine_level][fine_node_id]])] = nbr_length
    end
    return mixed_nbr_connection
end


function build_connection_count!(coarse_node_id,
    coarse_level,
    fine_node_id,
    fine_level,
    mixed_nbr_connection,
    finer_nbrs_dict, #fine neighbors of coarse node in and including fine node
    num_levels,
    gbl,#graphs by level
    idtp)
    coarse_name = idtp[coarse_level][coarse_node_id]
    if length(finer_nbrs_dict) == 0
        nbr_count = count(i->(idtp[fine_level][i][coarse_level]==coarse_name[coarse_level])
                        ,neighbors(gbl[num_levels].simple_graph, fine_node_id))#change to map based on other edge attribute, e.g. border length
        mixed_nbr_connection[Set([coarse_name,idtp[fine_level][fine_node_id]])] = nbr_count
    else
        nbr_count = 0
        for finer_node_id in keys(finer_nbrs_dict)
            build_connection_count!(coarse_node_id,
                coarse_level,
                finer_node_id,
                fine_level+1,
                mixed_nbr_connection,
                finer_nbrs_dict[finer_node_id], #fine neighbors of u in and including v
                num_levels,
                gbl,
                idtp)
            nbr_count += mixed_nbr_connection[Set([coarse_name,idtp[fine_level+1][finer_node_id]])]
        end
        mixed_nbr_connection[Set([coarse_name,idtp[fine_level][fine_node_id]])] = nbr_count
    end
    return mixed_nbr_connection
end
"""
    get_attributes(nodes::Array{Any, 1})

*Returns* an array of dicts `attributes` of length `length(nodes)` where
the attributes of the `nodes[i]` is at `attributes[i]` as a dictionary.
"""
function get_attributes(nodes::Array{Any,1})
    attributes = Array{Dict{String,Any}}(undef, length(nodes))
    for (index, node) in enumerate(nodes)
        attributes[index] = node
    end
    return attributes
end

""""""
function sum_attribute(
    graph::MultiLevelGraph,
    node_set::Dict{Tuple{Vararg{String}}, Any},
    attribute_name::String,
    level::Int=1
)::Real
    attribute = 0
    for node in keys(node_set)
        simple_graph = graph.graphs_by_level[level].simple_graph
        node_id = graph.partition_to_ids[level][node]
        if !intact_huh(node_set[node])
            attribute += sum_attribute(graph, node_set[node],
                                       attribute_name, level+1)
        else
            full_graph = graph.graphs_by_level[level]
            attribute += full_graph.node_attributes[node_id][attribute_name]
        end
    end

    return attribute
end

""""""
function induced_subgraph_edges(graph::MultiLevelGraph, vlist::Array{Int,1})::Array{Int,1}
    print("TODO (and think about)")
end

""""""
function get_subgraph_population(graph::MultiLevelGraph, nodes::BitSet)::Int
    print("TODO")
end

