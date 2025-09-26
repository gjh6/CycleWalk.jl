""""""
function run_metropolis_hastings!(
    partition::LinkCutPartition,
    proposal::Union{Function,Vector{Tuple{T, Function}}},
    measure::Measure,
    steps::Union{Int,Tuple{Int,Int}},
    rng::AbstractRNG;
    writer::Union{Writer, Nothing}=nothing,
    output_freq::Int=250,
    run_diagnostics::RunDiagnostics=RunDiagnostics(),
    prestepf::Function=(x...)->nothing,
    prestepargs::Tuple=(),
    output_initial::Bool=true,
    weight::Union{Float64, MutableFloat}=1.0#,
    # poststepf!::Function=(x...)->nothing,
    # poststepargs::Tuple=()
) where T <: Real
    # Want this, but need to redefine it: precompute_node_tree_counts!(partition)
    check_proposals_weights(proposal)
    
    initial_step, final_step = set_step_bounds(steps)
    if (initial_step == 0 || initial_step == 1) && output_initial
        output(partition, measure, initial_step, 1, writer)
    end

    for step = initial_step:final_step
        prestepf(step, prestepargs...)
        proposal!, proposal_index = get_random_proposal(proposal, rng)
        proposal_diagnostics = get_proposal_diagnostics(run_diagnostics, 
                                                        proposal!)
        p, update = proposal!(partition, rng, diagnostics=proposal_diagnostics)
        if p == 0
            output(partition, measure, step, output_freq, writer, 
                   run_diagnostics; weight=weight)
            continue
        end

        p *= get_delta_energy(partition, measure, update)
        update_acceptance_ratio_diagnostic!(proposal_diagnostics, p)

        if rand(rng) < p
            update_partition!(partition, update)
        end
        output(partition, measure, step, output_freq, writer, 
               run_diagnostics; weight=weight)
    end
end

""""""
function update_partition!(
    partition::LinkCutPartition,
    update::Update{T}
) where T <: Int
    # @show distPair, links, cuts

    # @show "debugging 1"
    # for node_ind = 1:partition.graph.num_nodes
    #     node = partition.lct.nodes[node_ind]
    #     r = find_root!(node).vertex
    #     # if !(r in partition.district_roots)
    #     #     @show node_ind, r
    #     # end
    #     @assert r in partition.district_roots
    #     @assert haskey(partition.roots_to_district, r)
    # end

    # @show cuts, links
    # for cut in cuts
    #     @show find_root!(partition.lct.nodes[cut[1]]).vertex
    #     @show find_root!(partition.lct.nodes[cut[2]]).vertex
    # end
    # for link in links
    #     @show find_root!(partition.lct.nodes[link[1]]).vertex
    #     @show find_root!(partition.lct.nodes[link[2]]).vertex
    # end

    for cut in update.cuts
        evert!(partition.lct.nodes[cut[1]])
        cut!(partition.lct.nodes[cut[2]])
    end
    for link in update.links
        evert!(partition.lct.nodes[link[1]])
        link!(partition.lct.nodes[link[1]], partition.lct.nodes[link[2]])
    end

    distPair = update.changed_districts
    new_roots = (find_root!(partition.lct.nodes[update.cuts[1][1]]),
                 find_root!(partition.lct.nodes[update.cuts[1][2]]))
    # swap if needed
    l11node = partition.lct.nodes[update.links[1][1]]
    l11dist_cur = partition.node_to_dist[update.links[1][1]]
    r11_new = find_root!(l11node)
    r11_root_ind_new = (r11_new != new_roots[1]) + 1
    r11_root_ind_cur = (l11dist_cur != distPair[1]) + 1
    if (r11_root_ind_new == r11_root_ind_cur) ⊻ !update.swap_link11
        new_roots = (new_roots[2], new_roots[1])
    end


    # modify roots
    old_root1 = partition.district_roots[distPair[1]]
    old_root2 = partition.district_roots[distPair[2]]
    partition.district_roots[distPair[1]] = new_roots[1].vertex
    partition.district_roots[distPair[2]] = new_roots[2].vertex
    delete!(partition.roots_to_district, old_root1)
    delete!(partition.roots_to_district, old_root2)
    partition.roots_to_district[new_roots[1].vertex] = distPair[1]
    partition.roots_to_district[new_roots[2].vertex] = distPair[2]

    # new_cross_d_edg 
    del_keys = [k for k in keys(partition.cross_district_edges) 
                if (distPair[1] in k || distPair[2] in k) && 
                    !haskey(update.new_cross_d_edg, k)]
    for dk in del_keys
        delete!(partition.cross_district_edges, dk)
    end
    for (dists, edgeset) in update.new_cross_d_edg
        partition.cross_district_edges[dists] = edgeset
    end

    partition.node_to_dist .= partition.node_to_dist_update

    partition.identifier += 1
    for (key, eData) in partition.energy_data
        # @show update.changed_districts
        # @show "showing eData:", eData
        update_energy_data!(eData, partition, update)
        # @show "showing eData:", eData
    end
    
    # tpop = 0
    # for di = 1:partition.num_dists
    #     r = partition.lct.nodes[partition.district_roots[di]]
    #     pop = sum_cc(r, partition, partition.graph.pop_col)
    #     @show di, pop
    #     tpop += pop
    # end
    # @show tpop, partition.graph.total_pop
    # @show partition.district_roots
    # @show "debugging 2"
    # for node_ind = 1:partition.graph.num_nodes
    #     node = partition.lct.nodes[node_ind]
    #     r = find_root!(node).vertex
    #     # if !(r in partition.district_roots)
    #     #     @show node_ind, r
    #     # end
    #     @assert r in partition.district_roots
    #     @assert haskey(partition.roots_to_district, r)
    # end
    # @assert tpop == partition.graph.total_pop
end


""""""
@inline function set_step_bounds(steps::Union{Tuple{T,T}, T}) where T<:Int
    if typeof(steps)<:Tuple
        return steps
    else
        return 1, steps
    end
end


""""""
@inline function get_random_proposal(
    proposal::Union{Function, Vector{Tuple{T, Function}}},
    rng::AbstractRNG
) where T <: Real
    if !(typeof(proposal) <: Vector)
        return proposal, 1
    end
    proposal_weights = [proposal[i][1] for i = 1:length(proposal)]
    index = findfirst(cumsum(proposal_weights) .> rand(rng))
    return proposal[index][2], index
end


""""""
@inline function check_proposals_weights(
    proposal::Union{Function, Vector{Tuple{T, Function}}}
) where T<:Real
    if !(typeof(proposal) <: Vector)
        return
    end
    weight_sum = sum(proposal[i][1] for i = 1:length(proposal))
    if weight_sum != 1
        throw(
            ArgumentError(
                "Chance of choosing a proposal must sum to 1",
            ),
        )
    end
end

function post_step!(
    proposal_diagnostics::Union{ProposalDiagnostics, Nothing}, 
    acceptance_ratio::Float64, 
    partition::LinkCutPartition, 
    measure::Measure,
    writer::Writer,
    update::Update{T}
) where T<:Int
    update_acceptance_ratio_diagnostic!(proposal_diagnostics, acceptance_ratio)
    # for report_func in values(writer.map_output_data)
    # for (desc, report_func) in writer.map_output_data
    #     @show desc, report_func∈measure.scores, report_func
    # end
end

