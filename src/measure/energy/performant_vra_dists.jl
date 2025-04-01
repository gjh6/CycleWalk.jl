mutable struct PerformantVRAData <: AbstractEnergyData
    identifier::Int64
    votes::Vector{Vector{Vector{Vector{Float64}}}}
    votes_update::Vector{Vector{Vector{Vector{Float64}}}}
end

function PerformantVRAData(
    partition::LinkCutPartition,
    elections::Vector{E}
) where {E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    identifier = partition.identifier - 1
    votes = Vector{Vector{Vector{Vector{Float64}}}}(undef, length(elections))
    vts_ud = Vector{Vector{Vector{Vector{Float64}}}}(undef, length(elections))
    for ei = 1:length(elections)
        num_elec_stages = length(elections[ei])
        votes[ei] = Vector{Vector{Vector{Float64}}}(undef, num_elec_stages)
        vts_ud[ei] = Vector{Vector{Vector{Float64}}}(undef, num_elec_stages)
        for si = 1:length(elections[ei])
            votes[ei][si] = Vector{Vector{Float64}}(undef, partition.num_dists)
            vts_ud[ei][si] = Vector{Vector{Float64}}(undef, partition.num_dists)
            num_candidates = length(elections[ei][si])
            for di = 1:partition.num_dists
                votes[ei][si][di] = Vector{Float64}(undef, num_candidates)
                vts_ud[ei][si][di] = -ones(Float64, num_candidates)
            end
        end
    end
    return PerformantVRAData(identifier, votes, vts_ud)
end


""""""
function set_vra_election_data!(
    votes::Vector{Vector{Vector{Vector{Float64}}}}, 
    partition::LinkCutPartition, 
    elections::Vector{E},
    node_to_dist::Vector{Int64},
    districts::Union{Tuple{Vararg{T}}, Vector{T}}
) where {T <: Int, E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    graph = partition.graph

    # zero out elections in districts to be computed
    for ei = 1:length(elections)
        for si = 1:length(elections[ei])
            for di in districts
                votes[ei][si][di] .= 0
            end
        end
    end
    
    # compute elections in districts to be computed
    for ni = 1:graph.num_nodes
        di = node_to_dist[ni]
        if di ∉ districts
            continue
        end
        n_attr = graph.node_attributes[ni]
        for ei = 1:length(elections)
            for si = 1:length(elections[ei])
                for ci = 1:length(elections[ei][si])
                    candidate = elections[ei][si][ci]
                    votes[ei][si][di][ci] += n_attr[candidate]
                end
            end
        end
    end
end


""""""
function get_vra_performance_votes(
    partition::LinkCutPartition,
    elections::Vector{E},
    districts::Union{Tuple{Vararg{T}}, Vector{T}}
        =collect(1:partition.num_dists);
    update::Union{Update, Nothing} = nothing
) where {T <: Int, E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    if !haskey(partition.energy_data, (PerformantVRAData, (elections,)))
        partition.energy_data[(PerformantVRAData, (elections,))] = 
            PerformantVRAData(partition, elections)
    end
    vra_data = partition.energy_data[(PerformantVRAData, (elections,))]

    votes = nothing

    if update === nothing
        votes = vra_data.votes
        if partition.identifier != vra_data.identifier
            node_to_dist = partition.node_to_dist    
            set_vra_election_data!(votes, partition, elections, node_to_dist, 
                                   districts)
            vra_data.identifier = partition.identifier
        end

    else
        @assert partition.identifier == vra_data.identifier
        votes = vra_data.votes_update
        recursive_copy!(votes, vra_data.votes)
        node_to_dist = partition.node_to_dist_update
        set_vra_election_data!(votes, partition, elections, node_to_dist,
                               update.changed_districts)
    end

    @assert votes !== nothing
    return votes
end


""""""
function get_vra_performance_margins(
    partition::LinkCutPartition,
    elections::Vector{E},
    districts::Union{Tuple{Vararg{T}}, Vector{T}}
        =collect(1:partition.num_dists);
    update::Union{Update, Nothing} = nothing
) where {T <: Int, E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    
    votes = get_vra_performance_votes(partition, elections, districts; 
                                      update=update)

    margins = Vector{Vector{Vector{Float64}}}(undef, length(elections))
    for ei = 1:length(elections)
        margins[ei] = Vector{Vector{Float64}}(undef, length(elections[ei]))
        for si = 1:length(elections[ei])
            margins[ei][si] = Vector{Float64}(undef, partition.num_dists)
            for di = 1:partition.num_dists
                vra_cand_votes = votes[ei][si][di][1]
                max_votes = maximum(votes[ei][si][di])
                if max_votes == vra_cand_votes
                    margins[ei][si][di] = 0
                else
                    # @show length(margins[ei]), length(margins[ei][si])
                    margins[ei][si][di] = (max_votes-vra_cand_votes)/max_votes
                end
            end
        end
    end

    return margins
end

""""""
function get_performant_vra_score(
    partition::LinkCutPartition,
    elections::Vector{E},
    target_districts::Int64,
    districts::Union{Tuple{Vararg{T}}, Vector{T}}
        =collect(1:partition.num_dists);
    weights::Vector{Float64} = ones(Float64, length(elections)),
    update::Union{Update, Nothing} = nothing
)::Float64 where {T <: Int, E <: Tuple{Tuple{Vararg{String}}, 
                                       Tuple{Vararg{String}}}}
    @assert length(weights) == length(elections)
    margins = get_vra_performance_margins(partition, elections, districts;
                                          update=update)

    score = 0
    scores = zeros(Float64, partition.num_dists)
    for ei = 1:length(elections)
        scores .= 0
        for di = 1:partition.num_dists 
            # @show ei, di, length(margins), length(margins[ei][1]), length(margins[ei][2])
            if margins[ei][1][di] > 0 || margins[ei][2][di] > 0
                scores[di] = margins[ei][1][di] + margins[ei][2][di]
            end
        end
        sort!(scores)
        score += weights[ei]*sum(scores[1:target_districts])
    end
    return score
end

function build_performant_vra_score(
    graph::BaseGraph,
    elections::Vector{E};
    weights::Vector{Float64} = ones(Float64, length(elections)),
    target_districts::Union{Nothing, Int64} = nothing,
    num_dists::Union{Nothing, Int64} = nothing,
    total_pop_col::Union{Nothing, String} = nothing,
    mino_pop_col::Union{Nothing, String, Tuple{Vararg{String}}} = nothing
) where {E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    if target_districts === nothing
        @assert total_pop_col !== nothing
        @assert mino_pop_col !== nothing
        @assert num_dists !== nothing
        target_districts = get_target_vra_districts(graph, num_dists, 
                                                    total_pop_col, mino_pop_col)
    end
    f(p, d=collect(1:p.num_dists); update=nothing) = 
        get_performant_vra_score(p, elections, target_districts, d; 
                                 weights=weights, update=update)
    return f
end

""""""
function get_performant_vra_report(
    partition::LinkCutPartition,
    elections::Vector{E},
    target_districts::Int64,
    districts::Union{Tuple{Vararg{T}}, Vector{T}}
        =collect(1:partition.num_dists);
    weights::Vector{Float64} = ones(Float64, length(elections)),
    update::Union{Update, Nothing} = nothing
) where {T <: Int, E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    margins = get_vra_performance_margins(partition, elections, districts;
                                          update=update)

    # dist_counts = zeros(Int64, length(elections))
    # for ei = 1:length(elections)
    #     for di = 1:partition.num_dists 
    #         if margins[ei][1][di] == 0 && margins[ei][2][di] == 0
    #             dist_counts[ei] += 1
    #         end
    #     end
    # end
    # return [[target_districts]; dist_counts]
    dist_counts = Dict{Float64, Int64}()
    for di = 1:partition.num_dists 
        count = 0
        for ei = 1:length(elections)
            if margins[ei][1][di] == 0 && margins[ei][2][di] == 0
                count += weights[ei]
            end
        end
        if count > 0 
            if haskey(dist_counts, count)
                dist_counts[count] += 1
            else
                dist_counts[count] = 1
            end
        end
    end
    return [[target_districts]; dist_counts]
end


function build_performant_vra_report(
    graph::BaseGraph,
    elections::Vector{E};
    weights::Vector{Float64} = ones(Float64, length(elections)),
    target_districts::Union{Nothing, Int64} = nothing,
    num_dists::Union{Nothing, Int64} = nothing,
    total_pop_col::Union{Nothing, String} = nothing,
    mino_pop_col::Union{Nothing, String, Tuple{Vararg{String}}} = nothing
) where {E <: Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}}}
    if target_districts === nothing
        @assert total_pop_col !== nothing
        @assert mino_pop_col !== nothing
        @assert num_dists !== nothing
        target_districts = get_target_vra_districts(graph, num_dists, 
                                                    total_pop_col, mino_pop_col)
    end
    f(p, d=collect(1:p.num_dists); update=nothing) = 
        get_performant_vra_report(p, elections, target_districts, d; 
                                  weights=weights, update=update)
    return f
end

function get_target_vra_districts(
    graph::BaseGraph,
    num_dists::Int64,
    total_pop_col::String,
    mino_pop_col::String
)
    total_pop = operate_on_attribute(graph, total_pop_col, sum)
    mino_pop = operate_on_attribute(graph, mino_pop_col, sum)
    return Int(floor(mino_pop/total_pop*num_dists))
end


function update_energy_data!(
    eData::PerformantVRAData,
    partition::LinkCutPartition,
    update::Update{T}
) where {T<:Int}
    for ei = 1:length(eData.votes)
        for si = 1:length(eData.votes[ei])
            for di in update.changed_districts
                if minimum(eData.votes_update[ei][si][di]) == -1
                    set_all!(eData.votes_update, -1.0)
                    return
                end
            end
        end
    end
    eData.identifier = partition.identifier
    recursive_copy!(eData.votes, eData.votes_update)
    set_all!(eData.votes_update, -1.0)
end
