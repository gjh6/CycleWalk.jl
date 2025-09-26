mutable struct Measure
    weights::Dict{Function, Float64}
    scores::Set{Function}
    descriptions::Dict{Function, String}
    # TO DO: constraints::Dict{Type{T} where T<:AbstractConstraint, AbstractConstraint}
end

""""""
function Measure()
    scores = Set{Function}()
    weights = Dict{Function, Float64}()
    descriptions = Dict{Function, String}()
    return Measure(weights, scores, descriptions)
end

""""""
function push_energy!(
    measure::Measure,
    score::Function,
    weight::Real;
    desc::String=""
)
    @assert keys(measure.weights) == measure.scores
    @assert keys(measure.descriptions) == measure.scores
    push!(measure.scores, score)
    measure.weights[score] = weight
    if desc == ""
        desc = string(score)
    end
    measure.descriptions[score] = desc
end

""""""
function get_delta_energy(
    partition::LinkCutPartition, 
    measure::Measure,
    update::Update{T}
) where T <: Int
    score = 0
    changed_districts = update.changed_districts
    # @show get_energy(partition, measure, changed_districts)
    # @show get_energy(partition, measure, changed_districts, update)
    # @show changed_districts
    score += get_energy(partition, measure, changed_districts)
    score -= get_energy(partition, measure, changed_districts, update)
    # @show score, exp(score)
    return exp(score)
end

""""""
function get_energy(
    partition::LinkCutPartition, 
    measure::Measure,
    districts::Union{Tuple{Vararg{T}}, Vector{T}}
        = collect(1:partition.num_dists),
    update::Union{Update{T}, Nothing}=nothing
)::Float64 where T <: Int
    score = 0.0
    for energy in measure.scores
        weight = measure.weights[energy]
        if weight == 0
            continue
        end
        # tmp = energy(partition, districts; update=update)
        # @show weight, tmp
        score += weight*energy(partition, districts; update=update)
    end
    return score
end

# ################################## TO DO: EVERYTHING BELOW SHOULD GO ELSEWHERE
""""""
function get_cut_edge_count(partition::MultiLevelPartition)
    return get_cut_edge_weights(partition, "connections")
end

""""""
function get_cut_edge_perimeter(partition::MultiLevelPartition)
    edge_perimeter_col = partition.graph.graphs_by_level[1].edge_perimeter_col
    return get_cut_edge_weights(partition, edge_perimeter_col)
end

""""""
function get_cut_edge_sum(
    partition::LinkCutPartition;
    column::String="connections"
)
    graph = partition.graph
    total = 0
    for e in edges(graph.simple_graph)
        n1, n2 = src(e), dst(e)
        if partition.node_to_dist[n1] == partition.node_to_dist[n2]
            continue
        end
        total += graph.edge_attributes[Set([n1, n2])][column]
    end
    return total
end

