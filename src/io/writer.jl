
AtlasParam=Dict{String, Any}
MapParam=Dict{String, Any}

mutable struct Writer
    atlas::Atlas{AtlasParam}
    map_param::MapParam
    map_output_data::Dict{String, Function}
    output_districting::Bool
    node_map::Dict{Tuple{Vararg{String}}, Int}
    node_field::String
    weight_type::DataType
end



function Writer(
    measure::Measure,
    constraints::Dict{Type{T} where T<:AbstractConstraint, AbstractConstraint},
    partition::LinkCutPartition,
    output_file_path::String;
    output_districting=true,
    description::String="",
    time_stamp=string(Dates.now()),
    io_mode::String="w",
    additional_parameters::Dict{String, Any}=Dict{String,Any}(),
    weight_type::DataType=Int64
    # proposal_diagnostics::Dict=Dict()
)
    graph = partition.graph
    scores = collect(measure.scores)
    energies = [measure.descriptions[e] for e in scores]
    weights = [measure.weights[e] for e in scores]
    atlasParam=AtlasParam("energies"=>energies, "energy weights"=>weights,
                          "districts"=>partition.num_dists)

    if haskey(constraints, PopulationConstraint)
        min_pop = constraints[PopulationConstraint].min_pop
        max_pop = constraints[PopulationConstraint].max_pop
        atlasParam["population bounds"] = [min_pop, max_pop]
    end

    f = @__FILE__
    #projdir = "/"*joinpath(split(f, "/")[1:end-1])
    versionNumber = pkgversion(CycleWalk)
    versionString = string(versionNumber)
    atlasParam["package.version"] = "CycleWalk v"*versionString

    for (key,val) in additional_parameters
        atlasParam[key] = val
    end
    
    # to add to atlasParam
    # other constraints

    dir = dirname(output_file_path)
    # split_path = split(output_file_path, "/")
    # dir = join(split_path[1:length(split_path)-1], "/")
    if !isdir(dir)
        mkpath(dir)
    end

    atlasHeader = AtlasHeader(description, time_stamp, AtlasParam, MapParam;weightType=weight_type)
    io = smartOpen(output_file_path, io_mode)
    newAtlas(io, atlasHeader, atlasParam)

    atlas = Atlas{AtlasParam}(io, description, time_stamp, atlasParam, MapParam,weight_type)
    map_output_data = Dict{String, Function}()

    node_map = get_node_map(partition.node_col, partition)

    return Writer(atlas, MapParam(), map_output_data, output_districting, 
                  node_map, partition.node_col,weight_type)#, proposal_diagnostics)
end

mutable struct Batch_Writer
    writer::Writer
    mapChannel::Channel
    listenerTask::Task
end

function Batch_Writer(
    measure::Measure,
    constraints::Dict{Type{T} where T<:AbstractConstraint, AbstractConstraint},
    partition::LinkCutPartition,
    output_file_path::String;
    output_districting::Bool,
    description::String,
    time_stamp::String,
    io_mode::String,
    additional_parameters::Dict{String, Any},
    weight_type::DataType=Int64
    )
    writer = Writer(measure, constraints, partition, output_file_path,output_districting,
                    description,time_stamp,
                    io_mode,additional_parameters,weight_type)
    mapChannel = Channel{Map}(typemax(Int)) # buffer size of Int max
    listenerTask = @async begin
        
    end
    return Batch_Writer(writer, mapChannel, listenerTask)



end

function push_writer!(
    writer::Writer,
    get_data::Function; 
    desc::Union{String, Nothing}=nothing
)
    if desc == nothing
        desc = string(get_data)
    end
    writer.map_output_data[desc] = get_data
end

function close_writer(writer::Writer)
    close(writer.atlas.io)
end

function get_node_map(
    node_field::String, 
    partition::LinkCutPartition,
    node_map::Union{Nothing, Dict{Tuple{Vararg{String}}, Int}} = nothing
)
    if node_map === nothing
        node_map = Dict{Tuple{Vararg{String}}, Int}()
    end

    graph = partition.graph
    for ni = 1:graph.num_nodes
        node_name = Tuple([graph.node_attributes[ni][node_field]])
        node_map[node_name] = partition.node_to_dist[ni]
    end
    return node_map
end


function get_node_map!(
    writer::Writer, 
    partition::LinkCutPartition
)
    return get_node_map(writer.node_field, partition, writer.node_map)
end

""""""
function prepare_output_map_and_writer!(writer::Writer,
    partition::LinkCutPartition,
    step::Int,
    run_diagnostics::RunDiagnostics,
    weight::Union{MutableFloat, Float64},   # This looks like it only writes Floats now.
    count::Int
 )::Map
    for (desc, f) in writer.map_output_data
        writer.map_param[desc] = f(partition)
    end
    for (rd_desc, proposal_diagnostics) in values(run_diagnostics)
        for (pd_desc, proposal_diagnostic) in proposal_diagnostics
            desc = "("*join([rd_desc, pd_desc], ",")*")"
            writer.map_param[desc] = proposal_diagnostic.data_vec
        end
    end

    # @show writer.map_param
########## get_map() or something
    if !writer.output_districting
        d = Dict{Tuple{Vararg{String}}, Int}()
        map = Map{MapParam,writer.weight_type}("step"*string(step-count), d, weight, 
                            writer.map_param)
    else
        map = Map{MapParam,writer.weight_type}("step"*string(step-count), 
                            get_node_map!(writer, partition), 
                            weight, writer.map_param)
    end
    return map
end
##########

""""""
function output(
    partition::LinkCutPartition,
    measure::Measure,
    step::Int,
    output_freq::Int,
    writer::Union{Writer, Nothing},
    run_diagnostics::RunDiagnostics=RunDiagnostics();
    weight::Union{MutableFloat, Float64}=1.0,   # This looks like it only writes Floats now.
    count::Int=0
)
    if writer == nothing || mod(step, output_freq) != 0
        return
    end
    if typeof(weight) == MutableFloat
        weight = weight.value
    end
    ##
    map=prepare_output_map_and_writer!(writer,partition, step,
            run_diagnostics,weight,count)
    ##

    try
        addMap(writer.atlas.io, map)
    catch e
        @show writer.map_param
        @show count
        # println("partition.node_to_district: ", partition.node_to_district)
        println("Could not add map to atlas")
        @assert false
    end
    reset_diagnostics!(run_diagnostics)
end