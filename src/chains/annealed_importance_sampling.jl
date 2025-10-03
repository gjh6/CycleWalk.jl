function track_weight_and_modify_measure!(
    cur_step::Int,
    partition::LinkCutPartition,
    total_steps::Int,
    weight::MutableFloat,
    measure::Measure,
    modify_measure!::Function
)
    # println("weight before", weight.value)
    e1 = get_energy(partition, measure)
    modify_measure!(measure, cur_step, total_steps)
    e2 = get_energy(partition, measure)
    weight.value += e2-e1
    # println("weiweight.value)
    #@show weight.value, cur_step, total_steps, values(measure.weights), e1, e2
end

function run_annealed_importance_sampling!(
    partition::LinkCutPartition,
    proposal::Union{Function,Vector{Tuple{T, Function}}},
    measure::Measure,
    modify_measure!::Function,
    total_steps::Int,
    base_steps_per_sample::Int,
    steps_per_annealing::Int,
    rng::AbstractRNG;
    writer::Union{Writer, Async_Writer, Nothing}=nothing,
    run_diagnostics::RunDiagnostics=RunDiagnostics()
) where T <: Real
    base_measure = deepcopy(measure)
    modify_measure!(base_measure, 0, 1)

    outer_steps = Int(total_steps/base_steps_per_sample)
    @show outer_steps
    @sync begin
        for ii = 1:outer_steps
            run_metropolis_hastings!(partition, proposal, base_measure, 
                                    base_steps_per_sample, rng)
            partition_to_anneal = deepcopy(partition)
            log_weight = MutableFloat(0.0)
            println("starting annealing :",ii)
            thread_measure = deepcopy(base_measure)
            modify_measure!(thread_measure, 0, 1)
            thread_rng=deepcopy(rng)
            # push!(block,@spawn run_metropolis_hastings_return_map!())
            run_metropolis_hastings!(partition_to_anneal, proposal, measure, 
                                    steps_per_annealing, thread_rng; writer=writer,
                                    output_freq=steps_per_annealing,
                                    output_initial=false,
                                    prestepf=track_weight_and_modify_measure!,
                                    prestepargs=(partition, steps_per_annealing, 
                                                log_weight, thread_measure, 
                                                modify_measure!),
                                    weight=log_weight)
            # break
        end
    end
end

function run_metropolis_hastings_return_map!(
    partition::LinkCutPartition,
    proposal::Union{Function,Vector{Tuple{T, Function}}},
    measure::Measure,
    steps::Union{Int,Tuple{Int,Int}},
    rng::AbstractRNG;
    writer::Union{Writer, Nothing}=nothing,
    run_diagnostics::RunDiagnostics=RunDiagnostics(),
    prestepf::Function=(x...)->nothing,
    prestepargs::Tuple=(),
    weight::Union{Float64, MutableFloat}=1.0
)::Map where T <: Real
    try
        run_metropolis_hastings!(partition, proposal, measure, 
                                        steps, rng;
                                        prestepf=prestepf,
                                        prestepargs=prestepargs,
                                        weight=weight,output_initial=false)
        map=prepare_output_map_and_writer!(writer,partition, steps,
                run_diagnostics,weight,count)
        return map
    catch e
        @error "Error during MCMC run: $e"
    end
end

function run_batch_annealed_importance_sampling!(
    partition::LinkCutPartition,
    proposal::Union{Function,Vector{Tuple{T, Function}}},
    measure::Measure,
    modify_measure!::Function,
    total_steps::Int,
    base_steps_per_sample::Int,
    steps_per_annealing::Int,
    rng::AbstractRNG;
    writer::Union{Async_Batch_Writer, Nothing}=nothing,
    run_diagnostics::RunDiagnostics=RunDiagnostics(),
    batch_size::Union{Int64,Symbol}=7
) where T <: Real
    base_measure = deepcopy(measure)
    modify_measure!(base_measure, 0, 1)

    outer_steps = Int(total_steps/base_steps_per_sample)
    @show outer_steps
    @sync begin
        for ii = 1:Int(ceil(outer_steps/batch_size))
            jobs = []
            for jj=1:batch_size
                run_metropolis_hastings!(partition, proposal, base_measure, 
                    base_steps_per_sample, rng)
                partition_to_anneal = deepcopy(partition)
                log_weight = MutableFloat(0.0)
                println("starting annealing :(",ii,",",jj,")")
                thread_measure = deepcopy(base_measure)
                modify_measure!(thread_measure, 0, 1)
                thread_rng=deepcopy(rng)
                t= @spawn run_metropolis_hastings_return_map!(
                        partition_to_anneal, proposal, thread_measure, 
                        steps_per_annealing, thread_rng; writer=writer.writer,
                        prestepf=track_weight_and_modify_measure!,
                        prestepargs=(partition_to_anneal, steps_per_annealing, 
                                        log_weight, thread_measure, modify_measure!),
                        weight=log_weight
                    )
                #@show t
                push!(jobs,t)
                #@show jobs
            end
            #@show jobs
            println("pushing block to channel")
            put!(writer.mapChannel,(ii,jobs))
        end
    end
end
