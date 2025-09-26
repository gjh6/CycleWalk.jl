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
    @show weight.value, cur_step, total_steps, values(measure.weights), e1, e2
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
    writer::Union{Writer, Nothing}=nothing,
    run_diagnostics::RunDiagnostics=RunDiagnostics()
) where T <: Real
    base_measure = deepcopy(measure)
    modify_measure!(base_measure, 0, 1)

    outer_steps = Int(total_steps/base_steps_per_sample)
    @show outer_steps
    for ii = 1:outer_steps
        run_metropolis_hastings!(partition, proposal, base_measure, 
                                 base_steps_per_sample, rng)
        partition_to_anneal = deepcopy(partition)
        log_weight = MutableFloat(0.0)
        @show "starting annealing"
        modify_measure!(measure, 0, 1)
        run_metropolis_hastings!(partition_to_anneal, proposal, measure, 
                                 steps_per_annealing, rng; writer=writer,
                                 output_freq=steps_per_annealing*2,
                                 output_initial=false,
                                 prestepf=track_weight_and_modify_measure!,
                                 prestepargs=(partition, steps_per_annealing, 
                                              log_weight, measure, 
                                              modify_measure!),
                                 weight=log_weight)
        # break
    end
end
