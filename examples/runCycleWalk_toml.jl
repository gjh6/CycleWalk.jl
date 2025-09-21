# Example script to run CycleWalk with parameters from a TOML file
#
# julia runCycleWalk_toml.jl toml/param_ct.toml
# julia runCycleWalk_toml.jl toml/param_grid10x10.toml
# julia runCycleWalk_toml.jl toml/param_hex10x10.toml
# julia runCycleWalk_toml.jl toml/param_grid4x4.toml --thread_id 200 --two_cycle_walk_frac .1 --cycle_walk_steps 1e4

import Pkg
Pkg.activate(".")
Pkg.instantiate()

using RandomNumbers
using CycleWalk

using UnPack, TOML, ArgMacros
include("runtimeParameters.jl") # see this file for parsing commandline args and passed toml data. 
                                # It also set some variables derived from the toml data or commandline args


base_graph = BaseGraph(pctGraphPath, pop_col, inc_node_data=node_data,
                       area_col=area_col, node_border_col=node_border_col, 
                       edge_perimeter_col=edge_perimeter_col)
graph = MultiLevelGraph(base_graph, geo_units)

constraints = initialize_constraints()
add_constraint!(constraints, PopulationConstraint(graph, num_dists, pop_dev))

rng = PCG.PCGStateOneseq(UInt64, rng_seed)
initial_partition = MultiLevelPartition(graph, constraints, num_dists; 
                                        rng=rng);

partition = LinkCutPartition(initial_partition, rng);

cycle_walk = build_two_tree_cycle_walk(constraints)
internal_walk = build_one_tree_cycle_walk(constraints)
proposal = [(two_cycle_walk_frac, cycle_walk), 
            (1.0-two_cycle_walk_frac, internal_walk)]

measure = Measure()
for fnct_str in measure_scores
    if fnct_str=="get_log_spanning_forests"
        push_energy!(measure, get_log_spanning_forests, gamma) 
    elseif fnct_str=="get_isoperimetric_score"
        push_energy!(measure, get_isoperimetric_score, iso_weight)
    else
        fnct = getfield(LiftedTreeWalk, Symbol(statistic))
        push_energy!(measure, fnct)
    end
end

@show output_file_path
writer = Writer(measure, constraints, partition, output_file_path; 
                additional_parameters=ad_param)
for stat in writer_stats
    fnct = getfield(LiftedTreeWalk, Symbol(stat))
    push_writer!(writer, fnct)
end

run_diagnostics = RunDiagnostics()

if params["run"]["run_diagnostics"]==true
    push_diagnostic!(run_diagnostics, cycle_walk, AcceptanceRatios(), 
                    desc = "cycle_walk")
    push_diagnostic!(run_diagnostics, cycle_walk, CycleLengthDiagnostic())
    push_diagnostic!(run_diagnostics, cycle_walk, DeltaNodesDiagnostic())
end
println("\n\nRunning CycleWalk with $steps steps, output every $outfreq steps")
println("Using $two_cycle_walk_frac fraction of two_cycle_walk proposals")
println("Using $gamma gamma and $iso_weight iso_weight")
println("Using $num_dists districts and $pop_dev population deviation")

run_metropolis_hastings!(partition, proposal, measure, steps, rng,
                         writer=writer, output_freq=outfreq,
                         run_diagnostics=run_diagnostics);
                        
close_writer(writer)
println("CycleWalk completed. Output written to $output_file_path")

