## Run:
# julia runCycleWalk_ct.jl

## Activate the CycleWalk environment  and load necessary packages
import Pkg
push!(LOAD_PATH, "..")
# Pkg.activate(".")
# Pkg.instantiate()

using RandomNumbers
using CycleWalk

twocycle_frac = 0.1
base_gamma = 0.0 # 0 is spanning forest measure, 1 is partition
target_gamma = 1.0 
base_iso_weight = 0.0 # weight on the sum of isoperimetric ratios; i.e. Polsby-Popper
target_iso_weight = 0.3
num_dists = 5
rng = PCG.PCGStateOneseq(UInt64, 4541901234)
pop_dev = 0.02 # population deviation (fraction from ideal)
cycle_walk_steps = 10^3/twocycle_frac
steps = Int(cycle_walk_steps/twocycle_frac)

total_steps = steps
base_steps_per_sample = Int(1000/twocycle_frac)
steps_per_annealing = 100# Int(10_000/twocycle_frac)

## build graph
pctGraphPath = joinpath("data","ct","CT_pct20.json")
nodeData = Set(["COUNTY", "NAME", "POP20", "area", "border_length"]);
graph = Graph(pctGraphPath, "POP20", "NAME"; inc_node_data=nodeData,
              area_col="area", node_border_col="border_length", 
              edge_perimeter_col="length")

## build partition
constraints = initialize_constraints()
add_constraint!(constraints, PopulationConstraint(graph, num_dists, pop_dev))
partition = LinkCutPartition(graph, constraints, num_dists; rng=rng, 
                             verbose=true);

## build proposal
cycle_walk = build_two_tree_cycle_walk(constraints)
internal_walk = build_one_tree_cycle_walk(constraints)
proposal = [(twocycle_frac, cycle_walk), 
            (1.0-twocycle_frac, internal_walk)]

## build measure
measure = Measure()
push_energy!(measure, get_log_spanning_forests, target_gamma) # add spanning forests energy
push_energy!(measure, get_isoperimetric_score, target_iso_weight) # add isoperimetric score energy

## establish output name and path
atlasName = "annealedImportanceSampling_2cyclefrac_"*string(twocycle_frac)
atlasName *= "_gamma"*string(target_gamma)
atlasName *= "_iso"*string(target_iso_weight)
atlasName *= ".jsonl.gz" # or just ".jsonl" for an uncompressed output
output_file_path = joinpath("output","ct", atlasName) # add output directory to path

## establish writer to which the output will be written
ad_param = Dict{String, Any}("popdev" => pop_dev) # specific info to write
ad_param["target_gamma"] = target_gamma
ad_param["target_iso_weight"] = target_iso_weight
ad_param["base_gamma"] = base_gamma
ad_param["base_iso_weight"] = base_iso_weight
writer = Writer(measure, constraints, partition, output_file_path; 
                additional_parameters=ad_param)
push_writer!(writer, get_log_spanning_trees) # add spanning trees count to writer
push_writer!(writer, get_log_spanning_forests) # add spanning forests count to writer
push_writer!(writer, get_isoperimetric_scores) # add isoperimetric scores to writer

## optional run diagnostics
# run_diagnostics = RunDiagnostics()
# push_diagnostic!(run_diagnostics, cycle_walk, AcceptanceRatios(), 
#                  desc = "cycle_walk")
# push_diagnostic!(run_diagnostics, cycle_walk, CycleLengthDiagnostic())
# push_diagnostic!(run_diagnostics, cycle_walk, DeltaNodesDiagnostic())

## run MCMC sampler
function modify_measure!(
    m::Measure,
    cur_step::Int,
    final_step::Int
)
    run_frac = cur_step/final_step
    m.weights[get_log_spanning_forests] = target_gamma*run_frac
    m.weights[get_isoperimetric_score] = target_iso_weight*run_frac
end

println("running annealed importance sampling; outputting here: "*output_file_path)
run_annealed_importance_sampling!(partition, proposal, measure, modify_measure!, 
                                  total_steps, base_steps_per_sample, 
                                  steps_per_annealing, rng; writer=writer)
close_writer(writer) # close atlas