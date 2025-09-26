name = "small square test graph (unweighted), 4 districts, pop=8, gamma∈{0,1}"
@testset "$name" begin
    constraints = initialize_constraints()
    add_constraint!(constraints, PopulationConstraint(4,4))

    observed_cuts = get_observed_cut_edges(small_square_graph, constraints, 4,
                                           cycle_steps=100_000)

    # check that the observed districts cut counts are correct 
    @test length(values(observed_cuts)) == 4

    # test distribution (calculated explicitly)
    ce_to_count = Dict(8=>256,10=>(128+64+32),11=>96,12=>78)
    steps = sum(values(observed_cuts))
    @test all([is_close(v/steps, ce_to_count[k]/654) for (k,v) in observed_cuts])

    ###################

    measure = Measure()
    push_energy!(measure, get_log_spanning_forests, 1.0) 

    observed_cuts = get_observed_cut_edges(small_square_graph, constraints, 4,
                                           measure, cycle_steps=100_000)

    # check that the observed districts cut counts are correct 
    @test length(values(observed_cuts)) == 4

    # test distribution (calculated explicitly)
    ce_to_count = Dict(8=>1,10=>14,11=>24,12=>78)
    steps = sum(values(observed_cuts))
    @test all([is_close(v/steps, ce_to_count[k]/117) for (k,v) in observed_cuts])
end
