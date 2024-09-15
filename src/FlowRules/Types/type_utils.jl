function force(M::jutulModel{D, T}, w::jutulForce{D, T}, tstep::Vector{T};
    ρOil::T=T(ρOil), ρH2O::T=T(ρH2O), g::T=T(10.0)) where {D, T}

    ## set up well information
    cell_loc = [Int.(round.(w.loc[i] ./ M.d[1:length(w.loc[1])])) for i = 1:length(w.loc)]
    Is = [setup_well(CartesianMesh(M), M.K, [cell_loc[i]], name = w.name[i]) for i = 1:length(w.loc)]
    ctrls = [w.name[i]==:Injector ? InjectorControl(TotalRateTarget(w.irate), [1.0, 0.0], density = ρOil) : ProducerControl(BottomHolePressureTarget(50*bar)) for i = 1:length(w.loc)]
    controls = Dict()
    for i = 1:length(w.loc)
        controls[w.name[i]] = ctrls[i]
    end
    return Is, controls
end

function force(M::jutulModel{D, T}, w::jutulVWell{D, T}, tstep::Vector{T};
    ρOil::T=T(ρOil), ρH2O::T=T(ρH2O), g::T=T(10.0)) where {D, T}

    ## set up well information
    cell_loc = [Int.(round.(w.loc[i] ./ M.d[1:length(w.loc[1])])) for i = 1:length(w.loc)]
    heel = [isnothing(w.startz) ? 1 : Int(div(w.startz[i], M.d[end])) for i = 1:length(w.loc)]
    toe = [isnothing(w.endz) ? M.n[end] : Int(div(w.endz[i], M.d[end])) for i = 1:length(w.loc)]
    Is = [setup_vertical_well(CartesianMesh(M), M.K, cell_loc[i]..., name = w.name[i]; heel = heel[i], toe = toe[i]) for i = 1:length(w.loc)]
    ctrls = [w.name[i]==:Injector ? InjectorControl(TotalRateTarget(w.irate), [1.0, 0.0], density = ρOil) : ProducerControl(BottomHolePressureTarget(50*bar)) for i = 1:length(w.loc)]
    controls = Dict()
    for i = 1:length(w.loc)
        controls[w.name[i]] = ctrls[i]
    end
    return Is, controls
end

function setup_well_model(M::jutulModel{D, T}, f::Union{jutulForce{D, T}, jutulVWell{D, T}}, tstep::Vector{T};
    visOil::T=T(visOil), visH2O::T=T(visH2O), ρOil::T=T(ρOil), ρH2O::T=T(ρH2O), g::T=T(10.0)) where {D, T}

    ### set up well controls
    Is, controls = force(M, f, tstep; ρOil=ρOil, ρH2O=ρH2O, g=g)    

    ### set up model, parameters
    sys = ImmiscibleSystem((VaporPhase(), AqueousPhase()), reference_densities = [ρOil, ρH2O])
    domain_spec = reservoir_domain(CartesianMesh(M), porosity = M.ϕ, permeability = M.K)
    domain = discretized_domain_tpfv_flow(domain_spec)
    model_parameters = Dict(:Reservoir => Dict(:PhaseViscosities=> [visOil, visH2O]))
    model, parameters = setup_reservoir_model(domain_spec, sys, wells = Is, parameters=model_parameters)
    select_output_variables!(model.models.Reservoir, :all)
    ρ = ConstantCompressibilityDensities(p_ref = 150*bar, density_ref = [ρOil, ρH2O], compressibility = [1e-4/bar, 1e-6/bar])
    replace_variables!(model, PhaseMassDensities = ρ)
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelPerm(sys, [2.0, 2.0], [0.1, 0.1], 1.0))
    for x ∈ keys(model.models)
        Jutul.select_output_variables!(model.models[x], :all)
    end

    ### forces
    forces = setup_reservoir_forces(model, control = controls)

    ### initial state
    Z = repeat((1:M.n[end])*M.d[end], inner = prod(M.n[1:2]))
    p0 = ρH2O * g * (Z .+ M.h) # rho * g * h
    state0 = setup_reservoir_state(model, Pressure = p0, Saturations = [0.0, 1.0])

    return model, parameters, state0, forces
end

function source(M::jutulModel{D, T}, f::jutulSource{D, T}; ρOil::T=T(ρOil)) where {D, T}
    model = simple_model(M; ρOil=ρOil)
    cell_loc = [Int.(round.(f.loc[i] ./ M.d)) for i = 1:length(f.loc)]
    cell = [sum([(cell_loc[i][d]-1) * prod(M.n[1:d-1]) for d = length(cell_loc[i]):-1:1]) + 1 for i = 1:length(cell_loc)]
    src  = [SourceTerm(cell[i], f.irate[i] * ρOil, fractional_flow = [T(f.irate[i] > 0), T(1)-T(f.irate[i] > 0)]) for i = 1:length(f.loc)]
    return setup_forces(model, sources = src)
end

function simple_model(M::jutulModel{D, T}; ρOil::T=T(ρOil), ρH2O::T=T(ρH2O)) where {D, T}
    sys = ImmiscibleSystem((VaporPhase(), AqueousPhase()))
    g = CartesianMesh(M.n, M.d .* M.n)
    domain_spec = reservoir_domain(g, porosity = M.ϕ, permeability = M.K)
    G = discretized_domain_tpfv_flow(domain_spec)
    model = SimulationModel(domain_spec, sys, output_level = :all)
    model.primary_variables[:Pressure] = JutulDarcy.Pressure(minimum = -Inf, max_rel = nothing)
    ρ = ConstantCompressibilityDensities(p_ref = 150*bar, density_ref = [ρOil, ρH2O], compressibility = [1e-4/bar, 1e-6/bar])
    replace_variables!(model, PhaseMassDensities = ρ)
    replace_variables!(model, RelativePermeabilities = BrooksCoreyRelPerm(sys, [2.0, 2.0], [0.1, 0.1], 1.0))
    return model
end
