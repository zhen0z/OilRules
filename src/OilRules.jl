__precompile__()
module OilRules

    export OilRulesPATH, Darcy, md

    using LinearAlgebra
    using Statistics
    using JutulDarcy
    using Jutul
    using Optim
    using Flux
    using ChainRulesCore
    import Jutul: JutulGeometry, get_facepos, compute_face_trans, compute_half_face_trans, expand_perm
    import Jutul: SimulationModel, select_output_variables!
    import Jutul: optimization_targets, variable_mapper, optimization_limits, print_parameter_optimization_config
    import Jutul: objective_opt!, gradient_opt!, objective_and_gradient_opt!
    import Base: +, -, *, /, ==
    import Base: display, length, size, getindex, setindex!, IndexStyle, vec, firstindex, lastindex
    import LinearAlgebra: norm, dot
    import ChainRulesCore: rrule

    OilRulesPATH = dirname(pathof(OilRules))

    visOil = 1e-2
    visH2O = 1e-3
    ρOil = 850
    ρH2O = 1e3
    bar = 1e5

    const Darcy = 9.869232667160130e-13
    const md = Darcy * 1e-3
    
    const day = 24*3600.0

    include("PropertyConversion/PropertyConversion.jl")
    include("FlowRules/FlowRules.jl")

end # module OilRules
