"""
```
CurdiaReis{T} <: AbstractRepModel{T}
```

The `CurdiaReis` type defines the structure of the simple New Keynesian DSGE
model described in 'Bayesian Estimation of DSGE Models' by Sungbae An and Frank
Schorfheide.

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model
  parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed
  as a function of elements of `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model
  parameters and steady-states to their indices in `parameters` and
  `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readable names to row and
column indices in the matrix representations of of the measurement equation and
equilibrium conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column in
  the measurement and equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in
  the measurement and equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a
  column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium
  condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states to
  their columns in the measurement and equilibrium condition equations. These
  are added after `gensys` solves the model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the
  model's measurement equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"curdia_reis\", cached
  here for filepath computation.

* `subspec::String`: The model subspecification number, indicating that some
  parameters from the original model spec (\"ss0\") are initialized
  differently. Cached here for filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation
  without changing the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as
  Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`,
  settings from `m.test_settings` are used in place of those in `m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model units.
  DSGE.jl will fetch data from the Federal Reserve Bank of
  St. Louis's FRED database; all other data must be downloaded by the
  user. See `load_data` and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct CurdiaReis{T} <: AbstractRepModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,Int}        #
    endogenous_states_augmented::OrderedDict{Symbol,Int}   #
    observables::OrderedDict{Symbol,Int}                   #
    pseudo_observables::OrderedDict{Symbol,Int}            #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::CurdiaReis) = "Julia implementation of model defined in 'Correlated Disturbances and U.S. Business Cycles' by Vasco Cúrdia and Ricardo Reis: CurdiaReis, $(m.subspec)"

"""
`init_model_indices!(m::CurdiaReis)`

Arguments:
`m:: CurdiaReis`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::CurdiaReis)
    # Endogenous states
    endogenous_states = collect([
        :w_t, :y_t, :γ_t, :c_t, :n_t, :k_t])

    # Exogenous shocks
    exogenous_shocks = collect([
        :a_sh, :g_sh])

    # Expectations shocks
    expected_shocks = collect([])

    # Equilibrium conditions
    equilibrium_conditions = collect([
        :eq_w, :eq_y, :eq_γ, :eq_n, :eq_euler, :eq_a, :eq_g])

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
    for (i,k) in enumerate(pseudo_observables);          m.pseudo_observables[k]          = i end
end


function CurdiaReis(subspec::String="ss0";
                       custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                       testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # initialize empty model
    m = CurdiaReis{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set settings
    model_settings!(m)
    default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)
    init_pseudo_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)

    return m
end

parg(n, v, μ, τ, description, tex_label) =
    parameter(n, float(v), (1e-20, 1e5), (1e-20, 1e5), ModelConstructors.Untransformed(),
        GammaAlt(float(μ), float(τ)), fixed=false, description=description, tex_label=tex_label)

partn(n, μ, σ, description, tex_label) =
    parameter(n, μ, (1e-20, 1-1e7), (1e-20, 1-1e7), ModelConstructors.Untransformed(),
        Normal(μ, σ), fixed=false, description=description, tex_label=tex_label)

function par_ig(n, ν, τ, description, tex_label)
    d = RootInverseGamma(ν, τ)
    parameter(n, mode(d), (1e-20, 1e5), (1e-20, 1e5), ModelConstructors.Untransformed(),
        d, fixed=false, description=description, tex_label=tex_label)
end

"""
```
init_parameters!(m::CurdiaReis)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::CurdiaReis)
    # Initialize parameters
    m <= parg(:γ, 2/3, 1.4, 1.0, "γ: Intertemporal substitution.", "\\gamma")

    m <= parg(:θ, 5, 5.5, 2.0, "θ: Labor elasticity.", "\\theta")

    m <= partn(:Φ_A, .75, .14, "ΦA: AR1 of A process", "\\Phi_A")

    m <= partn(:Φ_G, .4, .18, "ΦA: AR1 of G process", "\\Phi_G")

    m <= par_ig(:σ_A, 1.94, 1e-4, "σ_A: Standard deviation of shocks to technology.", "\\sigma_A")

    m <= par_ig(:σ_G, 1.53, 0.28, "σ_G: Standard deviation of shocks to technology.", "\\sigma_G")
end

"""
```
steadystate!(m::CurdiaReis)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::CurdiaReis)
    return m
end

function model_settings!(m::CurdiaReis)
    default_settings!(m)

    # Data
    m <= Setting(:data_id, 0, "Dataset identifier")
    m <= Setting(:cond_full_names, [:obs_output, :obs_hours],
        "Observables used in conditional forecasts")
    m <= Setting(:cond_semi_names, [:obs_output, :obs_hours],
        "Observables used in semiconditional forecasts")

    # Metropolis-Hastings
    m <= Setting(:mh_cc, 0.27,
                 "Jump size for Metropolis-Hastings (after initialization)")

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
end

function shock_groupings(m::CurdiaReis)
    tfp = ShockGroup("a", [:a_sh], RGB(1.0, 0.55, 0.0)) # darkorange
    gov = ShockGroup("g", [:g_sh], RGB(0.70, 0.13, 0.13)) # firebrick
    # pol = ShockGroup("pol", vcat([:rm_sh], [Symbol("rm_shl$i") for i = 1:n_anticipated_shocks(m)]),
    #                  RGB(1.0, 0.84, 0.0)) # gold
    # det = ShockGroup("dt", [:dettrend], :gray40)

    # return [tfp, gov, pol, det]
    return [tfp, gov]
end
