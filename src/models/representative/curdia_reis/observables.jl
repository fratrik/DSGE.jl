function init_observable_mappings!(m::CurdiaReis)

    observables = OrderedDict{Symbol,Observable}()
    population_mnemonic = get(get_setting(m, :population_mnemonic))

    ############################################################################
    ## 1. Real GDP Growth
    ############################################################################
    gdp_fwd_transform = function (levels)
        # FROM: Level of GDP (from FRED)
        # TO: Quarter-to-quarter percent change of real GDP per capita

        oneqtrpctchange(levels[!, :PRS85006163])
    end

    gdp_rev_transform = loggrowthtopct_annualized_percapita

    observables[:obs_output] = Observable(:obs_output, [:PRS85006163__FRED],
                                       gdp_fwd_transform, gdp_rev_transform,
                                       "Ouput Growth", "Nonfarm Business Ouput Growth Per Capita")

    ############################################################################
    ## 2. Hours
    ############################################################################

    hours_fwd_transform = function (levels)
        # FROM: CPI urban consumers index (from FRED)
        # TO: Annualized quarter-to-quarter percent change of CPI index
        levels[!, :temp] = percapita(m, :HOANBS, levels)

        oneqtrpctchange(levels[!,:temp])
    end

    hours_rev_transform = loggrowthtopct_annualized

    observables[:obs_hours] = Observable(:obs_hours, [:HOANBS__FRED, population_mnemonic],
                                        hours_fwd_transform, hours_rev_transform,
                                        "Hours Worked",
                                        "Nonfarm Business Hours Per Capita")

    m.observable_mappings = observables
end
