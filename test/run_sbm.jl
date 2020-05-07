function update(model, toposort, n)
    # update the vertical
    for i in eachindex(model.vertical)
        model.vertical[i] = Wflow.update_before_lateralflow(model.vertical[i])
        model.lateral.recharge[i] = model.vertical[i].recharge  * model.lateral.dl[i]
        model.lateral.zi[i] = model.vertical[i].zi
    end

    Wflow.update(model.lateral, model.network, toposort, n)

    sbm = Vector{Wflow.SBM}(undef, n)
    for i in eachindex(model.vertical)
        sbm[i] = Wflow.SBM{Float64,model.vertical[i].nlayers,model.vertical[i].nlayers + 1}(
            maxlayers = model.vertical[i].maxlayers,
            nlayers = model.vertical[i].nlayers,
            riverfrac = model.vertical[i].riverfrac,
            cfmax = model.vertical[i].cfmax,
            tt = model.vertical[i].tt,
            tti = model.vertical[i].tti,
            ttm = model.vertical[i].ttm,
            whc = model.vertical[i].whc,
            w_soil = model.vertical[i].w_soil,
            cf_soil = model.vertical[i].cf_soil,
            θₛ = model.vertical[i].θₛ,
            θᵣ = model.vertical[i].θᵣ,
            kv₀ = model.vertical[i].kv₀,
            kvfrac = model.vertical[i].kvfrac,
            m = model.vertical[i].m,
            hb = model.vertical[i].hb,
            soilthickness = model.vertical[i].soilthickness,
            act_thickl = model.vertical[i].act_thickl,
            sumlayers = model.vertical[i].sumlayers,
            infiltcappath = model.vertical[i].infiltcappath,
            infiltcapsoil = model.vertical[i].infiltcapsoil,
            maxleakage = model.vertical[i].maxleakage,
            waterfrac = model.vertical[i].waterfrac,
            pathfrac = model.vertical[i].pathfrac,
            altitude = model.vertical[i].altitude,
            rootingdepth = model.vertical[i].rootingdepth,
            rootdistpar = model.vertical[i].rootdistpar,
            capscale = model.vertical[i].capscale,
            et_reftopot = model.vertical[i].et_reftopot,
            sl = model.vertical[i].sl,
            swood = model.vertical[i].swood,
            kext = model.vertical[i].kext,
            c = model.vertical[i].c,
            e_r = model.vertical[i].e_r,
            cmax = model.vertical[i].cmax,
            canopygapfraction = model.vertical[i].canopygapfraction,
            lai = model.vertical[i].lai,
            canopystorage = model.vertical[i].canopystorage,
            snow = model.vertical[i].snow,
            snowwater = model.vertical[i].snowwater,
            tsoil = model.vertical[i].tsoil,
            actinfilt = model.vertical[i].actinfilt,
            recharge = model.vertical[i].recharge,
            transpiration = model.vertical[i].transpiration,
            soilevap = model.vertical[i].soilevap,
            interception = model.vertical[i].interception,
            ae_openw_r = model.vertical[i].ae_openw_r,
            ae_openw_l = model.vertical[i].ae_openw_l,
            actevapsat = model.vertical[i].actevapsat,
            actevap = model.vertical[i].actevap,
            ustorelayerdepth = model.vertical[i].ustorelayerdepth,
            transfer = model.vertical[i].transfer,
            satwaterdepth = model.vertical[i].satwaterdepth,
            actinfiltsoil = model.vertical[i].actinfiltsoil,
            actinfiltpath = model.vertical[i].actinfiltpath,
            excesswater = model.vertical[i].excesswater,
            excesswatersoil = model.vertical[i].excesswatersoil,
            excesswaterpath = model.vertical[i].excesswaterpath,
            zi = model.lateral.zi[i],
            exfiltsatwater = model.lateral.exfiltwater[i]
        )
    end

    model = Wflow.Model(model.network, model.lateral, sbm, model.clock, nothing, nothing)

    for i in eachindex(model.vertical)
        model.vertical[i] = Wflow.update_after_lateralflow(model.vertical[i])
    end
    # update the clock
    model.clock.iteration += 1
    model.clock.time += model.clock.Δt

    return model
end

model = Wflow.initialize_sbm_model(staticmaps_moselle_path, leafarea_moselle_path)
toposort = Wflow.topological_sort_by_dfs(model.network)
n = length(toposort)
model = update(model, toposort, n)
