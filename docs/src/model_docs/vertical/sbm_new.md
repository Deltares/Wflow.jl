# wflow\_sbm (@id wflow_sbm_desc)

## Introduction

The soil part of wflow\_sbm model concept is largely based on the Topog\_sbm mode, but 
has had considerable changes over time. The main changes in wflow\_sbm are:

- The unsaturated zone can be split up in different layers
- The addition of evapotranspiration losses
- The addition of capillary rise
- Wflow\_sbm routes water over a D8 network instead of an element network based on contour lines and trajectories

For the lateral components of this wflow\_sbm model water is routed over a D8 network, and
the kinematic wave approach is used for river, overland and lateral subsurface flow. This is
described in more detail in the section [Kinematic wave].

An overview of the different processes and fluxes in the wflow_sbm model:

![wflow_sbm model](../../images/wflow_sbm_soil.png)

## Precipitation 
The division between solid and liquid precipitation (snowfall and rainfall, respectively) is 
performed based on the air temperature. If the temperature is below a threshold temperature 
(`tt`), precipitation will fall as snow. A interval parameter (`tti`) defines the range 
over which precipitation is partly falling as snow, and partly as rain. Snowfall is added
to the snowpack, where it is subject to melting and refreezing (see the section on 
[snow and glaciers](@ref snow)). The amount of rainfall is subject to 
[interception](@ref interception), and ultimately becomes available for [evaporation](@ref evap) 
and/or [soil processes](@ref soil).

![snowfall](../../images/snowfall.png) 

## [Rainfall interception](@id interception)

Two different interception models are available: the analytical Gash model, and the modified
Rutter model. The simulation timestep defines which interception model is used, where daily
(or larger) timesteps use the Gash model, and timesteps smaller than daily use the modified 
Rutter model.

### The analytical Gash model (Gash, 1979)
The Gash model generally assumes that there is one rainfall event per day. It calculates the
amount of water required to fully saturate the canopy. The amount of stemflow is taken as 
a fraction (`0.1 * canopygapfraction`) of the precipitation. Throughfall is the precipitation 
that is not intercepted by the leaves and stems. If the interception is larger than the 
potential evaporation, the amount above potential evaporation is added to the throughfall, 
as the concept is based on a daily timestep (assuming that the canopy storage is emptied at 
the end of the timestep). 

### The modified Rutter model (Rutter et al., 1971)

The amount of stemflow is taken as a fraction (`0.1 * canopygapfraction`) of the 
precipitation. Throughfall equals to the amount of water that cannot be stored by the 
canopy, plus the rainfall that is not captured by the canopy. There is no flow from water
on the leaves to the ground, when the canopy storage is below the maximum canopy storage 
(`cmax`). Water can evaporate from the canopy storage, taken as the minimum between potential 
evaporation and the current storage. The "left-over" potential evaporation (if any) is 
returned as output.

## [Evaporation](@id evap)

The wflow\_sbm model assumes the input to be potential evaporation. A multiplication factor (`et_reftopot`, set to 1 by default) is present to correct the input evaporation if required. 

The potential evaporation left over after interception and open water evaporation (rivers and water bodies) is split in potential soil evaporation and potential transpiration based on the canopy gap fraction (assumed to be identical to the amount of bare soil). 

## [Snow and glaciers](@id snow)

The snow and glacier model is described in [Snow and glaciers](@ref snow_and_glac). Both
options can be enabled by specifying the following in the TOML file:

```toml
[model]
snow = true
glacier = true
```

## [Soil processes](@id soil)

### The SBM soil water accounting scheme
A detailed description of the Topog\_SBM model has been given by Vertessy (1999). Briefly:
the soil is considered as a bucket with a certain depth (``z_{t}`` [mm]), divided into a
saturated store (``S`` [mm]) and an unsaturated store (``U`` [mm]). The top of the ``S``
store forms a pseudo-water table at depth ``z_{i}`` [mm] such that the value of ``S`` at any
time is given by:
```math 
    S=(z_{t}-z_{i})(\theta_{s}-\theta_{r})
```
where ``\theta_{s}`` [-] and ``\theta_{r}`` [-] are the saturated and residual soil water
contents, respectively. 

The unsaturated store ``U`` is subdivided into storage (``U_{s}`` [mm]) and deficit
(``U_{d}`` [mm]):
```math
    U_{d}=(\theta_{s}-\theta_{r})z_{i}-U\\
    U_{s}=U-U_{d}
```
The saturation deficit (``S_{d}`` [mm]) for the soil profile as a whole is defined as: 

```math
    S_{d}=(\theta_{s}-\theta_{r})z_{t}-S
```

All infiltrating water that enters the ``U`` store first. The unsaturated layer can be
split-up in different layers, by providing the thickness [mm] of the layers in the TOML
file. The following example specifies three layers (from top to bottom) of 100, 300 and 800
mm:

```toml
[model]
thicknesslayers = [100, 300, 800]
```
The code checks for each grid cell the specified layers against the `soilthickness` [mm],
and adds or removes (partly) layer(s) based on the `soilthickness`.

Assuming a unit head gradient, the transfer of water (``st`` [mm t``^{-1}``]) from a ``U``
[mm] store layer is controlled by the saturated hydraulic conductivity ``K_{sat}`` [mm
t``^{-1}``] at depth ``z`` \[mm\] (bottom layer) or ``z_{i}`` [mm], the effective saturation
degree of the layer, and a Brooks-Corey power coefficient (parameter ``c``) based on the
pore size distribution index ``\lambda`` (Brooks and Corey, 1964):

```math
    st=K_{\mathit{sat}}\left(\frac{\theta-\theta_{r}}{\theta_{s}-\theta_{r}}\right)^{c}\\~\\
    c=\frac{2+3\lambda}{\lambda}
```
When the unsaturated layer is not split-up into different layers, it is possible to use the
original Topog\_SBM vertical transfer formulation, by specifying in the TOML file:

```toml
[model]
transfermethod = true
```
The transfer of water from the ``U`` [mm] store to the ``S`` [mm] store (``st`` [mm
t``^{-1}``]) is in that case controlled by the saturated hydraulic conductivity ``K_{sat}``
[mm t``^{-1}``] at depth ``z_{i}`` [mm] and the ratio between ``U`` [mm] and ``S_{d}``
[mm]: 

```math
    st=K_{\mathit{sat}}\frac{U_{s}}{S_{d}}
```
Saturated conductivity (``K_{sat}`` [mm t``^{-1}``]) declines with soil depth (``z`` [mm])
in the model according to: 

```math
    K_{sat}=K_{0}e^{(-fz)}
```
where ``K_{0}`` [mm t``^{-1}``] is the saturated conductivity at the soil surface and ``f``
is a scaling parameter [mm``^{-1}``].

The plot below shows the relation between soil depth ``z`` and saturated hydraulic
conductivity ``K_{sat}`` for different values of ``f``.

```@setup plot
    using Printf
    using CairoMakie
```

```@example plot
    let                                                                                     # hide
        fig = Figure(resolution = (800, 400))                                               # hide
        ax = Axis(fig[1, 1], xlabel = "Kₛₐₜ [mm/day]", ylabel = "-z [mm]")                  # hide

        z = 0:5.0:1000                                                                      # hide
        ksat = 100.0                                                                        # hide
        f = 0.6 ./ collect(50:150.0:800)                                                    # hide

        for fi in f                                                                         # hide
            lines!(ax, ksat .* exp.(-fi .* z), -z, label = @sprintf("f = %.2e", fi))        # hide
        end                                                                                 # hide

        Legend(fig[1, 2], ax, "f")                                                          # hide
        fig                                                                                 # hide
    end                                                                                     # hide
```

### Bare soil evaporation

If there is only one soil layer present in the wflow\_sbm model, the bare soil evaporation is scaled according to the wetness of the soil layer. The fraction of bare soil is assumed to be equal to the fraction not covered by the canopy (`conapygapfraction`). When the soil is fully saturated, evaporation is set to equal the potential evaporation. When the soil is not fully saturated, actual evaporation decrease linearly with decreasing soil moisture values, as indicated by the figure below.

![soil_evap](../../images/soil_evap.png) 

When more soil layers are present, soil evaporation is only provided from the upper soil layer, and soil evaporation is split in evaporation from the unsaturated store and evaporation from the saturated store. Water is first evaporated water from the unsaturated store. The remaining potential soil evaporation can be used for evaporation from the saturated store, but only when the water table is present in the upper soil layer. 

### Transpiration

The fraction of wet roots is determined using a sigmoid fuction (see figure below). The parameter `rootdistpar` defines the sharpness of the transition between fully wet and fully dry roots. The returned wetroots fraction is multiplied by the potential evaporation (and limited by the available water in saturated zone) to get the transpiration from the saturated part of the soil. This is implemented using the following code (`i` refers to the index of the vector that contains all active cells within the spatial model domain):
```julia
    # transpiration from saturated store
    wetroots = scurve(sbm.zi[i], rootingdepth, 1.0, sbm.rootdistpar[i])
    actevapsat = min(pottrans * wetroots, satwaterdepth)
    satwaterdepth = satwaterdepth - actevapsat
    restpottrans = pottrans - actevapsat
```

![soil_wetroots](../../images/soil_wetroots.png) 

The remaining potential evaporation is used to extract water from the unsaturated store. The maximum allowed extration of the unsaturated zone is determined based on the fraction of the unsaturated zone that is above the rooting depth, see conceptual figure below. This is implemented using the following code:
```julia
    if ust # whole_ust_available = true
        availcap = ustorelayerdepth * 0.99
    else
        if usl > 0.0
            availcap = min(1.0, max(0.0, (rootingdepth - sumlayer) / usl))
        else
            availcap = 0.0
        end
    end
    maxextr = availcap * ustorelayerdepth 
```

![soil_unsatevap](../../images/soil_unsatevap.png) 

!!! note
    When `whole_ust_available` is set to true in the TOML file, almost the complete unsaturated storage (99%) is available for transpiration, independent of the `rootingdepth`.

    ```toml
    [model]
    whole_ust_available = true
    ```

Next, a root water uptake reduction model is used to calculate a reduction coefficient as a function of soil water pressure. This concept is based on the concept presented by Feddes et al. (1978). This concept defines a reduction coefficient `a` as a function of soil water pressure (`h`). Four different levels of `h` are defined: `h2`, `h3`, and `h4` are defined as fixed values, and `h1` can be defined as input to the model (defaults to -10 cm). `h1` represents the air entry pressure, `h2` represents field capacity, `h3` represents the point of critical soil moisture content, and `h4` represents the wilting point. The current soil water pressure is determined following the concept defined by Brooks and Corey (1964): 

```math
    \frac{(\theta-\theta_r)}{(\theta_s-\theta_r)} =  \Bigg\lbrace{\left(\frac{h_b}{h}\right)^{\lambda}, h > h_b \atop 1 , h \leq h_b}
```
where ``h`` is the pressure head [cm], ``h_b`` is the air entry pressure head [cm], and
``\theta``, ``\theta_s``, ``\theta_r`` and ``\lambda`` as previously defined.

Whenever the current soil water pressure drops below `h4`, the root water uptake is set to zero. The root water uptake is at ideal conditions whenever the soil water pressure is above `h3`, with a linear transition between `h3` and `h4`. In the original concept, root water uptake is set be reduced when soil water pressures are above field capacity (`h2`). This is not inplemented in wlow\_sbm, as the assumption from the original concept does not apply to crops. 

![soil_rootwateruptake](../../images/soil_rootwateruptake.png) 

### Infiltration

The water available for infiltration is taken as the rainfall including meltwater. Infiltration is determined seperately for the compacted and non compacted areas, as these have different infitration capacities. Naturally, only the water that can be stored in the soil can infiltrate. If not all water can infiltrate, this is added as excess water to the runoff routing scheme. 

The infiltrating
water is split in two parts, the part that falls on compacted areas and the part that falls
on non-compacted areas. The maximum amount of water that can infiltrate in these areas is
calculated by taking the minimum of the maximum infiltration rate (`infiltcapsoil` [mm
t``^{-1}``] for non-compacted areas and `infiltcappath` [mm t``^{-1}``] for compacted areas)
and the amount of water available for infiltration `avail_forinfilt` [mm t``^{-1}``]. The
water that can actual infiltrate `infiltsoilpath` [mm t``^{-1}``] is calculated by taking
the minimum of the total maximum infiltration rate (compacted and non-compacted areas) and
the remaining storage capacity.

Infiltration excess occurs when the infiltration capacity is smaller then the throughfall
and stemflow rate. This amount of water (`infiltexcess` [mm t``^{-1}``]) becomes overland
flow (infiltration excess overland flow). Saturation excess occurs when the (upper) soil
becomes saturated and water cannot infiltrate anymore. This amount of water `excesswater`
[mm t``^{-1}``] becomes overland flow (saturation excess overland flow).

If snow processes are modelled, the infiltration capacity is reduced when the soil is frozen (or near freezing point). A infiltration correction factor is defined as a S-curve with the shape as defined below. A parameter (`cf_soil`) defines the base factor of infiltration when the soil is frozen. The soil temperature is calculated based on the soil temperature on the previous timestep, and the temperature difference between air and soil temperature weighted with a factor (`w_soil`, which defaults to 0.1125).

The near surface soil temperature is modelled using a simple equation (Wigmosta et al.,
2009):

```math
T_s^{t} = T_s^{t-1} + w  (T_a - T_s^{t-1})  
```
where ``T_s^{t}`` [``\degree``C] is the near-surface soil temperature at time ``t``, ``T_a``
[``\degree``C] is air temperature and ``w`` [-] is a weighting coefficient determined
through calibration (default is 0.1125 for daily timesteps).

A reduction factor (`cf_soil` [-], default is 0.038) is applied to the maximum infiltration
rate (`infiltcapsoil` and `infiltcappath`), when the following model settings are specified
in the TOML file:

```toml
[model]
soilinfreduction = true
snow = true
```
If `soilinfreduction` is set to `false`, water is allowed to infiltrate the soil, even if the soil is frozen.

A S-curve (see plot below) is used to make a smooth transition (a c-factor (``c``) of 8.0 is
used):

```math 
    b = \frac{1.0}{(1.0 - cf\_soil)}\\~\\
    soilinfredu = \frac{1.0}{b + exp(-c (T_s - a))} + cf\_soil\\~\\
    a = 0.0\\
    c = 8.0
```

![soil_frozeninfilt](../../images/soil_frozeninfilt.png) 

### Capillary rise

The actual capillary rise `actcapflux` [mm t``^{-1}``] is determined using the following
approach: first the saturated hydraulic conductivty `ksat` [mm t``^{-1}``] is determined at
the water table ``z_{i}``; next a potential capillary rise `maxcapflux` [mm t``^{-1}``] is
determined from the minimum of `ksat`, actual transpiration `actevapustore` [mm t``^{-1}``]
taken from the ``U`` store, available water in the ``S`` store (`satwaterdepth` [mm]) and
the deficit of the ``U`` store (`ustorecapacity` [mm]), as shown by the following code
block:

```julia
    maxcapflux = max(0.0, min(ksat, actevapustore, ustorecapacity, satwaterdepth))
```

Then the potential rise `maxcapflux` is scaled using the water table depth `zi` and a
maximum water depth `cap_hmax` [mm] beyond which capillary rise ceases as follows in the
code block below (`i` refers to the index of the vector that contains all active cells
within the spatial model domain):

```julia
    if sbm.zi[i] > rootingdepth
        capflux =
            maxcapflux *
            pow(1.0 - min(sbm.zi[i], sbm.cap_hmax[i]) / (sbm.cap_hmax[i]), 2.0)
    else
        capflux = 0.0
    end
```

If the roots reach the water table (`rootingdepth` ``\ge`` `sbm.zi`), `capflux` is set to
zero and thus setting the capillary rise `capflux` to zero.

Finally, the capillary rise `capflux` is limited by the unsaturated store deficit (one or
multiple layers), calculated as follows in the code block below (`i` refers to the index of
the vector that contains all active cells within the spatial model domain, and `k` refers to
the layer position):

```julia
    usl[k] * (sbm.θₛ[i] - sbm.θᵣ[i]) - usld[k]
```

where `usl` [mm] is the unsaturated layer thickness, `usld` is the `ustorelayerdepth` \[mm\]
(amount of water in the unsaturated layer), and ``\theta_{s}`` and ``\theta_{r}`` as
previously defined.

The calculation of the actual capillary rise `actcapflux` is as follows in the code block
below (`i` refers to the index of the vector that contains all active cells within the
spatial model domain, and `k` refers to the layer position):

```julia
    actcapflux = 0.0
    netcapflux = capflux
    for k = n_usl:-1:1
        toadd =
            min(netcapflux, max(usl[k] * (sbm.θₛ[i] - sbm.θᵣ[i]) - usld[k], 0.0))
        usld = setindex(usld, usld[k] + toadd, k)
        netcapflux = netcapflux - toadd
        actcapflux = actcapflux + toadd
    end
```

In case of multiple unsaturated layers (`n_usl` ``>`` 1), the calculation of the actual
capillary rise starts at the lowest unsaturated layer while keeping track of the remaining
capillary rise `netcapflux` [mm t``^{-1}``].

### Leakage

If the `maxleakage` (mm/day) input model parameter is set > 0, water is lost from the saturated zone and runs out of the model.

## Open water

Part of the water available for infiltration is diverted to the open water, based on the fractions of river and lakes of each grid cell. The amount of evaporation from open water is taken assumed to be equal to potential evaporation (if sufficient water is available).