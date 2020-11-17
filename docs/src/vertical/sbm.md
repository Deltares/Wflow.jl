# SBM vertical concept

The SBM vertical concept has its roots in the Topog\_SBM model but has had
considerable changes over time. The main differences are:

- The unsaturated zone can be split-up in different layers
- The addition of evapotranspiration losses
- The addition of a capillary rise

The sections below describe the working of the SBM vertical concept in more
detail.

## Snow
Snow modelling is enabled by specifying the following in the TOML file:

```
[model]
snow = true
```

The snow model is described in [Snow modelling](@ref)

## Soil
### Infiltration
If the surface is (partly) saturated the throughfall and stemflow that falls
onto the saturated area is added to the river runoff component (based on
fraction rivers, `riverfrac`) and to the overland runoff component (based on
open water fraction (`waterfrac`) minus `riverfrac`). Infiltration of the
remaining water is determined as follows:

The soil infiltration capacity can be adjusted in case the soil is frozen, this
is optional and can be set in the TOML file as follows:

```
[model]
# Enable iterations of the kinematic wave
soilinfreduction = true 
```
The remaining storage capacity of the unsaturated store is determined.  The
infiltrating water is split in two parts, the part that falls on compacted areas
and the part that falls on non-compacted areas. The maximum amount of water that
can infiltrate in these areas is calculated by taking the minimum of the maximum
infiltration rate (`infiltcapsoil` for non-compacted areas and `infiltcappath`
for compacted areas) and the water on these areas. The water that can actual
infiltrate is calculated by taking the minimum of the total maximum infiltration
rate (compacted and non-compacted areas) and the remaining storage capacity.

Infiltration excess occurs when the infiltration capacity is smaller then the
throughfall and stemflow rate. This amount of water (`infiltexcess`) becomes
overland flow (infiltration excess overland flow). Saturation excess occurs when
the (upper) soil becomes saturated and water cannot infiltrate anymore. This
amount of water (`excesswater` and `exfiltwater`) becomes overland flow
(saturation excess overland flow).

### The SBM soil water accounting scheme
A detailed description of the Topog\_SBM model has been given by Vertessy
(1999). Briefly: the soil is considered as a bucket with a certain depth
(``z_{t}``), divided into a saturated store (``S``) and an unsaturated store
(``U``), the magnitudes of which are expressed in units of depth. The top of the
``S`` store forms a pseudo-water table at depth ``z_{i}`` such that the value of
``S`` at any time is given by:
```math 
    S=(z_{t}-z_{i})(\theta_{s}-\theta_{r})
```
where ``\theta_{s}`` and ``\theta_{r}`` are the saturated and residual soil
water contents, respectively. 

The unsaturated store (``U``) is subdivided into storage (``U_{s}``) and deficit
(``U_{d}``) which are again expressed in units of depth:
```math
    U_{d}=(\theta_{s}-\theta_{r})z_{i}-U\\
    U_{s}=U-U_{d}
```
The saturation deficit (``S_{d}``) for the soil profile as a whole is defined
as: 

```math
    S_{d}=(\theta_{s}-\theta_{r})z_{t}-S
```

All infiltrating water that enters the ``U`` store first. The unsaturated layer
can be split-up in different layers, by providing the thickness [mm] of the
layers in the TOML file. The following example specifies three layers (from top
to bottom) of 100, 300 and 800 mm: 

```
[model]
thicknesslayers = [100, 300, 800]
```
The code checks for each grid cell the specified layers against the
`soilthickness`, and adds or removes (partly) layer(s) based on the
`soilthickness`.

Assuming a unit head gradient, the transfer of water (``st``) from a ``U`` store
layer is controlled by the saturated hydraulic conductivity ``K_{sat}`` at depth
``z`` (bottom layer) or :math:`z_{i}`, the effective saturation degree of the
layer, and a Brooks-Corey power coefficient (parameter ``c``) based on the pore
size distribution index ``\lambda`` (Brooks and Corey (1964)):

```math
    st=K_{\mathit{sat}}\left(\frac{\theta-\theta_{r}}{\theta_{s}-\theta_{r}}\right)^{c}\\~\\
    c=\frac{2+3\lambda}{\lambda}
```
When the unsaturated layer is not split-up into different layers, it is possible
to use the original Topog\_SBM vertical transfer formulation, by specifying in
the TOML file:

```
 [model]
    transfermethod = true
```
The transfer of water from the ``U`` store to the ``S`` store (``st``) is in
that case controlled by the saturated hydraulic conductivity ``K_{sat}`` at
depth ``z_{i}`` and the ratio between ``U`` and ``S_{d}``: 

```math
    st=K_{\mathit{sat}}\frac{U_{s}}{S_{d}}
```
Saturated conductivity (``K_{sat}``) declines with soil depth (``z``) in the
model according to: 

```math
    K_{sat}=K_{0}e^{(-fz)}
```
where ``K_{0}`` is the saturated conductivity at the soil surface and ``f`` is a
scaling parameter [mm``^{-1}``]. The scaling parameter :math:`f` is defined by:

```math
f=\frac{\theta_{s}-\theta_{r}}{M} 
```
with ``\theta_{s}`` and ``\theta_{r}`` as defined previously and ``M`` [mm]
representing a model parameter.

## Transpiration and soil evaporation
The potential eveporation left over after interception and open water
evaporation (rivers and water bodies) is split in potential soil evaporation and
potential transpiration based on the canopy gap fraction (assumed to be
identical to the amount of bare soil).

For the case of one single soil layer, soil evaporation is scaled according to:

    soilevapunsat = potsoilevap * min(1.0, saturationdeficit / sbm.soilwatercapacity[i])
 

As such, evaporation will be potential if the soil is fully wetted and it
decreases linear with increasing soil moisture deficit.

For more than one soil layer, soil evaporation is only provided from the upper
soil layer (often 100 mm) and soil evaporation is split in evaporation from the
unsaturated store and evaporation from the saturated store. First water is
evaporated water from the unsaturated store. Then the remaining potential soil
evaporation can be used for evaporation from the saturated store. This is only
possible, when the water table is present in the upper soil layer (very wet
conditions). Both the evaporation from the unsaturated store and the evaporation
from the saturated store are limited by the minimum of the remaining potential
soil evaporation and the available water in the unsaturated/saturated zone of
the upper soil layer. Also for multiple soil layers, the evaporation (both
unsaturated and saturated) decreases linearly with decreasing water
availability.

The original Topog\_SBM model does not include transpiration or a notion of
capillary rise. In SBM transpiration is first taken from the ``S`` store if the
roots reach the water table ``z_{i}``. If the ``S`` store cannot satisfy the
demand the ``U`` store is used next. First the number of wet roots is determined
(going from 1 to 0) using a sigmoid function as follows:

```math
    wetroots = 1.0/(1.0 + e^{-rootdistpar (zi - rootindepth)})
```

Here the sharpness parameter `rootdistpar` (by default a large negative value,
-500.0) determines if there is a stepwise output or a more gradual output
(default is stepwise). `zi` [mm] is the level of the water table in the grid
cell below the surface, `rootindepth` [mm] is the maximum depth of the roots
below the surface. For all values of `zi` smaller that `rootindepth` a value of
1 is returned if they are equal a value of 0.5 is returned if `zi` is larger
than the `rootindepth` a value of 0 is returned. The returned `wetroots`
fraction is multiplied by the potential evaporation (and limited by the
available water in saturated zone) to get the transpiration from the saturated
part of the soil:

        # transpiration from saturated store
        wetroots = scurve(sbm.zi[i], a = rootingdepth, c = sbm.rootdistpar[i])
        actevapsat = min(pottrans * wetroots, satwaterdepth)
        satwaterdepth = satwaterdepth - actevapsat
        restpottrans = pottrans - actevapsat

Next the remaining potential evaporation is used to extract water from the
unsaturated store. The fraction of roots (`availcap`) that cover the unsaturated
zone for each soil layer is used to calculate the potential root water
extraction rate (`maxextr`):

    maxextr = availcap * ustorelayerdepth

When `whole_ust_available` is set to true in the TOML file as follows, the
complete unsaturated storage is available for transpiration:

```
[model]
whole_ust_available = true
```

Next, the Feddes root water uptake reduction model (Feddes et al. (1978)) is
used to calculate a reduction coefficient as a function of soil water pressure.
Soil water pressure is calculated following Brooks and Corey (1964):

```math
    \frac{(\theta-\theta_r)}{(\theta_s-\theta_r)} =  \Bigg\lbrace{\left(\frac{h_b}{h}\right)^{\lambda}, h > h_b \atop 1 , h \leq h_b}
```
where ``h`` is the pressure head (cm), ``h_b`` is the air entry pressure head,
and ``\theta``, ``\theta_s``, ``\theta_r`` and ``\lambda`` as previously
defined.

Feddes (1978) described a transpiration reduction-curve for the reduction
coefficient ``\alpha``, as a function of ``h``.

Below, the function used in SBM, that calculates actual transpiration from the
unsaturated zone layer(s).

```@docs
Wflow.acttransp_unsat_sbm(rootingdepth, ustorelayerdepth, sumlayer, restpotevap, sum_actevapustore, c, usl, θₛ, θᵣ, hb, ust::Bool = false)
```

Capillary rise is determined using the following approach: first ``K_{sat}`` is
determined at the water table ``z_{i}``; next a potential capillary rise is
determined from the minimum of the ``K_{sat}``, the actual transpiration taken
from the ``U`` store, the available water in the ``S`` store and the deficit of
the ``U`` store. Finally the potential rise is scaled using the distance between
the roots and the water table using:

```math
CSF=CS/(CS+z_{i}-RT)
```
in which ``CSF`` is the scaling factor to multiply the potential rise with,
``CS`` is a model parameter (default = 100) and ``RT`` the rooting depth. If the
roots reach the water table (``RT > z_{i}``) ``CS`` is set to zero thus setting
the capillary rise to zero.

## Leakage

If the `maxleakage` parameter is set > 0, water is lost from the saturated zone
and runs out of the model.

## Soil temperature
The near surface soil temperature is modelled using a simple equation (Wigmosta
et al., 2009):

```math
T_s^{t} = T_s^{t-1} + w  (T_a - T_s^{t-1})  
```
where ``T_s^{t}`` is the near-surface soil temperature at time ``t``, ``T_a`` is
air temperature and ``w`` is a weighting coefficient determined through
calibration (default is 0.1125 for daily timesteps).

A reduction factor (`cf_soil`, default is 0.038) is applied to the maximum
infiltration rate (`infiltcapsoil` and `infiltcappath`), when the following
model settings are specified in the TOML file:

```
[model]
soilinfreduction = true
snow = true
```
A S-curve is used to make a smooth transition (a c-factor (``c``) of 8.0 is
used):

```math 
    b = \frac{1.0}{(1.0 - cf\_soil)}\\~\\
    soilinfredu = \frac{1.0}{b + exp(-c (T_s - a))} + cf\_soil\\~\\
    a = 0.0\\
    c = 8.0
```
