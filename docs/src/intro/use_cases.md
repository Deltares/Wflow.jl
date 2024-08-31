# Case studies

## [Wflow models for the Meuse and Rhine](@id case_rws)

Reliable hydrological models for the Rhine and the Meuse river basins are necessary for
short-term forecasting of river flows and long-term predictions for strategic water
management planning. In collaboration with Rijkswaterstaat, Deltares is developing a new
line of models for the Rhine and the Meuse basins. The models will be used for forecasting
and to estimate the impact of climate change on water resources and extreme streamflow. In
the model development, we aim to improve hydrological predictions by including relevant
processes in the model schematization. The modularity of the Wflow framework is ideal for
this as we can easily evaluate the combination of different vertical and lateral model
components. For example, the local inertial routing for river and overland flow enables us
to consider retention of water in the floodplains, which is likely to improve extreme
streamflow predictions.

![fig_case_rws](../images/case_rhine_meuse.png)

## [Operational flood forecasting in Australia](@id case_flifs)

In Australia, there was a need for high-resolution, fast and accurate rainfall-runoff models
to provide boundary conditions for a fast and detailed flood inundation model (SFINCS). The
domain of the flood model covers the entire North and East Coast of Australia. Although
many gauging stations are available to provide real-time information, many rivers are not
covered. For these locations, Wflow\_sbm models are used to provide this real-time
information. Additionally, these models are used to provide projections for potential future
scenarios. Using the HydroMT library, all Wflow\_sbm models were automatically built. The
high level of flexibility in spatial and temporal resolution, combined with the physics-based nature
of the concept, makes Wflow\_sbm particularly suitable for ungauged basins. Furthermore, the
model is detailed and computationally efficient enough for coupling with the fast
flood inundation model SFINCS.

![fig_case_flifs](../images/case_flifs_1.png)

The results of this proof of concept are very promising. Technically, we were able to
quickly set up the Wflow\_sbm models, couple them to the flood inundation models (SFINCS), and
run the models operationally under the Delft-FEWS platform. Model validation was conducted
for two basins by comparing the results of Wflow\_sbm against observations and the
results of calibrated URBS models. This validation demonstrated that the uncalibrated Wflow\_sbm
model results were already quite satisfactory, especially given the complex nature of these
basins, which include several small and large reservoirs. We could also showed the potential
for further calibration by adjusting the KsatHorFrac parameter.

Reference: De Kleermaeker, S., Leijnse, T., Morales, Y., Druery, C., Maguire, S.,
2022. Developing a real-time data and modelling framework for operational flood inundation
      forecasting in Australia. In Hydrology & Water Resources Symposium 2022 (HWRS 2022):
      The Past, the Present, the Future. Engineers Australia.
      https://search.informit.org/doi/10.3316/informit.916755150845355

![fig_case_flifs](../images/case_flifs_2.png)

## [Simulating plastic transport in Thailand](@id case_mfa)

For the Pollution Control Board of the Government of Thailand and the World Bank, we
supported a material flow analysis of plastics in Thailand using wflow. Plastic pollution
is a growing global issue. Plastic waste enters rivers and is transported to the ocean
where it persists and threatens the health of the ocean, seas and coasts. The initial
movement of plastic waste is in many cases triggered by runoff from (heavy) rainfall and
transported by water flow towards small streams and rivers. Therefore there is strong
relation to rainfall-runoff processes, which can be modeled using high-resolution
rainfall-runoff models.

In this study we applied the Wflow\_sbm model in combination with a fate-and-transport and
water quality model (DelWaq) to simulate the movement of plastics through five large river
basins and on three island and coastal zones (Krabi, Phuket, and Ko Samui; see screenshot of the
model below) in Thailand. Together with our partners Panya Consultants and HII, we were able
to identify hotspots of plastic pollution, estimate how much plastic waste would end up in the Gulf of
Thailand and recommend priority areas for reducing plastic waste reaching the sea.

![fig_case_mfa](../images/case_mfa_1.png)

The Wflow\_sbm models for the five large basins were calibrated. The presence of large dams and
reservoirs complicated calibration, but with the input for the dam operation,
the model performance for these basins could be largely improved. The figure below shows the
calibrated model results for the Chao Phraya, just upstream of Bangkok. The input from the
hydrological Wflow\_sbm model was used as input for the fate and transport model to assess
the amount of plastic transported to the ocean.

![fig_case_mfa](../images/case_mfa_3.png)

Link to World Bank report:
[https://www.worldbank.org/en/country/thailand/publication/plastic-waste-material-flow-analysis-for-thailand](https://www.worldbank.org/en/country/thailand/publication/plastic-waste-material-flow-analysis-for-thailand)
