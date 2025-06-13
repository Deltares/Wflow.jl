module LandSurfaceTemperature

using Wflow

using Wflow
using NCDatasets
using Parameters

# after Devi Purnamasari et al. 2025 
# https://doi.org/10.5194/hess-29-1483-2025
#@kw from Parameters package 
#a struct or "object" is a collection of fields
#a field is a variable that belongs to a struct
#@with_kw is a macro that creates a struct with fields and default values
#kw is short for keyword arguments part of Parameters package
@with_kw struct LSTModel
    albedo::Vector{Float64}  # Surface albedo e.g. LSA-SAF
    dssf::Vector{Float64}    # Downward shortwave radiation e.g. LSA-SAF
    ra::Vector{Float64}      # Aerodynamic resistance
    ts::Vector{Float64}      # Land surface temperature
    #constants...
    sigma::Float64 = 5.67e-8  # Stefan-Boltzmann constant for blackbody radiation
    cp::Float64 = 1005.0      # Specific heat capacity of air in J/kg/K
    rhoa::Float64 = 1.225     # Density of air in kg/m3
    # lambda as a fn of Ta: 2501 - 2.375*Ta
    # lambda::Float64 = 2.45e6  # Latent heat of vaporization in J/kg 
    hc::Float64 = 2.0        # crop height in m
    k::Float64 = 0.41        # von Kármán constant
end

function initialize_lst(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    dt::Second
)
    RS_in, albedo, crop_height = read_lst_inputs(dataset, config, indices)
    n = length(indices)
    ts = fill(0.0, n)  # Initial surface temperature
    
    return LSTModel(
        albedo = albedo,
        dssf = RS_in,
        ra = crop_height,
        ts = ts
    )
end