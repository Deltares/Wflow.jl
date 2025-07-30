using PkgToSoftwareBOM
using UUIDs
using Pkg

Pkg.activate("Wflow")

wflowLisence = SpdxLicenseExpressionV2("MIT")
organization = SpdxCreatorV2("Organization", "Deltares", "software@deltares.nl")
packageInstructions = spdxPackageInstructions(;
    spdxfile_toexclude = ["Wflow.spdx.json"],
    originator = organization,
    declaredLicense = wflowLisence,
    copyright = "Copyright (c) 2025 Deltares <software@deltares.nl>",
    name = "Wflow",
)

dependencies = Pkg.project().dependencies;
spdxData = spdxCreationData(;
    Name = "Wflow.jl",
    Creators = [organization],
    NamespaceURL = "https://github.com/Deltares/Wflow.jl/Wflow.spdx.json",
    rootpackages = dependencies,
    find_artifactsource = true,
    packageInstructions = Dict{UUID, spdxPackageInstructions}(
        Pkg.project().uuid => packageInstructions,
    ),
)

sbom = generateSPDX(spdxData)
writespdx(sbom, "Wflow.spdx.json")
