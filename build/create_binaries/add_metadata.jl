"""
Add the following metadata files to the newly created build:

- Build.toml
- Project.toml
- Manifest.toml
- README.md
- LICENSE
- dep_licenses/
"""

function add_metadata(project_dir, license_file, output_dir, git_repo, sbom_file)
    # save some environment variables in a Build.toml file for debugging purposes
    vars = ["BUILD_NUMBER", "BUILD_VCS_NUMBER"]
    dict = Dict(var => ENV[var] for var in vars if haskey(ENV, var))
    open(normpath(output_dir, "share/julia/Build.toml"), "w") do io
        TOML.print(io, dict)
    end

    # a stripped Project.toml is already added in the same location by PackageCompiler
    # however it is better to copy the original, since it includes the version and compat
    cp(
        normpath(project_dir, "Project.toml"),
        normpath(output_dir, "share/julia/Project.toml");
        force = true,
    )
    # the Manifest.toml always gives the exact version of Wflow that was built
    cp(
        normpath(git_repo, "Manifest.toml"),
        normpath(output_dir, "share/julia/Manifest.toml");
        force = true,
    )

    # put the LICENSE in the top level directory
    cp(license_file, normpath(output_dir, "LICENSE"); force = true)
    cp(normpath(project_dir, "README.md"), normpath(output_dir, "README.md"); force = true)
    open(normpath(output_dir, "README.md"), "a") do io
        # since the exact Wflow version may be hard to find in the Manifest.toml file
        # we can also extract that information, and add it to the README.md
        manifest = TOML.parsefile(normpath(git_repo, "Manifest.toml"))
        julia_version = manifest["julia_version"]
        version = TOML.parsefile(normpath(git_repo, "Wflow/Project.toml"))["version"]
        short_name =
            chomp(read(Cmd(`git rev-parse --abbrev-ref HEAD`; dir = git_repo), String))
        short_commit =
            chomp(read(Cmd(`git rev-parse --short=10 HEAD`; dir = git_repo), String))

        # get the release from the current tag, like `git describe --tags`
        # if it is a commit after a tag, it will be <tag>-g<short-commit>
        described_tag = chomp(
            read(Cmd(`git describe --tags --dirty --abbrev=7`; dir = git_repo), String),
        )
        tag = replace(described_tag, r"^v" => "")

        url = "https://github.com/Deltares/Wflow.jl/tree"
        version_info = """

        ## Version

        This build uses the Wflow.jl version mentioned below.

        ```toml
        release = "$tag"
        version = "$version"
        commit = "$url/$short_commit"
        branch = "$url/$short_name"
        julia_version = "$julia_version"
        ```"""
        println(io, version_info)
    end

    cp(sbom_file, normpath(output_dir, "Wflow.spdx.json"); force = true)

    # collect licences of all dependencies
    ctx = PackageCompiler.create_pkg_context(project_dir)
    license_dir = joinpath(output_dir, "dep_licenses")
    mkpath(license_dir)

    for (uuid, pkg_entry) in ctx.env.manifest.deps
        if isnothing(pkg_entry.tree_hash)
            # seems as though stdlib packages don't have a tree_sha. as there doesn't seem to be
            # a different way of detecting those, I'm going to assume that is characteristic
            continue
        end
        install_path =
            Pkg.Operations.find_installed(pkg_entry.name, uuid, pkg_entry.tree_hash)

        license = LicenseCheck.find_license(install_path)
        if !isnothing(license)
            license_file_path = joinpath(install_path, license.license_filename)
            cp(license_file_path, joinpath(license_dir, pkg_entry.name); force = true)
        end
    end
end
