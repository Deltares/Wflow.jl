"""
Add the following metadata files to the newly created build:

- Build.toml
- Project.toml
- Manifest.toml
- README.md
- LICENSE
"""
function add_metadata(project_dir, license_file, output_dir, git_repo)
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
        normpath(project_dir, "Manifest.toml"),
        normpath(output_dir, "share/julia/Manifest.toml");
        force = true,
    )

    # put the LICENSE in the top level directory
    cp(license_file, normpath(output_dir, "LICENSE"); force = true)
    cp(normpath(project_dir, "README.md"), normpath(output_dir, "README.md"); force = true)
    open(normpath(output_dir, "README.md"), "a") do io
        # since the exact Wflow version may be hard to find in the Manifest.toml file
        # we can also extract that information, and add it to the README.md
        manifest = TOML.parsefile(normpath(project_dir, "Manifest.toml"))
        if !haskey(manifest, "manifest_format")
            error("Manifest.toml is in the old format, run Pkg.upgrade_manifest()")
        end
        julia_version = manifest["julia_version"]
        version = TOML.parsefile(normpath(git_repo, "Project.toml"))["version"]
        repo = GitRepo(git_repo)
        branch = LibGit2.head(repo)
        commit = LibGit2.peel(LibGit2.GitCommit, branch)
        short_name = LibGit2.shortname(branch)
        short_commit = string(LibGit2.GitShortHash(LibGit2.GitHash(commit), 10))

        # get the release from the current tag, like `git describe --tags`
        # if it is a commit after a tag, it will be <tag>-g<short-commit>
        options = LibGit2.DescribeOptions(; describe_strategy = LibGit2.Consts.DESCRIBE_TAGS)
        result = LibGit2.GitDescribeResult(repo; options)
        suffix = "-dirty"
        foptions =
            LibGit2.DescribeFormatOptions(; dirty_suffix = Base.unsafe_convert(Cstring, suffix))
        GC.@preserve suffix tag = LibGit2.format(result; options = foptions)[2:end]  # skip v prefix

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
end
