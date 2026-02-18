import Pkg

const IS_INSTALLED = r"\s*Installed (.+)"
const DETAILS_BEGIN = """
<details>
<summary>
All package versions
</summary>

```
"""

"""
Update the Julia Manifest.toml and show the changes as well as outdated packages.
The output is written to a file that can be used as the body of a pull request.
"""
function (@main)(_)
    path = normpath(@__DIR__, "../.pixi/update-manifest-julia.md")
    redirect_stdio(; stdout = path, stderr = path) do
        println("Update the Julia Manifest.toml to get the latest dependencies.\n")
        println("__Changed packages__\n```")
        Pkg.update()
        println("```\n\n__Packages still outdated after update__\n```")
        Pkg.status(; outdated = true)
        println("```")
    end

    # The Pkg.update output first prints all installed package versions.
    # This is a lot, strip it out, sort it, and put it in a details tag at the end.
    installed_lines = String[]
    lines = readlines(path)
    return open(path, "w") do io
        for line in lines
            m = match(IS_INSTALLED, line)
            if m === nothing
                println(io, line)
            else
                push!(installed_lines, only(m.captures))
            end
        end

        println(io, DETAILS_BEGIN)
        sort!(installed_lines)
        foreach(line -> println(io, line), installed_lines)
        println("```\n\n</details>")
    end
end
