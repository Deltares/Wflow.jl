# Quick start

## Installation

This is a [Julia](https://julialang.org/) package, that can be installed using Julia's
built-in package manager, Pkg. Wflow is currently still unregistered, so the full URL is
needed when adding the package, `add Wflow` does not yet work, but `add
https://github.com/Deltares/Wflow.jl` does. In the rest of this section will provide more
detailed instructions, such that new Julia users can also install Wflow. If you already know
you can go ahead to the next section.

First download and install the [current stable release of
Julia](https://julialang.org/downloads/#current_stable_release). Please see [platform
specific instructions](https://julialang.org/downloads/platform/) for further installation
instructions and if you have trouble installing Julia.

If you are new to Julia, it may be a good idea to check out the [Getting Started section of
the Julia Manual](https://docs.julialang.org/en/v1/manual/getting-started/). Links to other
learning resources can be found at https://julialang.org/learning/.

Start an interactive session (also known as a read-eval-print loop or "REPL") by
double-clicking the Julia executable or running `julia` from the command line, in case you
have chosen to add it to your path. The default prompt will appear, in which Julia code can
be entered and run.

```julia-repl
julia> 1 + 2
3
```

Enter the Pkg REPL by pressing `]` from the Julia REPL. To get back to the Julia REPL, press
backspace or ^C. Just like in the Julia REPL, you can type `?` to get a help message, this
will list the available commands. To get more information on a specific command type `?`
followed by the command.

```julia-repl
pkg> ?

pkg> ?add
```

If you go back to the Julia REPL, note that the same works, and is a convenient way to
consult documentation for specific functions or types, or other objects.

```
help?> cos
```

Now go back to the Pkg REPL and install Wflow using:

```julia-repl
pkg> add https://github.com/Deltares/Wflow.jl
```

This can take a while, especially the first time, since compatible dependencies are also
automatically looked up and installed from the Pkg General registry. This command will track
the `master` branch, and will update to the latest commit on that branch when you run
`update`, or simply `up`, in the Pkg REPL. This will install Wflow in you home directory
under `.julia/packages/Wflow`. Note that packages installed under `packages` by `add` are
supposed to never be altered in that location, for Pkg and it's automatic dependency
handling to work well. If you want to make any changes to any of the files in the Wflow
directory, you need to do a development install. This can be done using:

```julia-repl
pkg> dev https://github.com/Deltares/Wflow.jl
```

This will clone the git repository, put it under your home directory in `.julia/dev/Wflow`,
and add the `Wflow` package to your project environment. Note that to receive updates, you
have to pull in the latest changes yourself using `git pull`.

Finally, go back the Julia REPL and try to load Wflow:

```julia-repl
julia> using Wflow
```

The first time this will take longer as any package that is new or changed needs to be
precompiled first, to allow faster loading on subsequent uses.

You have now succesfully installed Julia and Wflow. Before ending this section, we still
want to recommend a few tools that can make using and developing Julia code easier.

The first is `Revise.jl`. This package allows you to modify code and use the changes without
restarting Julia. Install it with `add Revise` from the Pkg REPL. Then create a file called
`.julia/config/startup.jl`, and put `using Revise` there. This will load Revise everytime
you start a Julia session.

There is a section on editors and IDEs for Julia on https://julialang.org/, scroll down to
see it. We use and recommend Microsoft's free and open source [Visual Studio
Code](https://code.visualstudio.com/). When combined with the [Julia
extension](https://www.julia-vscode.org/) it provides a powerful and interactive development
experience.

# Usage

To be added promptly.
