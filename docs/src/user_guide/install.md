# How to install

First, download and install the [current stable release of
Julia](https://julialang.org/downloads/#current_stable_release). 
If you have any issues installing Julia, please see [platform
specific instructions](https://julialang.org/downloads/platform/) for further instructions.

If you are new to Julia, it might be a good idea to check out the [Getting Started section of
the Julia Manual](https://docs.julialang.org/en/v1/manual/getting-started/).
You can also find additional learning resources at
[julialang.org/learning](https://julialang.org/learning/).


Wflow can be used in two different ways, depending on the required use of the code:

- If you want to stay up-to-date with the latest version, explore and modify the model code,
  and write your own Julia scripts around the wflow package, we recommend installing wflow
  as a [Julia](https://julialang.org/) package.
- If you don't need extra features, but just want to run simulations, a compiled executable
  version is available. This version includes a single executable, `wflow_cli`, which allows you to
  run the model via the command line.

Below we describe how to install both versions of wflow.

## Installing as Julia package

Wflow is a [Julia](https://julialang.org/) package that can be installed in several ways. 
Below, we show how to install wflow from Julia's package repository and
how to install the latest version from GitHub.

### Install from Julia's package repository

To access Julia's package manager, press `]` in the Julia REPL.  To get back to the Julia
REPL, press backspace or ^C.

!!! tip
    If you haven't used Julia in a while, it's a good idea to run `up` to update your
    packages.
    ```julia-repl
    pkg> up
    ```

To access Julia's package manager and install wflow, use:
```julia-repl
pkg> add Wflow
```

This process can take a while, especially on the first run, as compatible dependencies are
automatically resolved and installed from the Pkg General registry.

### Install from GitHub

You can also install wflow from the `master` branch on the repository as follows:

```julia-repl
pkg> add Wflow#master
```

This command tracks the `master` branch and updates to the latest commit on that
branch when you run `update`, or simply `up`, in the Pkg REPL. The `add` installs
wflow in your home directory under `.julia/packages/Wflow`. Note that packages installed
under `packages` by `add` should not be changed in the directory, 
as the change could disrupt Pkg's automatic dependency handling.

If you want to modify any files in the repository, you need to do
a development install. This can be done using:

```julia-repl
pkg> dev Wflow
```

This will clone the git repository, place it under your home directory in `.julia/dev/Wflow`,
and add the wflow package to your project environment. To receive updates, 
you'll need to pull the latest changes manually using `git pull`.

### Check installation of wflow

Finally, go back to the Julia REPL and try to load wflow:

```julia-repl
julia> using Wflow
```

The first time you do this, it may take longer as any new or changed packages need to be
precompiled to enable faster loading on subsequent uses. No error messages should
appear, which indicates that you have successfully installed wflow.

Before concluding this section, we recommend a few tools that can make using and
developing Julia code easier.

!!! tip
    There is a section on editors and IDEs for Julia on <https://julialang.org/>, scroll
    down to see it. We use and recommend Microsoft's free and open source [Visual Studio
    Code](https://code.visualstudio.com/). Combined with the [Julia
    extension](https://www.julia-vscode.org/) it provides a powerful and interactive
    development experience.

!!! tip
    If you plan to modify the code of wflow, we recommend installing the `Revise.jl`
    package. This package allows you to modify code and use the changes without restarting
    Julia. Install it with `add Revise` from the Pkg REPL. Then create a file called
    `.julia/config/startup.jl`, and put `using Revise` there. This will load Revise every
    time you start a Julia session.

## Installing the compiled executable

Binaries of `wflow_cli` can be downloaded from our website
[download.deltares.nl](https://download.deltares.nl/en/download/wflow/), and are currently
available for Windows. Download and install the `.msi` file. After installation, you will see
two folders in the installation directory. Only the `bin/wflow_cli` is used. The
`artifacts` folder contains binary dependencies such as netCDF.

```
artifacts\
bin\wflow_cli
```

To verify whether the installation was completed successfully, run `wflow_cli` with no
arguments in the command line. This will display the following message:

```
Usage: wflow_cli 'path/to/config.toml'
```


!!! note
    The old version of wflow, which was based on Python and PCRaster libraries, is also available for
    download from our website [download.deltares.nl](https://download.deltares.nl/en/download/wflow/).
    We recommend installing the Julia version, as this documentation is written to support
    this version.
