# How to install

As Julia is still a relatively new programming language, a basic "getting started with
Julia" guide can be found [here](https://docs.julialang.org/en/v1/manual/getting-started/).

Wflow can be used in two different ways, depending on the required use of the code:

- If you want to stay up-to-date with the latest version, explore and modify the model code,
  and write your own Julia scripts around the wflow package, we recommend installing wflow
  as a [Julia](https://julialang.org/) package.
- If you don't need extra features, but just want to run simulations, a complied executable
  version is available. This consists of a single executable, `wflow_cli`, allowing you to
  run the model via the command line.

Below, we describe how to install both versions of wflow. 

## Installing as Julia package

Wflow is a [Julia](https://julialang.org/) package, that can be installed through several 
different ways. Below we show how to install wflow from Julia's package repository, and
how to install the latest version from GitHub. 

### Install from Julia's package repository 

To access Julia's package manager, press `]` in the Julia REPL.  To get back to the Julia 
REPL, press backspace or ^C. 

!!! tip
    If you have not used Julia in a while, it can be a good idea to run `up` to update your
    packages.
    ```julia-repl
    pkg> up
    ```

Access Julia's package manager and install wflow using:
```julia-repl
pkg> add Wflow
```

This can take a while, especially the first time, since compatible dependencies are also
automatically looked up and installed from the Pkg General registry. 

### Install from GitHub

It is also possible to install wflow from the `master` branch as follows:

```julia-repl
pkg> add https://github.com/Deltares/Wflow.jl
```

This command will track the `master` branch, and will update to the latest commit on that
branch when you run `update`, or simply `up`, in the Pkg REPL. The use of `add` will install
wflow in you home directory under `.julia/packages/Wflow`. Note that packages installed
under `packages` by `add` are supposed to never be altered in that location, for Pkg and
it's automatic dependency handling to work well. 

If you want to make any changes to any of the files in the wflow directory, you need to do 
a development install. This can be done using:

```julia-repl
pkg> dev https://github.com/Deltares/Wflow.jl
```

This will clone the git repository, put it under your home directory in `.julia/dev/Wflow`,
and add the wflow package to your project environment. Note that to receive updates, you
have to pull in the latest changes yourself using `git pull`.

### Check installation of wflow

Finally, go back the Julia REPL and try to load wflow:

```julia-repl
julia> using Wflow
```

The first time this will take longer as any package that is new or changed needs to be
pre-compiled first, to allow faster loading on subsequent uses. No error messages should 
appear, indicating that you have now successfully installed wflow. 

Before ending this section, we still want to recommend a few tools that can make using and
developing Julia code easier.

!!! tip
    There is a section on editors and IDEs for Julia on <https://julialang.org/>, scroll
    down to see it. We use and recommend Microsoft's free and open source [Visual Studio
    Code](https://code.visualstudio.com/). When combined with the [Julia
    extension](https://www.julia-vscode.org/) it provides a powerful and interactive
    development experience.

!!! tip
    When planning to make changes to the code of wflow, we recommend installing the `Revise.jl`
    package. This package allows you to modify code and use the changes without restarting 
    Julia. Install it with `add Revise` from the Pkg REPL. Then create a file called
    `.julia/config/startup.jl`, and put `using Revise` there. This will load Revise every 
    time you start a Julia session.

## Installing the compiled executable

Binaries of `wflow_cli` can be downloaded from our website
[download.deltares.nl](https://download.deltares.nl/en/download/wflow/), and are currently
available for Windows. Download and install the `.msi` file. After installing you can see
two folders in the installation directory. It is only the `bin/wflow_cli` that is used. The
artifacts folder contains binary dependencies such as NetCDF.

```
artifacts\
bin\wflow_cli
```

Check whether the installation was performed successfully, run `wflow_cli` with no 
arguments in the command line will give the following message:

```
Usage: wflow_cli 'path/to/config.toml'
```


!!! note
    The old version of wflow, based on Python and PCRaster libraries is also available to
    download from our website [download.deltares.nl](https://download.deltares.nl/en/download/wflow/).
    We recommend installing the Julia version, as this documentation is written to support
    this version. 



