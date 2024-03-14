# Wflow ZMQ Server
Call [Wflow](https://github.com/Deltares/Wflow.jl) functions exposed through the Basic Model
Interface (BMI) implementation of Wflow and three additonal Wflow functions related to model
states and start time of the model in Unix time, using [ZeroMQ](https://zeromq.org/) with
TCP data transport. The Wflow ZMQ Server allows users to interact with a Wflow model from
many other programming languages with ZeroMQ bindings. An example is the use of
[OpenDA](https://openda.org/) Java software for data-assimilation. The coupling of OpenDA
and Wflow through a black box model approach is too slow because it requires a restart of
Wflow (multiple times), while this is not required with a client-server approach.

The Wflow ZMQ Server can be started through an interactive Julia version, using different
functions. First, start the `WflowServer` environment on startup, in the current directory:
```
julia --project=.
```
then start the Wflow ZMQ Server, listening on port 5556 with the `start` function as
follows:
```julia-repl
julia> using WflowServer
julia> WflowServer.start(5556)
```
or start the Wflow ZMQ Server, listening on port 5556 with the `main` function:
```julia-repl
julia> using WflowServer
julia> WflowServer.main(["port=5556"])
```
or start the Wflow ZMQ Server with the `main` function using the default port number 5555:
```julia-repl
julia> using WflowServer
julia> WflowServer.main()
```
Finally, it is also possible to start the Wflow ZMQ server directly from the command line
with the Julia script `run_server.jl` in the current directory, as follows:
```
julia --project=. run_server.jl --port 5556
```

## JSON
JSON is used for data serialization. The Wflow ZMQ Server maps a JSON message to a Julia
structure directly. Below examples of two messages that provide a function name `fn` and
other arguments required by the exposed Wflow functions (mostly through BMI):

```
# initialize a Wflow model through the configuration file `sbm_config.toml`
"""{"fn": "initialize", "config_file": "sbm_config.toml"}"""
# update model until time 86400.0
"""{"fn": "update_until", "time": 86400.0}"""
```

The [tests](/server/test/) provide further examples of the interface to the Wflow (BMI) and
server functions.
