# Wflow ZMQ Server
Call [Wflow](https://github.com/Deltares/Wflow.jl) functions exposed through the Basic Model
Interface (BMI) implementation of Wflow and three additonal Wflow functions related to model
states and start time of the model in Unix time, using [ZeroMQ](https://zeromq.org/) with
TCP data transport. The Wflow ZMQ Server allows users to interact with a Wflow model from
many other programming languages with ZeroMQ bindings. An example is the use of
[OpenDA](https://openda.org/) Java software for data-assimilation. The coupling of OpenDA
and Wflow through a black box model approach is too slow because it requires a restart of
Wflow (multiple times), while this is not required with a client-server approach.

Specify the `WflowServer` environment on startup, in the current directory:
```
julia --project=.
```
and start the Wflow ZMQ Server, listening on port `<port>` (default: 5555) as follows:
```julia-repl
julia> using WflowServer
julia> WflowServer.start() <port>
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
