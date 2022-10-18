using Thrift
import Thrift.process, Thrift.meta

global srvr

include("gen-jl/openda_bmi/openda_bmi.jl");
include("gen-jl/openda_bmi/openda_bmi_types.jl");

# create a client instance with our choice of protocol and transport
clnt_transport = TSocket(19999)
transport = TFramedTransport(clnt_transport)
proto = TBinaryProtocol(transport)
clnt = openda_bmi.BMIServiceClient(proto)

# open a connection
open(clnt_transport)

# invoke service
try
    openda_bmi.initialize(clnt, "..\\test\\sbm_config.toml")
catch ex
    !isa(ex, ModelException) && rethrow()
    println(ex.message)
end

try
    openda_bmi.update(clnt)
catch ex
    !isa(ex, ModelException) && rethrow()
    println(ex.message)
end

# close connection
close(clnt_transport)
