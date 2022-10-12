using Thrift
import Thrift.process, Thrift.meta

include("gen-jl/openda_bmi/openda_bmi.jl");

# create a client instance with our choice of protocol and transport
clnt_transport = TSocket(19999)
proto = TBinaryProtocol(clnt_transport)
clnt = openda_bmi.BMIServiceClient(proto)

# open a connection
open(clnt_transport)
# invoke service and print the result
openda_bmi.initialize(clnt, "..\\test\\sbm_config.toml")
# close connection
close(clnt_transport)
