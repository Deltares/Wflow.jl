using Thrift
import Thrift.process, Thrift.meta
using Wflow

include("gen-jl/openda_bmi/openda_bmi.jl");

# create a server instance with our choice of protocol and transport
srvr_processor = openda_bmi.BMIServiceProcessor()
srvr_transport = TServerSocket(19999)
#srvr = TSimpleServer(srvr_transport, srvr_processor, x->x, x->TBinaryProtocol(x), x->x, x->TBinaryProtocol(x))
srvr = TSimpleServer(srvr_transport, srvr_processor, x->TFramedTransport(x), x->TBinaryProtocol(x), x->TFramedTransport(x), x->TBinaryProtocol(x))

# start serving client requests
println("starting to serve requests...")
serve(srvr)
println("server stopped")