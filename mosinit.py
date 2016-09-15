from neuron import *

h("strdef simname, allfiles, simfiles, output_file, datestr, uname, osname, comment")
h.simname=simname = "mtlhpc"
h.allfiles=allfiles = "pyinit.py network.py"
h.simfiles=simfiles = "pyinit.py network.py"
h("runnum=1")
runnum = 1.0
h.datestr=datestr = "12aug9"
h.output_file=output_file = "data/12aug9.01"
h.uname=uname = "x86_64"
h.osname=osname="linux"
h("templates_loaded=0")
templates_loaded=0
h("xwindows=1.0")
xwindows = 1.0

h.xopen("nrnoc.hoc")
h.xopen("init.hoc")

from network import *

nq = anet.randruns(10000) # run the simulation with 10,000 random initializations

# print the apoptosis ratio, given the four inputs.
for i in xrange(2):
  for j in xrange(2):
    nq.verbose=0
    na = nq.select("TNF",i,"GF",j,"apop",1)
    ns = nq.select("TNF",i,"GF",j,"apop",0)
    print "TNF:",i,"GF:",j,"apop ratio:",na/(na+ns)
    nq.verbose=1
