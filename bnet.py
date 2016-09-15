# $Id: bnet.py,v 1.30 2012/08/07 15:25:30 samn Exp $ 

from neuron import h
from snutils import unique
try:
    import networkx as nx
    havenetworkx = True
except:
    print "Could not import networkx!"
    havenetworkx = False

# read names from file, one per line and return nqs
# lines that have 2 elements (name number) mean the node is 'special' and
# has a threshold for turning on (the number)
def readnames (fname):
    nqn = h.NQS("id","name","sthresh","input")
    nqn.strdec("name")
    idx = 0
    fp = open(fname)
    for line in fp.readlines():
        line = line.strip() # gets rid of trailing \n
        if line.startswith("//"): continue # skip comments
        line = line.split()
        sthresh = 0
        isinput = 0
        if len(line) > 1:
            if line[1].startswith("input"):
                isinput = 1
            else:
                sthresh = int(line[1])
        nqn.append(idx,line[0],sthresh,isinput)
        idx = idx + 1
    fp.close()
    return nqn

# lookup the id using the name of the node
def lookupid (nqn,name):
    nqn.verbose = 0.0
    idx = -1
    if 1 == nqn.select("name",h.SEQ,name):
        idx = int(nqn.getcol("id").x[0])
    nqn.tog("DB")
    nqn.verbose = 1.0
    return idx

# nqs with rules - reads from text file
def readrules (fname,nqn):
    fp = open(fname)
    nqr = h.NQS("ridx","vsrc","target","weight","vsrcstate","targname")
    nqr.odec("vsrc")
    nqr.odec("vsrcstate")
    nqr.strdec("targname")
    ridx = 0 # rule identifier
    for line in fp.readlines():
        line   = line.strip()
        if line.startswith("//"): continue # skip comments
        if len(line) < 1: continue
        ln     = line.split("->")
        src    = ln[0].strip()
        targ   = ln[1].strip()
        targid = lookupid(nqn,targ)
        [src, srcstate, weight] = src.split(":")
        src = src.split(",")
        srcstate = srcstate.split(",")
        weight = int(weight)
        vsrc   = h.Vector()
        vsrcstate = h.Vector()
        for sname in src:
            srcid = lookupid(nqn,sname)
            vsrc.append(srcid)
        for sstate in srcstate:
            if sstate == "ON":
                vsrcstate.append(1.0)
            else:
                vsrcstate.append(0.0)
        nqr.append(ridx,vsrc,targid,weight,vsrcstate,targ)
        ridx = ridx + 1
    fp.close()
    return nqr

# makenet: make the boolean network, set connections/start state/names
def makenet (nqn,nqr,vstart):
    bnet = h.BNET(int(nqn.v[0].size())) # make a boolean network with specified nodes
    nqn.tog("DB")
    for i in xrange(int(nqn.v[0].size())): # store the names in the BNET
        bnet.setnname(i,nqn.get("name",i).s)
    bnet.setsthresh(nqn.getcol("sthresh"))#set sthresh for all nodes (node is special iff sthresh > 0)
    # setup the connectivity
    nqr.tog("DB")
    for i in xrange(int(nqr.v[0].size())):
        vsrc = nqr.get("vsrc",i).o[0]
        targid = nqr.get("target",i).x[0]
        weight = nqr.get("weight",i).x[0]
        vsrcstate = nqr.get("vsrcstate",i).o[0]        
        bnet.setrule(vsrc,targid,weight,vsrcstate)
    bnet.setstart(vstart)   # set the starting state
    bnet.start()            # set the state to start
    return bnet


# bnet class: python wrapper for the BNET.mod ARTIFICIAL_CELL template
#  fnames is name of file with names of nodes and their properties
#  frules is name of file with rules
#  rseed is random number seed (optional)
class bnet:
    "boolean network - python wrapper over NEURON BNET.mod ARTIFICIAL_CELL template"
    def __init__(self,fnames,frules,rseed=1234):
        self.nqn = readnames(fnames)
        self.nqr = readrules(frules,self.nqn)
        self.vstate = h.Vector(self.nqn.v[0].size())
        self.vstart = h.Vector(self.nqn.v[0].size())#all nodes start @ 0 by default
        self.bn = makenet(self.nqn,self.nqr,self.vstart)
        self.lstate = []
        self.bn.getstate(self.vstate)
        self.makeNameDict() # map names to ids and vice-versa
        global havenetworkx
        if havenetworkx:
            self.makeNXG() # make a networkx Graph representation
        self.setspecial()
        self.setrand(rseed)

    # randvstart - pick random starting state using self.rdm
    #  avoids changing the special node starting states
    def randvstart (self):
        vtmp = h.Vector(self.bn.numnodes())
        vtmp.copy(self.vstart)
        self.rdm.discunif(0,1) # off (0) or on (1)
        self.vstart.setrand(self.rdm)
        # make sure special node starting states not modified
        for i in xrange(int(self.vspecial.size())):
            if self.vspecial.x[i]: self.vstart.x[i]=vtmp.x[i]            

    # setrand - setup the random number generator
    def setrand (self,rseed):
        self.rseed = rseed
        self.rdm = h.Random() # random number generator
        self.rdm.ACG(self.rseed)        

    # setup which nodes are special in 
    def setspecial (self):
        self.vspecial = h.Vector(self.bn.numnodes())
        self.nqn.tog("DB")
        for i in xrange(int(self.vspecial.size())):
            if self.nqn.getcol("input").x[i] != 0 or self.nqn.getcol("sthresh").x[i] > 0:
                self.vspecial.x[i] = 1

    # make a NX Graph object representation
    def makeNXG (self):
        self.nxg = nx.MultiDiGraph() # multiple parallel edges in a directed graph
        for n in self.idnames.keys():
            if type(n) == str:
                self.nxg.add_node(n)
        for idx in xrange(int(self.nqr.v[0].size())):
            vsrc = self.nqr.get("vsrc",idx).o[0]
            tn = self.idnames[int(self.nqr.get("target",idx).x[0])]
            w = self.nqr.get("weight",idx).x[0]
            for sdx in vsrc:
                sn = self.idnames[int(sdx)]
                self.nxg.add_edge(sn,tn,weight=w)

    # makeNameDict - make a dictionary with name -> id and id -> name mappings
    def makeNameDict (self):
        self.idnames = {}
        self.nqn.tog("DB")
        for i in xrange(int(self.nqn.v[0].size())):
            self.idnames[int(self.nqn.getcol("id").x[i])] = self.nqn.get("name",i).s
            self.idnames[self.nqn.get("name",i).s] = int(self.nqn.getcol("id").x[i])

    # get IDs of all sources -> targ
    def SRCIDs (self,targ):
        src = []
        if type(targ) == int:
            if self.nqr.select(-1,"target",targ) > 0:
                for idx in self.nqr.ind:
                    for sdx in self.nqr.get("vsrc",idx).o[0]:
                        src.append(int(sdx))
        elif type(targ) == str:
            if self.nqr.select(-1,"targname",h.SEQ,targ) > 0:
                for idx in self.nqr.ind:
                    for sdx in self.nqr.get("vsrc",idx).o[0]:
                        src.append(int(sdx))
        return unique(src)

    # get names of all source -> targ
    def SRCNames (self,targ):
        srcid = self.SRCIDs(targ)
        src = []
        for s in srcid: src.append( self.idnames[s] )
        return src

    # get IDs of all targets from src 
    def TARGIDs (self,src):
        targ = []
        nqr = self.nqr
        nqr.tog("DB")
        if type(src) == int:
            for idx in xrange(int(nqr.v[0].size())):
                vsrc = nqr.get("vsrc",idx).o[0]
                targid = int(nqr.get("target",idx).x[0])
                for sdx in vsrc:
                    if int(sdx) == src: targ.append(targid) # same source? save target
        elif type(src) == str:
            for idx in xrange(int(nqr.v[0].size())):
                vsrc = nqr.get("vsrc",idx).o[0]
                targid = int(nqr.get("target",idx).x[0])
                for sdx in vsrc:
                    if self.idnames[int(sdx)] == src: targ.append(targid) # same source? save target
        return unique(targ)

    # get names of all targets from src
    def TARGNames (self,src):
        targid = self.TARGIDs(src)
        targ = []
        for tt in targid: targ.append( self.idnames[tt] )
        return targ

    # advance the network by 1 timestep
    def advance (self):
        tt = self.bn.advancebn()
        self.bn.getstate(self.vstate) # update state vector
        return tt

    # check if a node is on(1) or off(0)
    def nodestate (self,name):
        idx = self.idnames[name]
        return self.vstate.x[idx]

    # run the boolean network bn for niters iterations and sets the
    # states in the global list of state vectors (lstate)
    def run (self,niters=1):
        if niters < 1: return
        vstate=h.Vector() # NB: vstate not class instance here
        self.bn.getstate(vstate)
        self.lstate.append(vstate) # initial state
        for i in xrange(niters):
            self.bn.advancebn() # iterate
            vstate=h.Vector()
            self.bn.getstate(vstate)            
            self.lstate.append(vstate)
        self.vstate.copy(self.lstate[-1]) # set vstate to last state

    # put self.bn into starting state
    def start (self):
        self.lstate = []
        self.bn.setstart(self.vstart)
        self.bn.start()

    # run the boolean network bn for niters iterations and sets the
    # states in the global list of state vectors (lstate)
    def startrun (self,niters=1):
        self.start()
        self.run(niters)

    # print states and rules
    def pr (self):
        self.bn.pr()

