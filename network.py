# $Id: network.py,v 1.34 2012/08/09 16:29:02 samn Exp $ 

from pyinit import *

h.load_file("pywrap.hoc") 

from bnet import *

#######################################
#      some utils to avoid the h.     #
vlk = h.vlk
Vector = h.Vector
NQS = h.NQS
gg = h.gg
ge = h.ge
gg()
g = h.g[0]
Random = h.Random
List = h.List
Matrix = h.Matrix
nqsdel = h.nqsdel
#######################################

class apopnet(bnet):
    "apoptosis boolean network"    
    def __init__ (self,fnames="apopnames.txt",frules="apoprules.txt"):
        bnet.__init__(self,fnames,frules)

    # run boolean network bn for niters using starting vector vstart
    # save animation to a gif using Graphviz and ImageMagick convert
    # requires Graphviz/ImageMagick installation and a directory called
    # frames in the current working directory to save output to. this
    # function will also write to a temporary file called __junk__.dot
    def myanim (self,niters,fgif):
        self.start()
        os.system("rm frames/frame*.gif") # get rid of old gif frames
        for i in xrange(niters):
            fn = "frames/frame." + ("%03d" % i) + ".gif"
            self.bn.graphviz("__junk__.dot",fn,"gif",1,15,10,20) # must have Graphviz installed
            self.bn.advancebn()
        os.system("convert -delay 200 -loop 0 frames/frame*.gif %s" % fgif) # uses ImageMagick convert

    # reachedapop - return True when the network has been run to an apoptotic state
    def reachedapop (self):
        idx = self.idnames["Apoptosis"]
        for i in xrange(len(self.lstate)):
            if self.lstate[i][idx] == 1: return True
        return False

    # survived - return True when network has not reached an apoptotic state
    def survived (self):
        return not self.reachedapop()

    # randruns - run for niters with different random starting states
    #  each run consists of nsteps advances
    #  returns an NQS with iteration number, starting/final state, and whether apoptotic
    def randruns (self,niters,nsteps=200,rseed=1234):
        nq=h.NQS("iter","TNF","GF","vstart","vfinal","apop")
        nq.odec("vstart")
        nq.odec("vfinal")
        nq.clear(niters*4) # reserve space
        self.setrand(rseed)
        for i in xrange(niters):
            if i % 100 == 0: print "init %d of %d" % (i, niters)
            self.randvstart() # randomize starting state
            for TNFState in xrange(2):
                self.vstart.x[self.idnames["TNF"]] = TNFState # sets TNF state
                for GFState in xrange(2):
                    self.vstart.x[self.idnames["GF"]] = GFState # sets GF state
                    self.startrun(nsteps) # start the network up and run it
                    nq.append(i,TNFState,GFState,self.vstart,self.lstate[-1],self.reachedapop())
        return nq

anet = apopnet() # make the apoptosis network

