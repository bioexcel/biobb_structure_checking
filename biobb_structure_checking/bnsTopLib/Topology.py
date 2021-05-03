import bnsTopLib.StructureWrapper

from Bio.PDB.NeighborSearch import NeighborSearch

class Topology():

    COVLNK = 2.0
    HBLNK  = 3.5
    INTDIST = 5.0
    BPTHRESDEF = 1.0

    def __init__(self,args):
        self.covLinkPairs=[]
        self.args = args
        self.conts={}
        self.intList={}
        self.interfPairs={}
        self.wcats = []
        self.bps = []
        self.bpsteps=[]
        self.HFSeqs=[]
        self.bpsref={}
        self.notInHF=set()

        
    def calcCovLinks(self, st):
        bckats = [] 
        for at in st.get_atoms():
            if '.'+at.id in bnsTopLib.StructureWrapper.backbone_link_atoms:
                bckats.append(at)

        nbsearch = NeighborSearch(bckats)

        for at1, at2 in nbsearch.search_all(Topology.COVLNK):
            if Topology.sameResidue(at1,at2):
                continue
            self.covLinkPairs.append (Topology.getOrderedResiduePair(at1,at2,self.args.useChains,self.args.useModels))
        if self.args.debug:
            print ("#DEBUG: covalently linked residue pairs")
            for r in sorted(self.covLinkPairs, key=lambda i: i[0].residue.index):
                print ("#DEBUG: ",r[0].resid(),r[1].resid())        

    def calcContacts(self):

        for ch1 in self.chList.getSortedSets():
            self.conts[ch1]={}
            self.intList[ch1]={}
            self.interfPairs[ch1]={}
            for ch2 in self.chList.getSortedSets():
                if ch2.ini <= ch1.ini:
                    continue
                if ch1.type != 'prot' and ch2.type !='prot':
                    continue
                self.conts[ch1][ch2]=[]
                self.interfPairs[ch1][ch2]=[]
                ats1 = []
                s1 = ch1._getResidues()
                s2 = ch2._getResidues()
                for r1 in s1:
                    for r2 in s2:
                        cont = s1[r1].getClosestContact(s2[r2],self.args.intdist)
                        if cont:
                            [at1,at2,d] = cont
                            [atom1,atom2] = Topology.getOrderedAtomPair(at1,at2,self.args.useChains,self.args.useModels)
                            [res1,res2] = Topology.getOrderedResiduePair(at1,at2,self.args.useChains, self.args.useModels)                            
                            self.interfPairs[ch1][ch2].append([res1,res2])
                        contp = s1[r1].getClosestPolarContact(s2[r2],Topology.HBLNK)
                        if contp:
                            [at1,at2,d] = contp
                            [atom1,atom2] = Topology.getOrderedAtomPair(at1,at2,self.args.useChains,self.args.useModels)
                            [res1,res2] = Topology.getOrderedResiduePair(at1,at2,self.args.useChains, self.args.useModels)                            
                            self.conts[ch1][ch2].append([atom1,atom2,d])
                self.intList[ch1][ch2] = bnsTopLib.ResidueSet.ResidueSetList(self.interfPairs[ch1][ch2])
        
    def checkExistsNA(self):
        hbAtoms = set()
        for hb in bnsTopLib.StructureWrapper.hbonds['WC']:
            for rat in hb:
                hbAtoms.add(rat)
        for ch in self.chList.getSortedSets():
            if ch.type != 'na':
                continue
            for r in ch.items:
                for at in r.residue.get_atoms():
                    if bnsTopLib.StructureWrapper.Atom(at).attype() in hbAtoms:
                        self.wcats.append(at)
        self.NAOk = len(self.wcats)
        return self.NAOk
        
    def calcBasePairs(self):
        wc_nbsearch = NeighborSearch(self.wcats)

        nhbs={}
        for at1, at2 in wc_nbsearch.search_all(bnsTopLib.Topology.HBLNK):
            if bnsTopLib.Topology.sameResidue(at1,at2):
                continue

            [atom1,atom2] = bnsTopLib.Topology.getOrderedAtomPair(at1,at2,self.args.useChains)

            res1 = bnsTopLib.StructureWrapper.Residue(atom1.at.get_parent(),self.args.useChains,self.args.useModels)
            res2 = bnsTopLib.StructureWrapper.Residue(atom2.at.get_parent(),self.args.useChains,self.args.useModels)

            if [atom1.attype(), atom2.attype()] not in bnsTopLib.StructureWrapper.hbonds['WC'] and \
                [atom2.attype(), atom1.attype()] not in bnsTopLib.StructureWrapper.hbonds['WC']:
                continue

            if [res1,res2] in self.covLinkPairs:
                continue

            if res1 not in nhbs:
                nhbs[res1] = {}

            if res2 not in nhbs[res1]:
                nhbs[res1][res2]=0

            nhbs[res1][res2]=nhbs[res1][res2]+atom1._hbscore(atom2)

            if self.args.debug:
                print ("#DEBUG: ", res1,res2,atom1._hbscore(atom2))

        if self.args.debug:
            print ("#DEBUG: HB count per pair of residues")
            for r1 in nhbs.keys():
                for r2 in nhbs[r1].keys():
                    print ("#DEBUG: ", r1, r2, nhbs[r1][r2])
        
        for r1 in nhbs.keys():
            maxv=0.
            pair=''
            for r2 in nhbs[r1].keys():
                if nhbs[r1][r2]>maxv:
                    pair=r2
                    maxv=nhbs[r1][r2]
            if maxv > self.args.bpthres:
                self.bps.append(bnsTopLib.StructureWrapper.BPair(r1,pair,maxv))

        for bp in sorted(self.bps):
            self.bpsref[bp.r1.residue.index]=bp
        
    def calcBPSteps(self):

        for bp1 in sorted(self.bps):        
            for bp2 in sorted(self.bps):
                if bp1 < bp2:
                    if bp1.r1.residue.index == bp2.r1.residue.index-1 and \
                    bp1.r2.residue.index == bp2.r2.residue.index+1:
                        self.bpsteps.append(bnsTopLib.StructureWrapper.BPStep(bp1,bp2))

        self.bpstpref={}

        for bpstp in sorted(self.bpsteps):
            self.bpstpref[bpstp.bp1.r1.residue.index]=bpstp
        
    def calcHelicalFrags(self):
        bpstepPairs=[]
        for bpst1 in sorted(self.bpsteps):
            for bpst2 in sorted(self.bpsteps):
                if bpst1 < bpst2:
                    if bpst1.bp1.r1.residue.index == bpst2.bp1.r1.residue.index-1:
                        bpstepPairs.append([bpst1,bpst2])

        self.HFragList=bnsTopLib.ResidueSet.BPSSetList(bpstepPairs)

        for ch in self.chList.getSortedSets():
            if ch.type != 'na':
                continue
            for r in ch.items:
                self.notInHF.add(r)

        for fr in self.HFragList.getSortedSets():
            frag={}
            for i in range (fr.ini,fr.fin+1):
                for bb in self.bpstpref[i].comps():
                    frag[bb]=1
            seq=[]
            for bp in sorted(frag.keys(), key=bnsTopLib.StructureWrapper.BPair.__index__):
                seq.append(bp.bpid())
                self.notInHF.remove(bp.r1)
                self.notInHF.remove(bp.r2)
                
            self.HFSeqs.append(seq)
        
            
        
    
    def json(self):
        jsondata = bnsTopLib.JSONWriter()   

        jsondata.set('useChains',self.args.useChains)
        jsondata.set('inputFormat', self.args.format)
        jsondata.set('id', self.args.id)

        jsondata.set('NChains', self.chList.n)
        for s in self.chList.getSortedSets():
            jsondata.appendData('chains',s.json())
        
        for r in sorted(self.covLinkPairs, key=lambda i: i[0].residue.index):
            jsondata.appendData('covLinks', [r[0].resid(True),r[1].resid(True)])
# Contacts & interfaces
        for ch1 in self.conts:
            for ch2 in self.conts[ch1]:
                for c in self.conts[ch1][ch2]:
                    jsondata.appendData('contacts', 
                        {
                            'ats':[c[0].atid(True), c[1].atid(True)],
                            'distance':float(c[2])
                        }
                    )
        
        intData={'cutoff':self.args.intdist}
        for ch1 in self.intList:
            for ch2 in self.intList[ch1]:
                intDataList=[]
                for s in self.intList[ch1][ch2].getSortedSets():
                    intDataList.append(','.join(s.getResidueIdList()))
                intData[str(ch1.id)+'-'+str(ch2.id)]=intDataList
        jsondata.set('interfaces',intData) 
        if self.NAOk:
            for bp in sorted(self.bps):
                jsondata.appendData('bpList', bp.json())
    
# Bpair steps from neighbour bps, relays on residue renumbering. TODO Check connectivity
            for bpstp in sorted(self.bpsteps):
                jsondata.appendData('bpStepList', bpstp.json())

# Continuous helical segments from stretches of overlapping bsteps
            for frs in self.HFSeqs:
                jsondata.appendData('HelicalFrags', frs)
            for r in self.notInHF:
                jsondata.appendData('ResiduesNotHF',r.resid(1))
        return jsondata
    
#===============================================================================    
    def sameResidue(at1,at2):
        return at1.get_parent() == at2.get_parent()

    def getOrderedResiduePair(at1,at2,useChains=False, useModels=False):
# Defining residues as res1 < res2
        res1 = bnsTopLib.StructureWrapper.Residue(at1.get_parent(), useChains,useModels)
        res2 = bnsTopLib.StructureWrapper.Residue(at2.get_parent(), useChains,useModels)
        if res1 < res2:
            return [res1,res2]
        else:
            return [res2,res1]

    def getOrderedAtomPair(at1,at2,useChains=False,useModels=False):
        atom1 = bnsTopLib.StructureWrapper.Atom(at1, useChains, useModels)
        atom2 = bnsTopLib.StructureWrapper.Atom(at2, useChains, useModels)
        if atom1 < atom2:
            return [atom1,atom2]
        else:
            return [atom2,atom1]
    

