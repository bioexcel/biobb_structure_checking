
import re

class ResidueSetList():
    def __init__(self, pairsList=None, debug=False):
        self.sets=[]
        self.n=0
        if (pairsList):
            self.addData(pairsList,debug)
            self.addMeta()

    def addMeta(self):
        chnum=1
        for s in self.getSortedSets():
            s.setMeta(chnum)
            chnum=chnum+1
        
    def addData(self, pairsList, debug=False):
        for [res1,res2] in pairsList:        
            i = self.find(res1)
            j = self.find(res2)
            if i == -1 and j == -1:
                s = ResidueSet()
                s.add(res1)
                s.add(res2)
                self.append(s)
            elif i != -1 and j != -1 and i != j:
                self.sets[i].union(self.sets[j])
                self.delete(j)
            elif j == -1:
                self.sets[i].add(res2)
            elif i == -1:
                self.sets[j].add(res1)
            if debug:
                for s in self.getSortedSets():
                    print ("#DEBUG:" ,s)
    
    def find(self,item):
        i = 0
        while i < self.n and item not in self.sets[i].items:
            i = i + 1
        if i == self.n:
            return -1
        else:
            return i

    def append(self,item):
        self.sets.append(item)
        self.n = len(self.sets)

    def delete(self,i):
        del self.sets[i]
        self.n = len(self.sets)

    def getSortedSets(self):
        return sorted(self.sets,key=lambda s: s.ini)
    

class ResidueSet():
    def __init__(self):
        self.ini=999999
        self.fin=0
        self.iniId=''
        self.finId=''
        self.items = set()
        self.type=''
        self.id=''
    
    def setMeta(self,id):
        if self.isProtein():
            self.type='prot'
        else:
            self.type='na'
        self.id=id
        
    def add(self,r):
        self.items.add(r)
        rn=r.__index__()
        if rn < self.ini:
            self.ini = rn
            self.inir = r.resid()
        if rn > self.fin:
            self.fin = rn
            self.finr = r.resid()
        
    def union(self,other):
        self.items = self.items.union(other.items)
        if other.ini < self.ini:
            self.ini  = other.ini
            self.inir = other.inir
        if other.fin > self.fin:
            self.fin = other.fin
            self.finr = other.finr

    def getSequence(self):
        seq=self._getResidues()
        ss=''
        for i in sorted(seq.keys()):
            ss=ss+seq[i]._getOneLetterResidueCode()
        return ss

    def getResidueIdList(self):
        seq=self._getResidues()
        seql = []
        for i in sorted(seq.keys()):
            rid = (seq[i].chain+str(seq[i].resNum)+"-"+seq[i]._getOneLetterResidueCode())
            if seq[i].useModels:
                seql.append(rid + "/"+ str(seq[i].model))
            else:
                seql.append(rid)
        return seql
    
    def isProtein(self):
        seq = self.getSequence()
        seq = re.sub('[ACTGUX]','',seq)
        if seq:
            return True
        else:
            return False

    def _getResidues(self):
        seq={}
        for r in self.items:
            seq[r.__index__()] = r
        #print (seq)
        return seq
    
    def __str__(self, mode=0):
        
        if mode == 0:
            return '%(id)s %(inir)s-%(finr)s (%(type)s):' % vars(self) + self.getSequence()
        elif mode == 1:
            return '%(id)s (%(inir)s-%(finr)s):' % vars(self) + ','.join(self.getResidueIdList())
        else:
            return ''
    
    def json(self):
        return  {
            'id': self.id,
            'iniRes' : str(self.inir), 
            'finRes' : str(self.finr), 
            'sequence' : self.getSequence(),
            'protein': self.isProtein(),
            'residueList' : self.getResidueIdList()
        }
    

class BPSSetList (ResidueSetList):
    def addData(self, pairsList, debug=False):
        for [res1,res2] in pairsList:        
            i = self.find(res1)
            j = self.find(res2)
            if i == -1 and j == -1:
                s = BPSSet()
                s.add(res1)
                s.add(res2)
                self.append(s)
            elif i != -1 and j != -1 and i != j:
                self.sets[i].union(self.sets[j])
                self.delete(j)
            elif j == -1:
                self.sets[i].add(res2)
            elif i == -1:
                self.sets[j].add(res1)
            if debug:
                for s in self.getSortedSets():
                    print ("#DEBUG:" ,s)
    
class BPSSet (ResidueSet):
    def add(self,bpst):
        self.items.add(bpst)
        rn=bpst.bp1.r1.residue.index
        if rn < self.ini:
            self.ini = rn
            self.inir = bpst
        if rn > self.fin:
            self.fin = rn
            self.finr = bpst
    
    def getSequence(self):
        seq=self._getResidues()
        ss=''
        for i in sorted(seq.keys()):
            ss=ss+seq[i].bp1.bpid()+","
        return ss
    
    def _getResidues(self):
        seq={}
        for bpst in self.items:
            seq[bpst.bp1.r1.residue.index] = bpst
        return seq
    
    def __str__(self):
        return str(self.inir.bp1.r1.resid(1)) + "-" + str(self.finr.bp1.r1.resid(1))+ ":  " + self.getSequence()
    
