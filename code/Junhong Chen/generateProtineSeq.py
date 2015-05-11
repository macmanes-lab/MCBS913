"""
Author: Junhong Chen

"""

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from sys import argv
import os

path = argv[1]



class CDS:
    def __init__(self,gff):
        
        self.data = dict()
        self.fname = gff
        
    def parse(self):
        
        file = open(self.fname,"r")
        for elem in file:
            if "CDS" in elem:
                tmp = elem.split()
                ind = tmp.index("CDS")
                if tmp[0] in self.data:
                    self.data[tmp[0]].append((int(tmp[ind+1]),int(tmp[ind+2]),tmp[ind+4]))
                else:
                    self.data[tmp[0]] = [(int(tmp[ind+1]),int(tmp[ind+2]),tmp[ind+4])]
                    
                
    def getContigName(self):
        
        return self.data.keys()
    
    
    def getContigNumber(self):
        
        return len(self.data)
    
    
    def getContigCDSIndex(self,name):
        
        if name in self.data:
            return self.data[name]
        else:
            print "No indices for that contig ID: ", name
            #return self.data[name.split(".")[0]]    
    
    
    def getContigCDSSize(self,name):
        
        return len(self.getContigCDSIndex(name))




class RefSeq:
    
    def __init__(self,fast):
        
        self.fname = fast
        self.data = dict()
        self.result = dict()
        self.cds = CDS(fast.split(".")[0]+".gff")
        
        
    def parse(self):
        
        fast = SeqIO.parse(open(self.fname),"fasta")
        
        for elem in fast:
            
            tmp = elem.id.split("|")[3]
            if tmp in self.data:
                print "ATTENTION: same contig ID in: " + self.fname
            else:
                self.data[tmp] = str(elem.seq)
                
    def getContigSeq(self,name):
        
        if name in self.data:
            return self.data[name]
        else:
            print "Can NOT find the contig: "+name
            
            
    def getContigData(self):
        return self.data
    
            
    def getContigID(self):
        
        return self.data.keys()
    
    def getContigCDSSize(self,name):
        
        return self.cds.getContigCDSSize(name)

        
    def translate(self,mode = IUPAC.ambiguous_dna):
        
        self.cds.parse()
        contig = self.data.keys()
        
        for name in contig:
            ind = self.cds.getContigCDSIndex(name)
            sq = self.data[name]
            ret = []
            for tup in ind:
                myseq = sq[tup[0]-1:tup[1]]
                
                #store Seq Object
                if tup[2] == "+":
                    ret.append(Seq(myseq, mode).translate())
                else:
                    ret.append(Seq(myseq, mode).reverse_complement().translate())
            self.result[name] = ret
            
        return self.result


    def getCDSSeq(self,name,index):
        sq = self.data[name]
        ind = self.cds.getContigCDSIndex(name)[index]
        print self.cds.getContigName();
        return sq[ind[0]-1:ind[1]]



def compareProtineSeq(path):
    
    refd = RefSeq(path+".fastd")
    refd.parse()

    refa = RefSeq(path+".fasta")
    refa.parse()
    
    
    refat = refa.translate()
    refdt = refd.translate()
    
    
    #print refat["NC_008752.1"][3]
    #print refdt["NC_008752.1"][3]
    
    #print refa.getCDSSeq("NC_008752.1",3)
    #print refd.getCDSSeq("NC_008752.1",3)
    
    
    id = refd.getContigID()
    
    ret = dict()
    
    
    for name in id:
        mis = []
        l = refa.getContigCDSSize(name)
        stat = 0
        
        for i in range(l):
            if refat[name][i] in refdt[name][i]:
                stat = stat + 1
            else:
                mis.append(i)
        
        ret[name] = (l,stat,mis)
        
    
    def sum(x):
        ret = 0.
        for el in x:
            ret = ret + el*1.
            
        return ret
    
    mis = [x[1] for x in ret.values()]
    tot = [x[0] for x in ret.values()]
    
    return sum(mis)/sum(tot)
    #return ret
    
    
 
 

def getFilesinCWD(path):

    if path[-1] is not "/":
        path = path + "/"
    
    ref = []    
    
    files = [f for f in os.listdir(path)]
    for i in range(1,5):
        for fo in files:
            f = fo.split(".")[0]
            if f not in ref and f.startswith(str(i)+"-"):
                ref.append(f)
                
    
    ret = [path+tp for tp in ref]

    return ret
 



def doCompare(path):
    
    fpath = getFilesinCWD(path)
    retp = [f.split("/")[-1] for f in fpath]
    ret = []
    for p in fpath:
        ret.append(compareProtineSeq(p))
                
    return retp,ret


if __name__ == "__main__":
    
    print doCompare(path)
    ##refa = RefSeq(path+".fasta")
    #refa.parse()
    #print refa.getCDSSeq("NC_008752",0)
