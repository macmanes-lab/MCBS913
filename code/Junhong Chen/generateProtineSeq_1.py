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
        
        self.name = "name"
        self.data = dict()
        self.fname = gff
        self.reverse = 0
        
    def parse(self):
        
        file = open(self.fname,"r")
        for elem in file:
            if "CDS" in elem:
                tmp = elem.split()
                ind = tmp.index("CDS")
                if self.name in self.data:
                    self.data[self.name].append((int(tmp[ind+1]),int(tmp[ind+2]),tmp[ind+4]))
                else:
                    self.data[self.name] = [(int(tmp[ind+1]),int(tmp[ind+2]),tmp[ind+4])]
                    
                if(tmp[ind+4] == '+'):
                    self.reverse = self.reverse + 1
                    
        self.length = len(self.data[self.name])
                
                    
                
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
        
        self.name = "name"
        self.fname = fast
        self.data = dict()
        self.result = dict()
        self.cds = CDS(fast.split(".")[0]+".gff")
        
        
    def parse(self):
        
        fast = SeqIO.parse(open(self.fname),"fasta")
        
        i = 0
        for elem in fast:
            i = i + 1
            tmp = elem.id.split("|")[3]
            if tmp in self.data:
                print "ATTENTION: same contig ID in: " + self.fname
            else:
                self.data[self.name] = str(elem.seq)
                
        #print "Seq in Reference: ",i
                
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
        #print self.cds.getContigName();
        return sq[ind[0]-1:ind[1]]
    
    
    def getCDS(self):
        return self.cds;



def compareProtineSeq(path):
    
    refd = RefSeq(path+".fastd")
    refd.parse()

    refa = RefSeq(path+".fasta")
    refa.parse()

    
    refat = refa.translate()
    refdt = refd.translate()
    
    #print "Match: CDS Index [71,1498], Annotation +"
    #print "----------------           Normal reference           ----------------------------"
    #print "Unambiguous CDS: "
    #print
    #print refa.getCDSSeq("name",0)
    #print
    #print "Tranlated protein sequence: "
    #print
    #print refat["name"][0]
    
    #print
    #print
    #print "----------------           Degenerate reference           -------------------------"
    #print "Ambiguous CDS"
    #print
    #print refd.getCDSSeq("name",0)
    #print
    #print "Tranlated protein sequence: "
    #print
    #print refdt["name"][0]
    
    
    #print
    #print
    #print "Mismatch: CDS Index [5603,7309], Annotation -"
    #print "----------------           Normal reference           ----------------------------"
    #print "Unambiguous CDS: "
    #print
    #print refa.getCDSSeq("name",3)
    #print
    #print "Reverse_complement of CDS"
    #print
    #print Seq(refa.getCDSSeq("name",3)).reverse_complement()
    #print
    #print "Tranlated protein sequence: "
    #print
    #print refat["name"][3]
    
    #print
    #print
    #print "----------------           Degenerate reference           -------------------------"
    #print "Ambiguous CDS"
    #print
    #print refd.getCDSSeq("name",3)
    #print
    #print "Reverse_complement of CDS"
    #print
    #print Seq(refd.getCDSSeq("name",3)).reverse_complement()
    #print
    #print "Tranlated protein sequence: "
    #print
    #print refdt["name"][3]    
    
    
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
    
    x,y = doCompare(path)
    
    for i in range(len(x)):
        print x[i]," ",y[i]
    
    #refa = RefSeq(path+".fasta")
    #refa.parse()
    #print refa.getCDSSeq("NC_008752",0)
