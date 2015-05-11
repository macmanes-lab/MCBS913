"""
Author: Junhong Chen

"""

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from sys import argv
import os

path = argv[1]






class RefSeq:
    
    def __init__(self,fast):
        
        self.fname = fast
        self.gname = fast.split(".")[0]+".gff"
        self.id = None
        self.seq = None
        self.cdsseq = []
        self.cdsinfo = []
        self.reverse = None
        
        
    def parse(self):
        
        fast = SeqIO.parse(open(self.fname),"fasta")       
        for elem in fast:
            self.id = elem.id.split("|")[3]
            self.seq = str(elem.seq)        
        
        i = 0;
        file = open(self.gname,"r")
        for elem in file:
            if "CDS" in elem:
                tmp = elem.split()
                ind = tmp.index("CDS")
                self.cdsinfo.append((int(tmp[ind+1]),int(tmp[ind+2]),tmp[ind+4]))
                self.cdsseq.append((self.seq[int(tmp[ind+1])-1:int(tmp[ind+2])],tmp[ind+4]))
                i = i  + 1

        file.close()
        
        
    
    def getCDSSize(self):
        return len(self.cdsinfo)
        


    def translate(self,mode = IUPAC.ambiguous_dna):
        
        ret = []
        self.reverse = 0
        for i in range(self.getCDSSize()):
            if self.cdsseq[i][1] == '+':
                ret.append((Seq(self.cdsseq[i][0], mode).translate(),'+'))
            else:
                ret.append((Seq(self.cdsseq[i][0], mode).reverse_complement().translate(),'-'))
                self.reverse = self.reverse + 1
        
        return ret
    
    
    def getCDSSeq(self):    
        return self.cdsseq
    
    
    def getReverseNumber(self):
        return self.reverse







def doCompare(fpath):
    
    refd = RefSeq(fpath+".fastd")
    refd.parse()

    refa = RefSeq(fpath+".fasta")
    refa.parse()
    
    refat = refa.translate()
    refdt = refd.translate()
    
    length = refa.getCDSSize()
    ret = []
    ssind = []
    xind = []
    #minusn = 0
    
    for i in range(length):
        if refat[i][0] in refdt[i][0]:
            
            tmp = str(refat[i][0])
            ss = None
            nox = None
            
            #if refat[i][1] == '-':
                #minusn = miunsn + 1
            
            if(len(tmp) > 0 and tmp[0] == 'M' and tmp[-1] == '*'):
                ss = True
            else:
                ss = False
                ssind.append(i)
                
            if 'X' in tmp:
                nox = False
                xind.append(i)
            else:
                nox = True
                
            ret.append([refat[i],ss,nox])
            
          
    return len(ssind)*1.0/len(ret),len(xind)*1.0/len(ret)#,minusn*1./len(ret)




def doItAll(path):
    
    fpath = getFilesinCWD(path)

    retp = [f.split("/")[-1] for f in fpath]
    ret = []
    for p in fpath:
        ret.append(doCompare(p))
                
    return retp,ret



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



if __name__ == "__main__":
    
    print doItAll(path)
    
