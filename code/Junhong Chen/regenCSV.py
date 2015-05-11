"""
Author: Junhong Chen

"""
import sys
import os


filename = sys.argv[1]
rootdir = sys.argv[2]


data = dict()


f = open(filename,"r")
next(f)

for line in f:
    
    tmp = line.strip()
    dl = tmp.split(",")
    head = dl[0]
    tail = dl[1::]
    name = head.split(".")[0]
    old_name = name
    if name.endswith("fd"):
        name = name[0:len(name)-3]
        tag = "degenerated reference"
    else:
        tag = "normal reference"
         
    if name not in data:
        data[name] = dict()
    
    data[name][tag] = tail
    
f.close()

ret = []
    
for dirName, subdirList, fileList in os.walk(rootdir):    
    for fname in fileList:
        if fname.endswith(".fna"):
            head = fname.split(".")[0]
            tmp = dirName + "/" + head
            tp = tmp.split("/")[6::]
            if head in data:
                res1 = tp + ["normal reference"]+data[head]["normal reference"]
                str1 = ",".join(res1)
                res2 = tp + ["degenerated reference"]+data[head]["degenerated reference"]
                str2 = ",".join(res2)
                ret.append(str1)
                ret.append(str2)
            
            
            
    


f = open(sys.argv[3],"w")

for line in ret:
    f.write(line+os.linesep)
    
f.close()
