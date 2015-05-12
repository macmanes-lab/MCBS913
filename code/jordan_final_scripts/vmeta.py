#!/usr/local/bin/python
"""
    Script: Vmeta
    Author: Jordan Ramsdell
    Version: 0.3333333...
"""

import re
import cPickle as pickle
from collections import namedtuple
from collections import Counter
import mmap
import pygal


# define some globals
GO = namedtuple('GO', ['id', 'desc'])
uniprot_dict = {}
go_dict = {}
pfam_dict = {}


class Contig(object):
    """
    :ivar pfam_features :
    :type pfam_features : list of [Feature]
    :ivar uniprot_features :
    :type uniprot_features : list of [Feature]

    """

    def __init__(self, name):
        self.name = name

        # parse name for more info
        self.length = int(re.search("length_(\d*)", name).group(1))
        #self.coverage = float(re.search("cov_([1-9\.]*)", name).group(1))
        self.go_ids = []
        self.uniprot_features = []
        self.pfam_features = []

class Feature(object):
    def __init__(self):
        self.start = 0
        self.stop = 0
        self.length = 0
        self.direction = "+"
        self.term = None


class PfamTerm(object):
    def __init__(self, name):
        self.name = name
        self.terms = []
    def add_term(self, term):
        self.terms.append(term)
    def search_terms(self, pattern):
        for i in self.terms:
            if pattern.search( i.desc ):
                yield i

class UniprotTerm(object):
    def __init__(self, name):
        self.gene_names = ""
        self.keywords = ""
        self.taxa = ""
        self.comments = ""
        self.go_terms = []
        self.name = name

class GoTerm(object):
    def __init__(self, id):
        self.id = id
        self.name = ""
        self.desc = bytearray()
        self.synonym = bytearray()
        self.is_a = bytearray()
        self.namespace = bytearray()


def pfam2go(filename):
    pfam_dict = {}
    f = open(filename, "r")
    for line in f:
        if line[0] == '!':
            continue

        id = re.search("Pfam:(.*?)\s(.*?)\s", line) # 1:pfam, 2:prot_name
        term = re.search("GO:(.*?)\s\;\sGO:(\d*)", line) # 1:go_id, 2: go_desc

        try:
            cur_term = pfam_dict[id.group(1)]
        except KeyError:
            cur_term = PfamTerm(id.group(2))
            pfam_dict[id.group(1)] = cur_term

        cur_term.add_term(GO(term.group(2), term.group(1)))

    f.close()

    return pfam_dict

def create_uniprot_database(filename):
    f = open(filename, "r+b")
    uniprot_dict = {}
    term = None

    mm = mmap.mmap(f.fileno(), 0, mmap.PROT_READ)
    while 1:
        line = mm.readline()
        if line == '':
            break
        header = line[0:2]
        if header == 'AC':
            name = re.match("AC\s*(.*?)\;", line).group(1)
            term = UniprotTerm(name)
            uniprot_dict[name] = term
        elif header == 'PE' and int(line[5:6]) > 2:
            uniprot_dict.pop(term.name)
            term = None
        if not term:
            continue

        if header == 'GN':
            term.gene_names = line[5:]
        elif header == 'OC':
            term.taxa += line[5:]
        elif header == 'CC':
            term.comments += line[5:]
        elif header == 'DR' and line[5:7] == 'GO':
            m = re.search("(.*?);\s(.*?);", line[12:])
            term.go_terms.append(GO(m.group(1), m.group(2)))
        elif header == 'KW':
            term.keywords += line[5:]

    mm.close()
    f.close()

    print len(uniprot_dict)
    f = open('mydata.txt', 'w+b')
    pickle.dump(uniprot_dict, f, 2)
    f.close()

def create_godict(filename):
    f = open(filename, "r+b")
    mm = mmap.mmap(f.fileno(), 0, mmap.PROT_READ)
    term = None
    go_dict = {}
    start_reading = False
    while 1:
        line = mm.readline()
        if line == '' or line == '[Typedef]\n':
            break
        if line == '[Term]\n':
            start_reading = True

        if not start_reading or line == "\n" or line[0] == '[':
            continue

        try:
            header = re.match("(.*?):", line).group(1)
            if header == 'id':
                name = re.search("GO:(\d*)", line).group(1)
                term = GoTerm(name)
                go_dict[name] = term
            elif header == 'name':
                term.name = line[6:]
            elif header == 'namespace':
                term.namespace = line[11:]
            elif header == 'def':
                term.desc = re.search('\"(.*?)\"', line).group(1)
            elif header == 'is_a':
                term.is_a.extend(line[6:])
            elif header == 'synonym':
                term.synonym.extend(re.search('\"(.*?)\"', line).group(1))

        except:
            print line

    mm.close()
    f.close()
    return go_dict


def load_gff(filename):
    """
    :param filename:
    :type filename : str
    :rtype : dict of [str,Contig]
    """
    f = open(filename, "r+b")
    contigs = {}
    mm = mmap.mmap(f.fileno(), 0, mmap.PROT_READ)

    # read through gff
    while 1:
        line = mm.readline()
        if line == '' or line[0] == '>':
            break
        if line[0] == "#":
            continue

        elements = line.split("\t")
        contig_name = elements[0]

        if not contig_name:
            continue

        contig = contigs.get(contig_name)
        if not contig:
            contig = Contig(contig_name)
            contigs[contig_name] = contig

        # get feature info
        feature_start = elements[3]
        feature_stop = elements[4]
        feature_direction = elements[6]
        feature_notes = elements[8]

        feature = Feature()
        feature.start = feature_start
        feature.stop = feature_stop
        feature.direction = feature_direction

        match = re.search("UniProtKB:(.*?);", feature_notes)
        if match:
            term = uniprot_dict.get(match.group(1))
            if term:
                feature.term = term
                contig.uniprot_features.append(feature)

        else:
            match = re.search("Pfam:(.*?)(\.|;)", feature_notes)
            if match:
                term = pfam_dict.get(match.group(1))
                if term:
                    feature.term = term
                    contig.pfam_features.append(feature)

    mm.close()
    f.close()
    return contigs


def get_taxa(contigs, level=3):
    """
    :type contigs : dict of [Contig]
    :type level : int
    :rtype : Counter
    """
    taxa_list = []
    for key in contigs:
        contig = contigs[key]
        for i in contig.uniprot_features:
            term = i.term
            try:
                taxa_list.append(term.taxa.split(";")[level])
            except:
                continue

    return Counter(taxa_list)


def get_goterms(contigs, pattern):
    """
    :type contigs : dict of [str,Contig]
    :type pattern : str
    :rtype : Counter
    """
    global go_dict
    go_list = []
    p = re.compile(pattern)
    for key in contigs:
        contig = contigs[key]
        for i in contig.uniprot_features:
            terms = i.term.go_terms
            for term in terms:
                goterm = go_dict.get(term.id)
                if not goterm:
                    continue
                match = p.search(goterm.name)
                if match:
                    go_list.append(goterm.name.rstrip())

        for i in contig.pfam_features:
            terms = i.term.terms
            for term in terms:
                goterm = go_dict.get(term.id)
                if not goterm:
                    continue
                match = p.search(goterm.name)
                if match:
                    go_list.append(goterm.name.rstrip())

    return Counter(go_list)

def compare_differences(counter1, counter2):
    """
    :param counter1:
    :param counter2:
    :type counter1 : Counter
    :type counter2 : Counter
    :return:
    """
    dif = {}
    for key in counter1:
        dif[key] = counter1[key]

    for key in counter2:
        if dif.get(key):
            dif[key] -= counter2[key]
        else:
            dif[key] = -1 * counter2[key]

    return dif

def compare_terms(counter1, counter2, outname):
    """
    :param counter1:
    :param counter2:
    :type counter1 : Counter
    :type counter2 : Counter
    """
    bar_chart = pygal.Bar()
    total_1 = reduce(lambda x,y:x+y, [counter1[z] for z in counter1 if
                                      counter2.get(z)])
    total_2 = reduce(lambda x,y:x+y, [counter2[z] for z in counter2 if
                                      counter1.get(z)])
    for i in counter1:
        if counter2.get(i):
            bar_chart.add( i, [float(counter1[i]) / total_1, float(counter2[i]) / total_2])
    bar_chart.render_to_file(outname + ".svg")


def make_piechart(counter, outname):
    """
    :type counter: Counter
    :type outname : str
    """

    pie = pygal.Pie()
    for i in counter:
        pie.add(i, [{'value': counter[i], 'label':i}])
    pie.render_to_file(outname + ".svg")


def initialize():
    global uniprot_dict, go_dict, pfam_dict
    print "Loading databases."
    with open('mydata.txt', 'r+b') as f:
        uniprot_dict = pickle.load(f)
    go_dict = create_godict('go.obo')
    pfam_dict = pfam2go('pfam2go')




#--------------------------- main --------------------------------------
def main():

    go_dict = create_godict(options.go)
    f = open('mydata.txt', 'r+b')
    uniprot_dict = pickle.load(f)
    f.close()
    pfam_dict = pfam2go(options.pfam)
    load_gff(options.gff, uniprot_dict, pfam_dict, go_dict)






#execute main if this script was run via command-line
if __name__ == "__main__":
    main()
