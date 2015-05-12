#!/usr/local/bin/python
"""
    Script: Prokka Genbank Parser
    Author: Jordan Ramsdell
    Desc:
        Used to parse through genbanks created by Prokka. Returns counts
        based on user-supplied regular expressions. Prints to standard out.

"""

import sys
import re
import mmap
from collections import Counter

def is_near_ends(genbank, pattern, size):
    counter = 0
    instances = 0
    pat_size = re.compile('(\d+)\.\.(\d+)')
    iterator = re.finditer('product=\"(.*?)\"', genbank, re.DOTALL)

    for m in iterator:
        if not re.search(pattern, m.group(1)):
            continue
        instances += 1
        start = m.start()

        m2 = pat_size.search(genbank, start - 400, start)
        try:
            m2.group(1)
        except:
            continue

        if int(m2.group(1)) <= 100 or abs(size - int(m2.group(2))) <= 100:
            counter += 1

    return counter, instances


#--------------------------- main --------------------------------------
def main():
    args = sys.argv

    # read genbank file
    f = open(args[1], "r")
    file1_lines = f.readlines()
    f.close()

    mydict1 = {}
    cur_locus = None
    for line in file1_lines:
        if "LOCUS" in line:
            cur_locus = line
            mydict1[line] = bytearray()
        elif cur_locus and not "note=" in line:
            mydict1[cur_locus].extend(line)

    # open regex file to parse through
    f = open(args[2], "r")
    reg_list = f.readlines()
    f.close()

    headers = ["regex", "count", "smaller_than_5000", "at ends", "average size"]
    print "\t".join(headers)

    # calculate statistics per regular expression
    for regex in reg_list:

        count = 0
        total_count = 0
        total = 0
        whole_size = 0
        transpose_count = 0
        grase_count = 0
        ends_count = 0
        trna_count = 0
        pol_count = 0
        repeat_count = 0
        transcription_count = 0
        sigma_count = 0
        atp_count = 0

        regex = regex.rstrip()
        p = re.compile(regex)
        for i in mydict1:
            genbank = mydict1[i]

            if re.search(p, genbank):
                sz = int(re.search("(\d*)\sbp", i).group(1))
                cend, cinst = is_near_ends(genbank, p, sz)
                ends_count += cend
                total_count += cinst
                transpose_count += genbank.count("ransposase")
                grase_count += genbank.count("grase")
                trna_count += genbank.count("tRNA-")
                pol_count += genbank.count("polymerase subunit")
                repeat_count += genbank.count("repeat")
                transcription_count += genbank.count("transcription")
                sigma_count += genbank.count("sigma")
                atp_count += genbank.count("ATP-binding")

                if sz <= 5000:
                    count += 1
                total += sz

        tot = 0 if not total_count else total / total_count
        # statlist = [regex, count, tot, total_count,
        #            transpose_count, grase_count, ends_count,
        #            trna_count, pol_count, repeat_count,
        #            transcription_count, sigma_count, atp_count]

        statlist = [regex, total_count, count, ends_count, tot]
        statlist = map(lambda x: str(x), statlist)

        print "\t".join(statlist)

# execute main if this script was run via command-line
if __name__ == "__main__":
    main()
