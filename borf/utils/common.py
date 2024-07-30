# contains reusable funcitons used throughout the experiments


# 1s cds chains of two transcripts and returns a set of statistics:
# 1. number match
#    - inframe
#    - outframe
# 2. number mismatch
# 3. match start
# 4. match stop

import os
import re
import random
import pyBigWig
import subprocess
import numpy as np
import pandas as pd
from enum import Enum, auto
from intervaltree import Interval, IntervalTree

from typing import Tuple,List

def it_eq(it1,it2):
    '''
    Evaluates equality two interval trees. Only compares the intervals, and not objects stored within
    '''
    if len(it1) != len(it2):
        return False
    
    for i1,i2 in zip(sorted(it1),sorted(it2)):
        if i1[0] != i2[0] or i1[1] != i2[1]:
            return False
        
    return True

class Types (Enum):
    Transcript = 1
    MRNA = 2
    UTR = 3
    UTR5p = 4
    UTR3p = 5
    Bundle = 6
    Gene = 7
    Exon = 8
    CDS = 9
    Intron = 10
    Other = 11

    @staticmethod
    def str2type(type_str: str):

        if type_str == "transcript":
            return Types.Transcript
        elif type_str == "mRNA":
            return Types.MRNA
        elif type_str == "UTR":
            return Types.UTR
        elif type_str == "UTR5p":
            return Types.UTR5p
        elif type_str == "UTR3p":
            return Types.UTR3p
        elif type_str == "bundle":
            return Types.Bundle
        elif type_str == "gene":
            return Types.Gene
        elif type_str == "exon":
            return Types.Exon
        elif type_str == "CDS":
            return Types.CDS
        elif type_str == "intron":
            return Types.Intron
        else:
            raise ValueError("unknown type: "+type_str)
        
    @staticmethod
    def type2str(type):
        if type == Types.Transcript:
            return "transcript"
        elif type == Types.MRNA:
            return "mRNA"
        elif type == Types.UTR:
            return "UTR"
        elif type == Types.UTR5p:
            return "UTR5p"
        elif type == Types.UTR3p:
            return "UTR3p"
        elif type == Types.Bundle:
            return "bundle"
        elif type == Types.Gene:
            return "gene"
        elif type == Types.Exon:
            return "exon"
        elif type == Types.CDS:
            return "CDS"
        elif type == Types.Intron:
            return "intron"
        else:
            raise ValueError("unknown type: "+str(type))


gff3cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"] # columns in a GFF3 file

def load_fasta_dict(fa_fname: str, rev: bool=False, upper: bool=False) -> dict:
    """
    loads a fasta file into a dictionary
    Args:
        fa_fname: fasta file name
        rev: if True, returns a dictionary with sequences as keys and names as values
        upper: if True, converts all sequences to upper case

    Returns: dictionary with names as keys and sequences as values
    """
    res = dict()
    with open(fa_fname, "r") as inFP:
        cur_nm = None

        for line in inFP:
            if line[0] == ">":
                cur_nm = line.strip()[1:].split()[0]
                assert cur_nm not in res, "duplicate record name: " + cur_nm
                res[cur_nm] = ""
            else:
                assert cur_nm is not None, "empty record name"
                res[cur_nm] += line.strip().upper()

    if rev:
        im = dict()
        for k, v in res.items():
            im[v] = im.get(v, []) + [k]

        res = im
    return res

def intersect(s1: Tuple[int,...], s2: Tuple[int,...]) -> Tuple[int,Tuple[int,int,int]]:
    """
    returns the intersection of two intervals
    Args:
        s1: first interval
        s2: second interval

    Returns: length of the intersection and the intersection itself
    """
    res = [0, -1, 0]
    tis = max(s1[0], s2[0])
    tie = min(s1[1], s2[1])
    if (tis <= tie):
        res[0] = tis
        res[1] = tie
        return (tie - tis) + 1, res
    return 0, res

def split(s1: Tuple[int,...], s2: Tuple[int,...]) -> Tuple[Tuple[int,int,int],Tuple[int,int,int],Tuple[int,int,int]]:
    """
    splits the first interval into three intervals: left, intersection, right
    Args:
        s1: first interval
        s2: second interval

    Returns: left, intersection, right
    """
    left = [0, -1, -1]
    right = [0, -1, -1]

    il, inter = intersect(s1, s2)
    if il > 0:
        if inter[0] > s1[0]:
            left[0] = s1[0]
            left[1] = inter[0] - 1
            left[2] = s1[2]
        if inter[1] < s1[1]:
            right[0] = inter[1] + 1
            right[1] = s1[1]
            right[2] = s1[2]
    else:
        if s1[0] < s2[0]:
            left = s1
        else:
            right = s1

    return left, inter, right


def slen(s: Tuple[int,...]) -> int:
    """
    returns the length of an interval
    Args:
        s: interval
    Returns: length of the interval
    """
    return (s[1] - s[0]) + 1


def clen(chain: List[Tuple[int,...]]) -> int:
    """
    returns the length of a chain of intervals
    Args:
        chain: chain of intervals
    Returns: length of the chain
    """
    res = 0
    for c in chain:
        res += slen(c)
    return res


def compare(i1: List[Tuple[int,...]], i2: List[Tuple[int,...]]) -> List[Tuple[int,int,int]]:
    """
    compares two chains of intervals
    Args:
        i1: first chain of intervals
        i2: second chain of intervals

    Returns: TODO:
    """
    intervals = []
    for i in i1:
        intervals.append([i[0], i[1]])
        intervals[-1].append(-1)
    for i in i2:
        intervals.append([i[0], i[1]])
        intervals[-1].append(1)
    intervals.sort()

    if len(i1) == 0 and len(i2) == 0:
        return []

    stack = []
    stack.append(intervals[0])
    for i in intervals[1:]:

        left, inter, right = split(stack[-1], i)
        if slen(right) > 0:
            assert slen(inter) == slen(i)  # must be intirely contained within
        else:
            tmp, inter2, right = split(i, stack[-1])
            if (slen(tmp) > 0):
                t2 = stack[-1]
                stack[-1] = tmp
                stack.append(t2)

            else:
                assert slen(tmp) <= 0, str(tmp) + "," + str(inter2) + "," + str(right)
            assert inter == inter2

        stack.pop()

        if slen(left) > 0:
            stack.append(left)
        if slen(inter) > 0:
            inter[2] = 0
            stack.append(inter)
        if slen(right) > 0:
            stack.append(right)

    return stack

# runs compare() funciton and labels all matches as in and out of frame accordingly
def compare_label_frame(chain1, chain2, strand):
    if chain2 is np.nan or len(chain2) == 0:
        [[x[0], x[1], -1] for x in chain1]
    if chain1 is np.nan or len(chain1) == 0:
        [[x[0], x[1], 1] for x in chain2]

    mod_chain = compare(chain1, chain2)

    if strand == "-":
        mod_chain.reverse()

    t_frame = 0
    q_frame = 0

    for mc in mod_chain:
        if (mc[2] == -1):  # extra positions in the query
            q_frame += slen(mc)
        elif (mc[2] == 1):  # template positions missing from the query
            t_frame += slen(mc)
        elif (mc[2] == 0):  # matching positions between query and template
            if (q_frame % 3 == t_frame % 3):
                mc[2] = 100  # inframe
            else:
                mc[2] = -100  # outframe
        else:
            print("wrong code")
            return

    return mod_chain


def compare_and_extract(chain1, chain2, strand):
    if chain2 is np.nan or len(chain2) == 0:
        return pd.Series([[[x[0], x[1], -1] for x in chain1], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
    if chain1 is np.nan or len(chain1) == 0:
        return pd.Series([[[x[0], x[1], 1] for x in chain2], -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])

    # 1. compute the total number of matching positions between query and template
    # 2. compute the number of matching positions in frame between query and template
    mod_chain = compare(chain1, chain2)

    c1len = clen(chain1)
    c2len = clen(chain2)

    if strand == "-":
        mod_chain.reverse()

    num_bp_extra = 0
    num_bp_missing = 0
    num_bp_inframe = 0
    num_bp_match = 0
    num_bp_outframe = 0

    t_frame = 0
    q_frame = 0

    for mc in mod_chain:
        if (mc[2] == -1):  # extra positions in the query
            num_bp_extra += slen(mc)
            q_frame += slen(mc)
        elif (mc[2] == 1):  # template positions missing from the query
            num_bp_missing += slen(mc)
            t_frame += slen(mc)
        elif (mc[2] == 0):  # matching positions between query and template
            num_bp_match += slen(mc)
            if (q_frame % 3 == t_frame % 3):
                num_bp_inframe += slen(mc)  # TODO: shouldn't this be stranded?
            else:
                num_bp_outframe += slen(mc)
        else:
            print("wrong code")
            return

    # compute lpi, ilpi, mlpi, etc
    lpi = int((100.0 * (float(c1len) / float(c2len))))
    ilpi = int((100.0 * (float(num_bp_inframe) / float(c2len))))
    mlpi = int((100.0 * (float(num_bp_match) / float(c2len))))

    match_start = chain1[0][0] == chain2[0][0] if strand == '+' else chain1[-1][1] == chain2[-1][1]
    match_end = chain1[-1][1] == chain2[-1][1] if strand == '+' else chain1[0][0] == chain2[0][0]

    return pd.Series(
        [mod_chain, c1len, c2len, match_start, match_end, num_bp_extra, num_bp_missing, num_bp_inframe, num_bp_match,
         num_bp_outframe, lpi, ilpi, mlpi])


def load_tid2aa(fname):
    tid2aa = dict()

    with open(fname, "r") as inFP:
        cur_tid = ""
        cur_aa = ""
        for line in inFP:
            if line[0] == ">":

                if not len(cur_tid) == 0:
                    tid2aa[cur_tid] = cur_aa

                cur_tid = line[1:].rstrip()
                cur_aa = ""
            else:
                cur_aa += line.rstrip()

        if not len(cur_tid) == 0:
            tid2aa[cur_tid] = cur_aa

    res = pd.DataFrame.from_dict(tid2aa, orient="index").reset_index()
    res.columns = ["tid", "aa"]
    return res


def merge(segs):
    segs.sort()
    res = [[segs[0][0], segs[0][1]]]
    for s in segs:
        prev = res[-1]
        if s[0] <= prev[1]:
            prev[1] = max(prev[1], s[1])
        else:
            res.append([s[0], s[1]])

    return res


def load_segments(fname, feature_type, strandless):
    res = dict({"+": dict(),
                "-": dict()})
    if strandless:
        res = dict()

    with open(fname, "r") as inFP:
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            if not lcs[2] == feature_type:
                continue

            if strandless:
                res.setdefault(lcs[0], set())
                res[lcs[0]].add((int(lcs[3]), int(lcs[4])))
            else:
                res[lcs[6]].setdefault(lcs[0], set())
                res[lcs[6]][lcs[0]].add((int(lcs[3]), int(lcs[4])))

    for k, v in res.items():
        if strandless:
            res[k] = merge(list(v))
        else:
            for k2, v2 in v.items():
                res[k][k2] = merge(list(v2))

    return res


def extract_from_comp(segs):  # separated "shared,left,right" into separate objects
    left = []
    shared = []
    right = []
    for s in segs:
        if s[2] == -1:
            left.append(s[:2])
        if s[2] == 1:
            right.append(s[:2])
        if s[2] == 0:
            shared.append(s[:2])

    return left, shared, right

# extract scores
def get_scores(scores_fname, ud, num_random=None):
    res = dict()
    min_score = 0
    max_score = 0

    bw = pyBigWig.open(scores_fname)

    for k, v in ud.items():
        res.setdefault(k, [])
        for seqid, v2 in v.items():
            if not seqid in bw.chroms():
                continue
            for r in v2:
                res[k].extend(bw.values(seqid, r[0], r[1] + 1))
        if len(res[k]) > 0:
            min_score = min(min_score, min(res[k]))
            max_score = max(max_score, max(res[k]))

    if num_random is not None:
        for k, v in res.items():
            if len(v) > num_random:
                res[k] = random.sample(v, num_random)

    return res, min_score, max_score
    
# get_mean_score
def mean_score(score_fname):
    means = []
    bw = pyBigWig.open(score_fname)
    for seqid, l in bw.chroms().items():
        means.append(bw.stats(seqid, 0, l - 1))
    return np.mean(means)


# extract sashimi and gtf for a specified set of transcripts based on several annotations
def extract_sashimi(sbin, cmp_gtf_fname, ref_gtf_fname, q_gtf_fname, out_base_fname, cmp_tid, qtids, title_str):
    out_gtf_fname = out_base_fname + ".gtf"
    out_svg_fname = out_base_fname + ".svg"

    with open(out_gtf_fname, "w+") as outFP:
        # first write out MANE
        with open(cmp_gtf_fname, "r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid == cmp_tid:
                    outFP.write(line)
        # next write out ORFanage
        with open(q_gtf_fname, "r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid in qtids or tid == qtids:
                    lcs[8] = "transcript_id \"ORFanage:" + tid + "\""
                    line = "\t".join(lcs) + "\n"
                    outFP.write(line)
        # lastly write out regular RefSeq
        with open(ref_gtf_fname, "r") as inFP:
            for line in inFP:
                lcs = line.split("\t")
                if not len(lcs) == 9:
                    continue
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid in qtids or tid == qtids:
                    outFP.write(line)

    sashimi_cmd = [sbin,
                   "--compare", cmp_tid,
                   "--title", title_str,
                   "--gtf", out_gtf_fname,
                   "-o", out_svg_fname]
    print(" ".join(sashimi_cmd))
    subprocess.call(sashimi_cmd)


def subset_gtf_by_seqid(in_gtf_fname, out_gtf_fname, seqids):
    with open(out_gtf_fname, "w+") as outFP:
        with open(in_gtf_fname, "r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                if not seqids == False and lcs[0] not in seqids:
                    continue

                outFP.write(line)


def subset_gtf(in_gtf_fname, out_gtf_fname, gids, tids):
    writing_tid = ""
    with open(out_gtf_fname, "w+") as outFP:
        with open(in_gtf_fname, "r") as inFP:
            for line in inFP:
                lcs = line.rstrip().split("\t")
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if lcs[2] == "transcript":
                    gid = lcs[8].split("gene_id \"", 1)[1].split("\"", 1)[0]
                    if not gids == False and gid in gids:
                        outFP.write(line)
                        writing_tid = tid
                        continue

                if not tids == False and tid in tids:
                    outFP.write(line)
                    continue

                # handle non transcript ffeatures without gene_id for whihc a transcript was found based on gene-id
                if writing_tid == tid:
                    outFP.write(line)
                    continue

# extracts num_pos coordinates from the chain
# if reverse - will extract from the end
def get_coords(chain, num_pos, reverse):
    res_coords = []

    tmp_chain = chain
    inc = 1
    if reverse:
        inc = -1
        tmp_chain = [[x[1], x[0]] for x in tmp_chain[::-1]]

    for c in tmp_chain:
        for i in range(c[0], c[1] + 1, inc):
            res_coords.append(i)
            if len(res_coords) >= num_pos:
                return res_coords


def contained_intervals(is1, is2, inverse=False):  # if inverse - return intervals of is2 hich contain intervals in is1
    res = []
    for i1 in is1:
        for i2 in is2:
            if i1[0] >= i2[0] and i1[1] <= i2[1]:
                if inverse:
                    res.append(i2)
                else:
                    res.append(i1)
                break
    return res

# extract attribute key values into dictionary
def extract_attributes(attribute_str: str) -> dict:
    gff = attribute_str.startswith("ID=") or attribute_str.startswith("Parent=")

    attrs = attribute_str.rstrip().rstrip(";").split(";")
    attrs = [x.strip() for x in attrs]
    attrs = [x.strip("\"") for x in attrs]
    attrs_dict = dict()
    sep = " \""
    if gff:
        sep = "="
    for at in attrs:
        k, v = at.split(sep)
        attrs_dict.setdefault(k, v)

    return attrs_dict

# renames attributes
def rename_attributes(attrs: dict, rename_dict=dict) -> dict:
    res_dict = {}
    for k, v in attrs.items():
        if k in rename_dict:
            res_dict[rename_dict[k]] = v
        else:
            res_dict[k] = v
    return res_dict

# converts attribute key values back into string
def to_attribute_string(attrs: dict, gff=False,
                        feature_type=None) -> str:
    order = ["ID", "Parent", "transcript_id", "gene_id", "gene_name", "gene_type", "db_xref", "description", "max_TPM",
            "sample_count", "assembly_id", "tag"]
    res = ""
    sep = " "
    quote = "\""
    end = "; "
    if gff:
        assert feature_type in ["gene", "transcript", "exon", "CDS"], "wrong type: " + str(feature_type)
        sep = "="
        quote = ""
        end = ";"

    for k in order:
        if k in attrs:
            if gff:
                assert ";" not in attrs[k], "invalid character in attribute: " + attrs[k]

            if gff and feature_type == "gene" and k == "transcript_id":
                continue
            elif gff and feature_type == "gene" and k == "gene_id":
                res += "ID=" + quote + attrs[k] + quote + end
            elif gff and feature_type == "transcript" and k == "transcript_id":
                res += "ID=" + quote + attrs[k] + quote + end
            elif gff and feature_type == "transcript" and k == "gene_id":
                res += "Parent=" + quote + attrs[k] + quote + end
            elif gff and feature_type in ["exon", "CDS"] and k == "transcript_id":
                res += "Parent=" + quote + attrs[k] + quote + end
            elif gff and feature_type in ["exon", "CDS"] and k == "gene_id":
                continue
            else:
                res += k + sep + quote + attrs[k] + quote + end

    # add any other attributes in sorted order
    for k in sorted(list(attrs)):
        if k not in order:
            if gff:
                assert ";" not in attrs[k], "invalid character in attribute: " + attrs[k]
            res += k + sep + quote + attrs[k] + quote + end

    if not gff:
        res = res.rstrip()
    if gff:
        res = res.rstrip(";")
    return res

def get_intervals(gtf_fname, feature="exon", invert=False):
    res_intervals = dict()  # seqids and strand as keys, lists of introns as values1

    intervals = {}
    with open(gtf_fname, 'r') as inFP:
        for line in inFP:
            if line[0] == "#":
                continue
            lcs = line.strip().split('\t')
            if lcs[2] == feature:
                tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
                if tid not in intervals:
                    intervals[tid] = {"seqname": lcs[0],
                                      "strand": lcs[6],
                                      "intervals": []}
                intervals[tid]["intervals"].append((int(lcs[3]), int(lcs[4])))

    for tid, idata in intervals.items():
        if invert:  # get introns
            for ii in range(1, len(idata["intervals"][1:]), 1):
                key = (idata["seqname"], idata["strand"])
                res_intervals.setdefault(key, dict())
                rit = (idata["intervals"][ii - 1][1] + 1, idata["intervals"][ii][0] - 1)
                res_intervals[key].setdefault(rit, set())
                res_intervals[key][rit].add(tid)
        else:
            for ii in range(len(idata["intervals"])):
                key = (idata["seqname"], idata["strand"])
                res_intervals.setdefault(key, dict())
                rit = (idata["intervals"][ii][0], idata["intervals"][ii][1])
                res_intervals[key].setdefault(rit, set())
                res_intervals[key][rit].add(tid)

    return res_intervals

# find longest ORF in each transcript
# what do we do if there is multiple ORFs of the same length? - just count for now. could also just skip those genes alltogether
def find_longest_orfs(seq):
    longest = []
    max_len = 0

    matches = re.finditer(r'(?=(ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)))', seq)
    for match in matches:
        result = match.group(1)
        coords = [match.start(), match.start() + len(result) - 1]

        if max_len < len(result):
            longest = [coords]
            max_len = len(result)
            continue
        if max_len == len(result):
            longest.append(coords)
            continue

    return longest


def trans2genome(chain, strand, zero_pos):
    chain_pos = -1
    left_to_stop = zero_pos
    found_pos = False
    if strand == '+':
        for i in range(len(chain)):
            clen = slen(chain[i])
            if left_to_stop < clen:  # found the segment with the stop codon
                chain_pos = chain[i][0] + left_to_stop
                found_pos = True
                break

            left_to_stop -= clen

        if not found_pos:  # return the last position
            chain_pos = chain[-1][1]

    else:
        for i in range(len(chain) - 1, -1, -1):
            clen = slen(chain[i])
            if left_to_stop < clen:  # found the cds segment with the stop codon
                chain_pos = chain[i][1] - left_to_stop
                found_pos = True
                break

            left_to_stop -= clen

        if not found_pos:  # return the last position
            chain_pos = chain[0][0]

    assert chain_pos >= 0, "unexpected chain_pos<0"
    return chain_pos


def cut_chain(chain, start, end):
    res = []
    for cs, ce in chain:
        new_cs = cs
        new_ce = ce
        if new_cs <= start and new_ce >= start:
            new_cs = start
        if new_ce >= end:
            new_ce = end
            res.append([new_cs, new_ce])
            break
        if new_ce < start or new_cs > end:
            continue
        res.append([new_cs, new_ce])
    return res
