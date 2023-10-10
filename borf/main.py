from typing import Iterator
import numpy as np
import argparse
import copy
import sys
import os

from classes.txgroup import Transcriptome, Gene, Bundle
from classes.transcript import Transcript

from utils.common import *

def compare_chains(chain1, chain2, strand):
    res = {
        "ilpi":0,
        "mlpi":0,
        "lpi":0,
        "num_bp_extra":0,
        "num_bp_missing":0,
        "num_bp_inframe":0,
        "num_bp_match":0,
        "num_bp_outframe":0,
        "match_start":False,
        "match_end":False
    }
    if chain2 is np.nan or len(chain2) == 0:
        return res
    if chain1 is np.nan or len(chain1) == 0:
        return res

    # 1. compute the total number of matching positions between query and template
    # 2. compute the number of matching positions in frame between query and template
    mod_chain = compare(chain1, chain2)

    c1len = clen(chain1)
    c2len = clen(chain2)

    if strand == "-":
        mod_chain.reverse()

    t_frame = 0
    q_frame = 0

    for mc in mod_chain:
        if (mc[2] == -1):  # extra positions in the query
            res["num_bp_extra"] += slen(mc)
            q_frame += slen(mc)
        elif (mc[2] == 1):  # template positions missing from the query
            res["num_bp_missing"] += slen(mc)
            t_frame += slen(mc)
        elif (mc[2] == 0):  # matching positions between query and template
            res["num_bp_match"] += slen(mc)
            if (q_frame % 3 == t_frame % 3):
                res["num_bp_inframe"] += slen(mc)  # TODO: shouldn't this be stranded?
            else:
                res["num_bp_outframe"] += slen(mc)
        else:
            print("wrong code")
            return

    # compute lpi, ilpi, mlpi, etc
    res["lpi"] = int((100.0 * (float(c1len) / float(c2len))))
    res["ilpi"] = int((100.0 * (float(res["num_bp_inframe"]) / float(c2len))))
    res["mlpi"] = int((100.0 * (float(res["num_bp_match"]) / float(c2len))))

    match_start = chain1[0][0] == chain2[0][0] if strand == '+' else chain1[-1][1] == chain2[-1][1]
    match_end = chain1[-1][1] == chain2[-1][1] if strand == '+' else chain1[0][0] == chain2[0][0]

    return res

def borf(args):

    out_fname = args.output.rstrip(".gtf").rstrip(".GTF").rstrip(".gff").rstrip(".GFF")
    out_fp = open(out_fname,"w+") # clustered version of the annotation file
    out_fp.write("locus\tall_transcripts\trepresentative\n")

    transcriptome = Transcriptome()
    transcriptome.build_from_file(args.gtf)
    if args.use_geneid:
        transcriptome.gid_sort()
    else:
        transcriptome.coordinate_sort()

    # logic:
    # 1. extract fasta sequences and find those that are complete ORFs
    # 2. cluster ORFs using ILPI between all ORFs (complete or incomplete)
    # 3. select representative ORF for each group this time making sure only complete ORFs are considered
    # 4. TODO: what to do if too few proteins (e.g. 2) are found in a cluster?

    for locus in transcriptome.gene_it() if args.use_geneid else transcriptome.bundle_it():
        locus_name = locus.get_gid() if args.use_geneid else locus.seqid+locus.strand+str(locus.start)+"-"+str(locus.end)
        
        strand = locus.strand
        coding_transcripts = [x for x in locus.object_it() if x.has_cds()]
        if len(coding_transcripts) == 0:
            continue
        elif len(coding_transcripts) == 1:
            out_fp.write(locus_name+"\t"+coding_transcripts[0].get_tid()+"\t"+coding_transcripts[0].get_tid()+"\n")
        else:
            mat_idxs = [] # transcript IDs that correspond to the rows/columns in the matrix
            dist_mat = np.zeros((len(coding_transcripts), len(coding_transcripts)))
            
            for i1,tx1 in enumerate(coding_transcripts):
                if not tx1.has_cds():
                    continue
                cds1 = sorted(tx1.cds)
                mat_idxs.append(tx1.get_tid())
                for i2,tx2 in enumerate(coding_transcripts):
                    if not tx2.has_cds():
                        continue
                    cds2 = sorted(tx2.cds)
                    comparison = compare_chains(cds1, cds2, strand)
                    dist_mat[i1,i2] = comparison["ilpi"] # TODO: could make it consider other parameters in addition to ILPI

            # find the best representative transcript
            # compute median of medians
            moms = np.median(np.sort(dist_mat, axis=1), axis=1)
            # get index of the best element
            rep_idx = np.argmin(moms)

            all_txs = ",".join([x.get_tid() for x in coding_transcripts])

            out_fp.write(locus_name+"\t"+all_txs+"\t"+coding_transcripts[rep_idx].get_tid()+"\n")


    out_fp.close()

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("-g",
                        "--gtf",
                        required=True,
                        type=str,
                        help="Annotation in a GFF/GTF format.")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        type=str,
                        help="Output base file name. Anythign after the last dot will be ignored.")
    parser.add_argument("--use_geneid",
                        action="store_true",
                        required=False,
                        help="If selected will use gene_id attribute to group transcripts. \
                            Without this flag transcripts are combined based on overlap instead.")
    
    parser.set_defaults(func=borf)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
