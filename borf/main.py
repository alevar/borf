from typing import Iterator
import numpy as np
import argparse
import copy
import sys
import os

from classes.txgroup import Transcriptome, Gene, Bundle
from classes.transcript import Transcript

def borf(args):

    out_fname = args.output.rstrip(".gtf").rstrip(".GTF").rstrip(".gff").rstrip(".GFF")
    out_gtf_fp = open(out_fname,"w+") # clustered version of the annotation file

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

    for locus in transcriptome.gene_it() if args.use_geneid else transcriptome.bundle_it():
        for tx in locus.object_it():
            print(tx)
            break
        break

    out_gtf_fp.close()

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
