"""Perform assembly based on debruijn graph."""

#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

# ===============
# IMPORT *
# ===============
# [A]
import argparse
# [C]
import csv
# [O]
import os
# [R]
import re
# [S]
import sys


__author__ = "ROUAUD Lucas"
__credits__ = __author__
__version__ = "1.0.0"
__maintainer__ = __author__
__email__ = "lucas.rouaud@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = f"{path} is a directory"
        else:
            msg = f"{path} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = f"{path} is a file"
        else:
            msg = f"{path} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__,
                                     usage=f"{sys.argv[0]} -h")
    parser.add_argument("-i", dest="genome_file", type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument("-g", dest="min_gene_len", type=int,
                        default=50, help="Minimum gene length to consider")
    parser.add_argument("-s", dest="max_shine_dalgarno_distance", type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument("-d", dest="min_gap", type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument("-p", dest="predicted_genes_file", type=str,
                        default=os.curdir + os.sep + "predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument("-o", dest="fasta_file", type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")

    return vars(parser.parse_args())


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string.
    """
    return fasta_file


def find_start(start_regex, sequence, start, stop):
    """Find the start codon.
    """
    return (start_regex, sequence, start, stop)


def find_stop(stop_regex, sequence, start):
    """Find the stop codon.
    """
    return (stop_regex, sequence, start)


def has_shine_dalgarno(shine_regex, sequence, start,
                       max_shine_dalgarno_distance):
    """Find a shine dalgarno motif before the start codon.
    """
    return (shine_regex, sequence, start, max_shine_dalgarno_distance)


def predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes.
    """
    return (sequence, start_regex, stop_regex, shine_regex, min_gene_len,
            max_shine_dalgarno_distance, min_gap)


def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions.
    """
    try:
        with open(predicted_genes_file, "wt", encoding="utf-8") as pred_g:
            predict_genes_writer = csv.writer(pred_g, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit(f"Error cannot open {predicted_genes_file}")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc,
                probable_genes_comp):
    """Write gene sequence in fasta format.
    """
    try:
        with open(fasta_file, "wt", encoding="utf-8") as fasta:
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(f">gene_{i + 1}{os.linesep}{2}"
                            f"{fill(sequence[gene_pos[0]-1:gene_pos[1]])}")
                trans = i

            trans += 1

            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(f">gene_{trans + 1 + j}{os.linesep}{2}"
                            f"{fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])}")
    except IOError:
        sys.exit(f"Error cannot open {fasta_file}")


def reverse_complement(kmer):
    """Get the reverse complement.
    """
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([complement[base] for base in kmer[::-1]])


if __name__ == "__main__":
    m_start_regex = re.compile("AT[TG]|[ATCG]TG")
    m_stop_regex = re.compile("TA[GA]|TGA")
    m_shine_regex = re.compile("A?G?GAGG|GGAG|GG.{1}GG")

    args = get_arguments()
