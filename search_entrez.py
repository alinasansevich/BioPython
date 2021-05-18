#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:04:43 2021

@author: alina

This program will search the Entrez database for a gene, using the accession number provided by the user.
The user may choose the results format (FASTA or GenBank full description file).

Resources:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
    http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec139
    http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3Aefetch
"""

import os
import argparse
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "alinasansevich@gmail.com"     # Always tell NCBI who you are

# Create parser
parser = argparse.ArgumentParser(description='Search the Entrez database with a gene sequence accession number.')
# Add positional argument
parser.add_argument('accession_num', metavar='accession_n', help='Gene sequence accession number')
# Add optional arguments
parser.add_argument('-f', '--fasta', action='store_true', help='Save sequence as FASTA file format (.fasta)')
parser.add_argument('-p', '--print', action='store_true', help='Print sequence in FASTA format in the Terminal')
parser.add_argument('-s', '--save', action='store_true', help='Save GenBank full description to file (.txt)')
args = parser.parse_args()

def get_genbank_file(accession_num):
    """
    Fetch GenBank full description data for accession_num and 
    print to Terminal or Save to file as .txt.
    
    accession_num : str
        Gene sequence accession number.
    """
    handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="gb", retmode="text")
    return handle.read()


def print_genbank(genbank):
    # save to file
    filename = '{}.txt'.format(args.accession_num)
    f = open(filename, "w")
    f.write(genbank)
    f.close()


def get_fasta(accession_num):
    """
    Fetch sequence in FASTA format for accession_num and save to file as .fasta
    
    accession_num : str
        Gene sequence accession number.
    """
    handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="fasta", retmode="text")
    # save to file
    filename = '{}.fasta'.format(accession_num)
    out_handle = open(filename, "w")
    in_handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="fasta")
    record = SeqIO.parse(in_handle, "fasta")
    SeqIO.write(record, out_handle, "fasta")
    in_handle.close()
    out_handle.close()


def print_fasta(accession_num):
    """
    Fetch sequence in FASTA format for accession_num and print to Terminal
    
    accession_num : str
        Gene sequence accession number.
    """
    handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="fasta", retmode="text")
    # print on terminal
    return handle.read()



if __name__ == '__main__':
    print('\n\n')
    genbank = get_genbank_file(args.accession_num)
    print(genbank)
    if args.save:
        print_genbank(genbank)
    if args.fasta:
        get_fasta(args.accession_num)
    if args.print:
        fasta = print_fasta(args.accession_num)
        print(fasta)
        


# NM_119948.4