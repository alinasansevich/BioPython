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

# Create the top-level parser
parser = argparse.ArgumentParser(description='Search the Entrez database with a gene sequence accession number.')
parser.add_argument('-a', '--accession', type=str, metavar='', required=True, help='Gene sequence accession number')
subparsers = parser.add_subparsers(help='Choose file format (FASTA or GenBank full description) \n and action to perform (Print to Terminal or Save to file)')

# Create the parser for the "choose format" command
parser_f = subparsers.add_parser('file_format', help='Choose file format (FASTA or GenBank full description)')
parser_f.add_argument('-f', '--fasta', action='store', required=True, help='FASTA file format')
parser_f.add_argument('-g', '--genbank', action='store', required=True, help='GenBank full description file')

# create the parser for the "choose action" command
parser_d = subparsers.add_parser('display', help='Choose how to display data: Print to Terminal or Save to file')
parser_d.add_argument('-p', '--print', action='store', required=True, help='Print in Terminal')
parser_d.add_argument('-s', '--save', action='store', required=True, help='Save to file (.txt)')

group = parser.add_mutually_exclusive_group()
group.add_argument('-q', '--quiet', action='store_true', help='print quiet')
group.add_argument('-v', '--verbose', action='store_true', help='print verbose')
args = parser.parse_args()


def get_genbank_file(accession_num):
    """
    Fetch GenBank full description data for accession_num and 
    print to Terminal or Save to file as .txt.
    
    accession_num : str
        Gene sequence accession number.
    """
    handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="gb", retmode="text")
   
    if args.print:
        print(handle.read())     # print on terminal
    elif args.save:
        f = open("demofile.txt", "w")    # save to file
        f.write(handle.read())
        f.close()


def get_fasta(accession_num):
    """
    Fetch sequence in FASTA format for accession_num and 
    print to Terminal or Save to file as .fasta
    
    accession_num : str
        Gene sequence accession number.
    """
    
    if args.print:
        handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="fasta", retmode="text")
        print(handle.read())   # print on terminal
    elif args.save:
        filename = '{}.fasta'.format(accession_num)
        out_handle = open(filename, "w")
        in_handle = Entrez.efetch(db="nucleotide", id=accession_num, rettype="fasta")
        record = SeqIO.parse(in_handle, "fasta")
        SeqIO.write(record, out_handle, "fasta")
        in_handle.close()
        out_handle.close()


if __name__ == '__main__':
    if args.quiet:
        if args.fasta:
            get_fasta(args.accession)
            print('Done.')
        elif args.genbank:
            get_genbank_file(args.accession)
            print('Done.')
    elif args.verbose:
        if args.fasta:
            print("Fetching data...")
            get_fasta(args.accession)
            print('Done.')
        elif args.genbank:
            print("Fetching data...")
            get_genbank_file(args.accession)
            print('Done.')
        
            
        







# https://www.reddit.com/r/learnpython/comments/69okp4/argparse_one_argument_mutually_exclusive_to_two/






















# NM_119948.4

# https://stackoverflow.com/questions/17909294/python-argparse-mutual-exclusive-group
# create the top-level parser
# parser = argparse.ArgumentParser(prog='PROG')
# parser.add_argument('--foo', action='store_true', help='help for foo arg.')
# subparsers = parser.add_subparsers(help='help for subcommand')

# create the parser for the "command_1" command
# parser_a = subparsers.add_parser('command_1', help='command_1 help')
# parser_a.add_argument('a', type=str, help='help for bar, positional')

# create the parser for the "command_2" command
# parser_b = subparsers.add_parser('command_2', help='help for command_2')
# parser_b.add_argument('-b', type=str, help='help for b')
# parser_b.add_argument('-c', type=str, action='store', default='', help='test')
