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
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "alinasansevich@gmail.com"     # Always tell NCBI who you are




##### First half ready:
# fetch the GenBank flat file: rettype="gb", retmode="text"
handle = Entrez.efetch(db="nucleotide", id="NM_119948.4", rettype="gb", retmode="text")
# print on terminal
print(handle.read())
# save to file
f = open("demofile.txt", "w")
f.write(handle.read())
f.close()




############ this works
# fetch FASTA: rettype="fasta", retmode="text"
handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="fasta", retmode="text")
# print on terminal
print(handle.read())



out_handle = open("NM_119948.4.fasta", "w")
in_handle = Entrez.efetch(db="nucleotide", id="NM_119948.4", rettype="fasta")
record = SeqIO.parse(in_handle, "fasta")
SeqIO.write(record, out_handle, "fasta")
in_handle.close()
out_handle.close()



