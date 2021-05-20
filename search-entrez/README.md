# search_entrez.py

This is a CLI program that searches the Entrez database with a sequence (nucleotide for now, I'll add protein searches later) accession number.
***
## What your project does



***
## How to install it

[Installation instructions for Biopython](https://biopython.org/wiki/Download)

Explain this:

    print("To make use of NCBI's E-utilities, NCBI requires you to specify your email address with each request.")
    Entrez.email = input("Please enter your email address: ")


***
## Example usage

    $ python3 search_entrez.py NM_119948.4D

_Positional arguments:_

accession_n &rarr; Gene sequence accession number &rarr; _EDIT to add protein_

_Optional arguments:_

-h, --help &rarr; Show help message and exit<br>
-f, --fasta &rarr; Save sequence as FASTA file format (.fasta)<br>
-p, --print &rarr; Print sequence in FASTA format in the Terminal<br>
-s, --save &rarr; Save GenBank full description to file (.txt)

***
## License and author info


continue following &rarr; https://dbader.org/blog/write-a-great-readme-for-your-github-project
