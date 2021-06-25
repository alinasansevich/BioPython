# search_entrez.py
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This is a CLI app that searches the Entrez database with a nucleotide sequence accession number.


## Table of contents
- [search_entrez.py](#search_entrezpy)
- [Table of contents](#table-of-contents)
- [Installation](#installation)
- [Usage](#usage)
- [TLDR](#tldr)
- [License](#license)

***

## Installation
[(Back to top)](#table-of-contents)

This project was developed using [Biopython 1.78](https://biopython.org/wiki/Download). You'll need to install this library before running the program.

***

## Usage
[(Back to top)](#table-of-contents)

On the terminal, navigate to the directory where you have `search_entrez.py` and then type or paste the following line:

    $ python3 search_entrez.py <your-accession-number-here>

The program will prompt you to enter a valid email address ([find out why here](https://www.ncbi.nlm.nih.gov/books/NBK25497/)):

    `To make use of NCBI's E-utilities, NCBI requires you to specify your email address with each request.`
    `Please enter your email address:`
    
It will then print to screen the full report for the provided accession number. The app can also print to screen the results in FASTA format, and/or save them in a file (.txt or .fasta). All available options are detailed below:

_Positional arguments:_

accession_n &rarr; Gene sequence accession number

_Optional arguments:_

-h, --help &rarr; Show help message and exit<br>
-f, --fasta &rarr; Save sequence as FASTA file format (.fasta)<br>
-p, --print &rarr; Print sequence in FASTA format in the Terminal<br>
-s, --save &rarr; Save GenBank full description to file (.txt)

## TLDR
[(Back to top)](#table-of-contents)

If you are just interested in copy/pasting your sequence in FASTA format from the terminal, execute the following line:

    $ python3 search_entrez.py <your-accession-number-here> -p

If you prefer to save it to file, execute this:

    $ python3 search_entrez.py <your-accession-number-here> -f

It will create a file in your directory named <your-accession-number-here>.fasta

Or you could do both:
    
    $ python3 search_entrez.py <your-accession-number-here> -f -p

***
    
## License
[(Back to top)](#table-of-contents)
    
This repository contains content developed by Alina Sansevich and is distributed under the MIT license.<br>
