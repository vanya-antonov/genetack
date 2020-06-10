# Introduction
* GeneTack       --  predicts frameshifts in relatively short DNA fragments with all the genes located on the direct strand only
* GeneTack-GM    --  genome-wide frameshift prediction (the input genome sequence should be long enough to estimate the HMM parameters)
* mod2hmmdef.pl  --  converts .mod file generated by GeneMarkS to the .hmm_def format used by GeneTack

# GeneTack-GM prerequisites
* NOTE: you can skip this step and proceed to Installation if you only want to analyze short sequences using the GeneTack tool
* Install the GeneMarkS program package from http://topaz.gatech.edu/GeneMark/license_download.cgi
    - Make sure that the following tools are in your $PATH: gm, gmsn.pl, gmhmmp

# Installation Instructions
* Unpack the distribution archive and 'cd genetack-X.XX'
* Type 'make' and you are good to go!
* For convenience you can copy the distribution folder to a desired location and add it to your $PATH

# Examples
* GeneTack test run:
```
./genetack -m examples/model.hmm_def -f examples/fragment.fasta
```
* GeneTack-GM test run (takes ~ 3 min):
```
./genetack_gm.pl examples/genome_1Mb.fasta
```

* mod2hmmdef.pl test run:
    - Firt, run GeneMarkS to produce a model file that is called "GeneMark_hmm.mod" by default:
    `gmsn.pl --clean examples/genome_1Mb.fasta`
    - Next, convert the generated "GeneMark_hmm.mod" file into "GeneMark_hmm.hmm_def" using a template file:
      `./mod2hmmdef.pl GeneMark_hmm.mod hmm_def_files/genetack_gm.hmm_def GeneMark_hmm.hmm_def`
    - Finally, the produced "GeneMark_hmm.hmm_def" model can be used to run GeneTack:
      `./genetack -m GeneMark_hmm.hmm_def -f examples/fragment.fasta`

The distribution includes 3 tools:

* GeneTack -- predicts frameshifts in relatively short DNA fragments with all the genes located on the direct strand only

* GeneTack-GM -- genome-wide frameshift prediction (the input genome sequence should be long enough to estimate the HMM parameters)

* mod2hmmdef.pl -- converts .mod file generated by GeneMarkS to the .hmm_def format.

# Citation
Antonov I, Borodovsky M. Genetack: frameshift identification in protein-coding sequences by the Viterbi algorithm. J Bioinform Comput Biol. 2010 Jun;8(3):535-51.
https://pubmed.ncbi.nlm.nih.gov/20556861
