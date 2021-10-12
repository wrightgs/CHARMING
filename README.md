# Overview

CHARMING.py (Codon HARMonizING) is an implementation of the codon harmonization algorithm outlined in the paper "CHARMING: Harmonizing synonymous codon usage to replicate a desired codon usage pattern." CHARMING replicates the codon usage patterns of a wild type (WT) sequence from its origin (native) species in a synonymous mutant expressed in a destination (host) organism. It can also be used to create synonymous mutants that match the codon usage pattern of the WT sequence in its origin species.


# Python version

CHARMING.py was written using Python version 3.6.10, and requires Scipy, NumPy, and Matplotlib.


# Arguments

CHARMING requires six arguments, with an optional seventh argument if the user wishes to harmonize to a destination species.

List of arguments in order:

1) seqFile - The path to a file containing the sequence you would like to harmonize in FASTA format. The first line of the file should be a header with a leading ">" and the sequence should follow on the subsequent line(s)

2) numMutants - The desired number of synonymous, harmonized mutants. This value must be an integer greater than 1.

3) windowSize - The size of the sliding window to use with the desired codon usage measure. This value must be a positive integer greater than 1. The recommended window size is 10.

4) measure - The desired codon calculator, must be either "MM" (for %MinMax) or "GM" (for geometric mean). See the associated paper for the distinction between %MinMax and the geometric mean calculation. 

5) outName - The name the user wishes to associate with the harmonization. This name will be added to the beginning of the file names of the files output by CHARMING. This argument should have no spaces and conform to the file name requirements of the system the script is run on.

6) codonTableOrigin - The path to the codon usage table file for the origin species. The file should have 64 lines, where each of the lines contains a codon followed by that codon's associated frequency, separated by whitespace. Therefore, the first three lines file should look like (note - codon order does not matter):
TTT 26.1
TCT 23.5
TAT 18.8
...

*7) codonTableDestination - (optional) The path to the codon usage table file for the destination species if harmonizing between two species. If this argument is omitted, CHARMING will create synonymous mutants of the input sequence that match the pattern of the sequence in the origin species. This file should be in the same format as codonTableOrigin.


# Output files


Once executed, CHARMING will output a minimum of three files to the working directory. Assuming that parameter 4 and 5 above were set to "MM" and "mySequence" respectively, the three files that would be output are "mySequence_MM_harmonized_sequences.txt", "mySequence_MM_harmonized_values.txt", and "mySequence_MM_harmonization_seq_1_plot.png" (one .png file will be created for each output harmonized mutant). 

File descriptions:

"mySequence_MM_harmonized_sequences.txt" contains the number of requested harmonized mutants, written in the FASTA format. Each output sequence has been assigned a number, and the header of each sequence contains that number, the GC content for that sequence, the net distance from the target values (Equation 1 in the corresponding paper), and the Pearson correlation with the target values. Sequences are sorted by net deviation to the target values in ascending order.

"mySequence_MM_harmonized_values.txt" contains the values of the chosen codon usage measure output for each harmonized sequence (in case the user wishes to analyze/visualize their harmonized mutants themselves). The first line contains the codon usage measure target values (i.e., the measure values for the WT sequence in the origin species). Each subsequent line contains the measure values for each output harmonized mutant. The numbers identifying each set of values correspond to the sequence numbers from mySequence_MM_harmonized_sequences.txt.

"mySequence_MM_harmonization_seq_1_plot.png" - a plot of the harmonized mutant's values vs the target values (labeled "WT"). One such plot will be created for each output harmonized mutant. The number after "seq" in the file name corresponds to the sequence numbers in the other two output files.



# Mini Tutorial

In the github repository with this script are files containing an E. coli gene (EcolTestSeq.txt), an E. coli codon usage table (EcolCUB.txt), and an S. cerevisiae table (ScerCUB.txt). Download these files and place them in the same folder as CHARMING.py. Then, navigate to this folder in your command line environment of choice. Ensure that you have Python 3 installed.

To run CHARMING to generate 3 output synonymous sequences, harmonizing from E. coli to S. cerevisiae using the %MinMax codon usage measure with window size 10, and saving the output files with names beginning with "tutorial", type:

python CHARMING.py EcolTestSeq.txt 3 10 MM tutorial EcolCUB.txt ScerCUB.txt

Once execution is complete, 5 total files should have been saved to your working directory containing info on 3 harmonized mutants!


If you would like to run CHARMING to generate 3 synonymous sequences for expression in E. coli instead of S. cerevisiae (keeping the other parameters the same), simply remove the last argument and type:

python CHARMING.py EcolTestSeq.txt 3 10 MM tutorial EcolCUB.txt



# CHARMING specifics

When harmonizing from one species to another while requesting multiple output sequences, CHARMING generates one sequence using Rodriguez initialization, one sequence harmonizing the WT sequence with no initialization, and all other output sequences by harmonizing randomly generated synonymous mutant sequences with no initialization. See the corresponding paper for more information.

When harmonizing mutants for expression in the origin organism (i.e., running CHARMING with no argument 7), all output sequences are the result of harmonizing randomly generated synonymous mutant sequences with no initialization. See the corresponding paper for more information.

Because some harmonized sequences are generated from randomly generated sequences (and because some randomly generated sequences are more difficult for CHARMING to harmonize well), CHARMING generates harmonized mutants for 10*N (where N is the requested number of output sequences) and outputs the N best harmonized mutants from this set. Here, harmonized mutant A is considered better than harmonized mutant B if the net deviation of A from the target values (Equation 1 in the corresponding paper) is less than the net deviation of B. 


