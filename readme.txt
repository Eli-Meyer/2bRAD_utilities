##############################################
2bRAD analysis pipeline version 2.0
This version includes modifications developed in 2015-2016:
	> Improved descrimination between similar loci in de novo reference building, using a ML tree approach
	> Mapping and basecalling with SAM file format, making it compatible with other mappers besides SHRiMP
	> Improved sensitivity for SNPs near the ends of alignments, using global alignments instead of local
Currently this version should be considered BETA. 
I make no guarantees about its accuracy, although its behaved well in my tests so far.
Please test it on real biological datasets and help me improve it by reporting any bugs you find!
eli.meyer@oregonstate.edu
##############################################

This repository contains a series of scripts used for our standard 2bRAD analysis pipeline.
You'll also need the "sequence utilities" and "sequence processing" repositories, and several 
other software packages or modules installed and in your path:
BioPerl			http://www.bioperl.org/wiki/Main_Page
SHRiMP			http://compbio.cs.toronto.edu/shrimp/
cross_match.manyreads	http://www.phrap.org/phredphrapconsed.html

The usage information for each script (usage.txt) can be reproduced by running each script without arguments.

Any questions or bugs? Please contact Eli Meyer: eli.meyer@science.oregonstate.edu.
