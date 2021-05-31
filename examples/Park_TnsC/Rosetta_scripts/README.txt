We mostly follow this tutorial: https://faculty.washington.edu/dimaio/files/rosetta_density_tutorial_aug18.pdf
Looking at case #5: completing partial models guided by experimental density data.

#to create fragment file (following step 5A)

~/src/Rosetta/main/source/bin/grower_prep.default.macosgccrelease -in:file:fasta input.fasta -fragsizes 3 9 -fragamounts 100 100

#generates the fragment files that you need for RosettaES

other than that, you need the mapfile (usually map.mrc in our scripts) , the fasta file (input.fasta), and the input pdb (input.pdb). The pdb should match the input fasta file in terms of sequence, but it can be missing segments. Each segment is referred to by its index, starting from 1. RosettaES will go through the polypeptide chain to see what segments are missing and rebuild them in sequence.

things to check:

1. make sure the databases/executables are pointing to the right path in runES.sh and runRosettaES.sh

2. make sure the mapresolution in runES.sh is correct (currently set to 3.2 Angstrom)

3. in order to run RosettaES, execute the runRosettaES.sh script in the run-time directory, referring to the loop you want to rebuild by its index (starting from the N-terminus). 

