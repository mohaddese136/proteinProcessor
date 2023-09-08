# proteinProcessor
this app reads a list of pdb files in dataset folder and generates following output:
1- a FASTA file for each protein 
2- a AA-SS.txt file that contains each proteins sequnce
3- a SUMMARY.txt file thant contain information below:
  -- Total number of processed PDBs
  -- Average of processed protein lengths
  -- Standard deviation of processed protein lengths
  -- Names of processed proteins
  -- probability of the presence of each of the twenty amino acids in each of the three states of the type II structure (Helix Sheet coil)
  -- probability of binary sequences of second type structure elements
