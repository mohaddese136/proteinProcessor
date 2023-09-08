# proteinProcessor
this app reads a list of pdb files in dataset folder and generates following output:<br>
1- a FASTA file for each protein <br>
2- a AA-SS.txt file that contains each proteins sequnce<br>
3- a SUMMARY.txt file thant contain information below:<br>
  -- Total number of processed PDBs<br>
  -- Average of processed protein lengths<br>
  -- Standard deviation of processed protein lengths<br>
  -- Names of processed proteins<br>
  -- probability of the presence of each of the twenty amino acids in each of the three states of the type II structure (Helix Sheet coil)<br>
  -- probability of binary sequences of second type structure elements<br>
