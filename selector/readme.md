This folder contains the file necessary to sort and protonates cystine molecules.


Histo_dist_chi3.pl creates a yaml which contains the n-n dist and Chi3 of each 
t of the specified file

selector.pl generates a list of unique ts based on the critera from Hist_dist_chi3 
and specified increments and puts the list in best

bndsort.pl sorts the ts listed in best into cyclic and non-cyclic cystine and prints 
xyzs of those ts into the sorted folder

p2rotonation.pl protonates the sorted xyz files and prints the resulting xyz
into the protonated folder

c2ons_dihe_five calculates the energy of the protonated energy with the five 
dihedrals constrained


