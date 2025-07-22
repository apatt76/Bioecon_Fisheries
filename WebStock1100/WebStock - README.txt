Each .csv file is a table with n lines and n+2 colomns, with n the number of species of the foodweb. 
The matrix n*n is the binary adjacency matrix representing the foodlinks between each pair of species.
The n+1 colomn is a logical vector such that, 1 if the species is a fish, 0 otherwise.
The n+2 colomn is a vector containing the initial biomass of each species 