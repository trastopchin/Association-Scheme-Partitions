# Tal Rastopchin and Gabby Masini
# March 8, 2020

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./reports/association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes of order 05-07
# as05 is a list of the schemes of order 5, and so on
Read("./reports/schemes_order05.gap");
Read("./reports/schemes_order06.gap");
Read("./reports/schemes_order07.gap");

# prints with line at end
Println := function(object)
	Print(object);
	Print("\n");
end;

# prints out a matrix with better formatting
PrintMatrix := function(R)
	local row;
	for row in R do
		Println(row);
	od;
	Print("\n");
end;

# the relation matrix of the scheme example we drew
R :=
[[0, 1, 1, 1, 3, 2],
 [1, 0, 3, 2, 1, 1],
 [1, 2, 0, 3, 1, 1],
 [1, 3, 2, 0, 1, 1],
 [2, 1, 1, 1, 0, 3],
 [3, 1, 1, 1, 2, 0]];

# double check this is an association scheme
Print("Is R an association scheme : ");
Println(IsAssociationScheme(R));

# print out the adjacency matrices of R
Println("Ajacency matrices of R : ");
adjMats := AdjacencyMatrices(R);
for M in adjMats do
	PrintMatrix(M);
od;

# we were trying to figure out how this function works
# R relation matrix, L closed subset
MyIsClosedSubset := function(R, L)
    local e, m;

    # sum of the adjacency matrices that have the relations in L
    # what will this look like ?
    e := SumOfAdjacencyMatrices(R, L);

    #  sum of entries of first row of e (row valency)
    m := Sum(e[1]);

    # we know the diagonal is all 1s bc L must contain the identity relation
    return (e^2 = m * e);
end;
