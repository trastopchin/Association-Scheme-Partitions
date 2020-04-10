# Tal Rastopchin and Gabby Masini
# April 7, 2020

# Read("enumerate_partitions.gap");

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes of order 5 and 16
# as05 is a list of the schemes of order 5
Read("./classified_schemes/schemes_order05.gap");
Read("./classified_schemes/schemes_order06.gap");
Read("./classified_schemes/schemes_order07.gap");
Read("./classified_schemes/schemes_order08.gap");
Read("./classified_schemes/schemes_order09.gap");
Read("./classified_schemes/schemes_order10.gap");

# we need thispackage
LoadPackage("datastructures");

# prints a new line
Newline := function()
	Print("\n");
end;

# prints an object with a new line
Println := function(object)
	Print(object);
  Print("\n");
end;

# prints out a matrix with better formatting
PrintMatrix := function(A)
	local row;
	for row in A do
		Println(row);
	od;
	Print("\n");
end;

# makes a zero vector of dimension n
MyZeroVector := function (n)
  local zeroVector, i;
  zeroVector := [];
  # is there a faster way of doing this?
  for i in [1..n] do
    Add(zeroVector, 0);
  od;
  return zeroVector;
end;

# enumerates our partition bases
EnumeratePartitionBases := function (n)
  local partitionBases, partitions, partition, partitionBasis, part, partitionVector, element;

  partitionBases := []; # the set of partition bases
  partitions := PartitionsSet([1..n]); # compute partitions

  # each partition gives us a partition basis
  for partition in partitions do
    partitionBasis := []; # the current set of j vectors

    for part in partition do
      partitionVector := MyZeroVector(n); # the current j vector

      for element in part do
        partitionVector[element] := 1;
      od;

      Add(partitionBasis, partitionVector);
    od;
    Add(partitionBases, partitionBasis);
  od;

  return partitionBases;
end;

# a scheme of order 5
# scheme that comes from Z_5 ?
M:=[[0, 1, 2, 2, 1],
	  [1, 0, 1, 2, 2],
	  [2, 1, 0, 1, 2],
	  [2, 2, 1, 0, 1],
	  [1, 2, 2, 1, 0]];
# a scheme of order 6 that comes from S_3
S_3 := [ [ 0, 1, 1, 2, 2, 1 ],
         [ 1, 0, 2, 1, 1, 2 ],
				 [ 1, 2, 0, 1, 1, 2 ],
				 [ 2, 1, 1, 0, 2, 1 ],
				 [ 2, 1, 1, 2, 0, 1 ],
				 [ 1, 2, 2, 1, 1, 0 ] ];

# Creates a basis for the module given a partition.
# 	partition, the set of partition vectors of the all
#			one vector of length n.
VectorSpaceBasisFromPartition := function(partition)
	local V, basis;
	V := VectorSpace(Rationals, partition);
	basis := Basis(V, partition);
	return basis;
end;

# Computes the coefficients of a vector sigma * j
# in terms of the basis.
# 	sigma, an nxn adjacency matrix of a scheme R
#		j, a length n vector of the partition basis
#		basis, the partition basis of our vector space (module?)
ComputeCoefficients := function(sigma, j, basis)
	local result;
	result := sigma * j;
	return Coefficients(basis, result);
end;

# Determines whether or not a specific partition is
# a good partition of the given scheme. Does this by
# performing an exhaustive search.
# 	R, the relational mattrx of our given scheme.
# 	partition, the partition we are inspecting.
IsGoodPartition := function(R, partition)
	local basis, j, sigma, coefficients;
	# create the vector space over the rationals with our partition basis
	basis := VectorSpaceBasisFromPartition(partition);
	for j in partition do
		for sigma in AdjacencyMatrices(R) do
			# we compute the coefficients of sigma * j in terms of basis
			coefficients := ComputeCoefficients(sigma, j, basis);

			if coefficients = fail then
				return false;
			fi;
		od;
	od;
	return true;
end;

#
ComputeGoodPartitions := function(R)
	local numGoodPartitions, partitions, partition, isGoodPartition;
  numGoodPartitions := 0;
	# enumerate potential partitions
	partitions := EnumeratePartitionBases(OrderOfScheme(R));

	# for each potential partition
	for partition in partitions do

		isGoodPartition := IsGoodPartition(R, partition);

		# if a good partition
		if isGoodPartition = true then
			numGoodPartitions := numGoodPartitions + 1;
			PrintMatrix(TransposedMat(partition));
		fi;
	od;

	Print("The number of good partitions is : ");
	Println(numGoodPartitions);
end;

ComputeGoodPartitions(S_3);
