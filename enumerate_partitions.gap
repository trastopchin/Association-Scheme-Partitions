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

ComputeGoodPartitions := function(R)
	local partitionBases, partitionBasis, V, basis, isGoodPartition, j, sigma, result, coefficients;
	# enumerate potential partition bases
	partitionBases := EnumeratePartitionBases(OrderOfScheme(R));

	# for each potential partition basis
	for partitionBasis in partitionBases do

		# create the vector space over the rationals with our partition basis
		V:= VectorSpace(Rationals, partitionBasis);
		basis := Basis(V, partitionBasis);

		# check if all numbers b_pi^l exist
		isGoodPartition := true;
		for j in partitionBasis do
			for sigma in AdjacencyMatrices(R) do
				result := sigma * j; # sigma_p * j_i

				# figure out what the are the b_pi^l
				coefficients := Coefficients(basis, result);

				if not (coefficients = fail) then
					# coefficients are the b_pi^l
				else
					isGoodPartition := false;
					break; # break inner loop
				fi;
			od;

			# break outer loop
			if not(isGoodPartition) then
				break;
			fi;
		od;

		# if a good partition
		if isGoodPartition = true then
			Println("The partition basis ");
			PrintMatrix(TransposedMat(partitionBasis));
			Println("Is a good partitioning of ");
			PrintMatrix(R);
		fi;
	od;
end;

ComputeGoodPartitions(M);
