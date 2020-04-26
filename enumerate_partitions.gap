# Tal Rastopchin and Gabby Masini
# April 21, 2020

# Read("enumerate_partitions.gap");

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes
# as05 is a list of the schemes of order 5
Read("./classified_schemes/schemes_order05.gap");
Read("./classified_schemes/schemes_order06.gap");
Read("./classified_schemes/schemes_order07.gap");
Read("./classified_schemes/schemes_order08.gap");
Read("./classified_schemes/schemes_order09.gap");
Read("./classified_schemes/schemes_order10.gap");
Read("./classified_schemes/schemes_order11.gap");

# prints a new line
Newline := function()
	Print("\n");
end;

# prints an object with a new line
Println := function(object)
	Print(object);
  Print("\n");
end;

# prints out a matrix A with better formatting
PrintMatrix := function(A)
	local row;
	for row in A do
		Println(row);
	od;
	Print("\n");
end;

# makes a zero vector of dimension n
# TODO: can we do this in constant time?
MyZeroVector := function (n)
  local zeroVector, i;
  zeroVector := [];
  for i in [1..n] do
    Add(zeroVector, 0);
  od;
  return zeroVector;
end;

# Takes a group G and converts it to an association scheme
# with 0 relation in diagonal.
# Prints an error if it is not an association scheme,
# otherwise returns the scheme
GroupToScheme := function(G)
	local elements, dictionary, n, M, i, j, zeroIndex, row, element, index;

	# the first element of elements is always the Identity
	elements := Elements(G);
	n := Length (elements);

	# we use this dictionary to index our elements
	dictionary :=  NewDictionary(elements[1], true);
	# construct our hash table of element, relation pairs
	# (we do this so our algorithm is O(n^2) and not O(n^3))
  for i in [0..(n-1)] do
		AddDictionary(dictionary, elements[i+1], i);
	od;

	# this is our relation matrix
	M := MyZeroVector(n);

	# construct our relation matrix
	for i in [1..n] do
		zeroIndex := 0;
		row := [];
		for j in [1..n] do
			# we multiply every element of G by every other element of G
			element := elements[i]*elements[j];
			# look up its index in the DictionaryByList
			index := LookupDictionary(dictionary, element);
			Add(row, index);
			# we determine where to add the vector so the zero is along
			# the relational matrix's diagonal
			if index = 0 then
				zeroIndex := j;
			fi;
		od;
		M[zeroIndex] := row;
	od;

	# assert that our result is a relational matrix
	if not(IsAssociationScheme(M)) then
		Print("ERROR\n");
	else
  	return M;
	fi;
end;

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
# 	sigma, an n x n adjacency matrix of a scheme R
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
# 	R, the relational matrix of our given scheme.
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


# Applies a permutation to a partition j vector
#		permutation, our input permutation
#		part, the j vector we are permuting
PermutePart := function(permutation, part)
	local result, i;
	result := ShallowCopy(part);
	for i in [1..Length(part)] do
		# i gets mapped to permuation applied to i
		result[i^permutation] := part[i];
	od;
	return result;
end;

# Returns true if and only if partition1 and partition2 are
# equivalent under the given automorphism
CheckEquivalent := function(automorphism, partition1, partition2)
	local part, permuted;
	# if every part of partition1 gets sent to a part belonging to
	# partition2 by the automorphism, then partion1 must be equivalent
	# to partion2.
	for part in partition1 do
		permuted := PermutePart(automorphism, part);
		if not(permuted in partition2) then
			return false;
		fi;
  od;
	return true;
end;

# Returns true if and only if partition1 and partition2 are
# equivalent under the given automorphism group, autoGroup
IsIsomorphic := function(autoGroup, partition1, partition2)
	local auto;
	# make sure partitions have same cardinality
	if not(Length(partition1) = Length(partition2)) then
		return false;
	fi;
	for auto in autoGroup do
		if CheckEquivalent(auto, partition1, partition2) then
			return true;
		fi;
	od;
	# if every automorphism does not send partition1 to partition2,
	# then they are not isomorphic
	return false;
end;

# Applies a permutation to every j vector of a partition
#		permutation, our input permutation
#		partition, the partition we are permuting
PermutePartition := function (permutation, partition)
  local result, i;
	result := [];
	for i in [1..Length(partition)] do
		Add(result, PermutePart(permutation, partition[i]));
	od;
	return result;
end;

# Converts a partitions of integers [1..n] to
# a set of j vectors of length n of 1s and 0s
# depending on the indices indicated by the partition
# returns a set of j vectors based on partition
PartitionToJVectors := function(partition, n)
	local jvectors, part, element, jvector;
	jvectors := [];

	for part in partition do
		jvector := MyZeroVector(n); # the current j vector

		for element in part do
			jvector[element] := 1;
		od;

		Add(jvectors, jvector);
	od;
	return jvectors;
end;

# Enumerates the ordered partitions
CreatePartitions := function(n)
	local v, current, temp, partitions, index, base;
	partitions := [];
	v := MyZeroVector(n);
	Add(partitions, ShallowCopy(v));
  base := n;

	while base > 0 do
		# add 1 in base n
		v[n] := v[n] + 1;
		index := n;

		# carry if we go over the base
		if v[index] > (base-1) then
			base := base-1;
		fi;
		while v[index] > base-1 do
			v[index] := 0;
			# carry
			index := index - 1;
			v[index] := v[index] + 1;
		od;

		Add(partitions, ShallowCopy(v));
	od;
	return partitions;
end;

# CreatePartitions := function(n)
# 	local v, current, temp, partitions, index, base, k;
# 	partitions := [];
# 	v := MyZeroVector(n);
# 	Add(partitions, ShallowCopy(v));
# 	base := n;
# 	while base > 0 do
# 		# add 1 in base n
# 		v[n] := v[n] + 1;
# 		index := n;
#
# 		# carry if we go over the base
# 		# if v[index] > (base-1) then
# 		# 	base := base-1;
# 		# fi;
# 		while v[index] > base-1 do
# 			Print(v);
# 			Print("\n");
# 			v[index] := 0;
# 			# carry
# 			index := index - 1;
# 			k := v[index] + 1;
# 			v[index] := k;
# 		od;
# 		base := base -1;
# 		Add(partitions, ShallowCopy(v));
# 	od;
# 	return partitions;
# end;
"
[ [ 0, 0, 1 ], [ 0, 0, 2 ], [ 0, 1, 0 ] ]

	000, 001, 002, 010, 011.

	The way I would do this is to start with all 0’s.
	(This corresponds to the partition where every cell is a single element).
	 Then, start increasing the last digit until it equals the number of 0’s
	 that precede it.  Then, make the last digit a 0 and the second-to-last
	 digit a 1.  Then, repeat with the last digit.  That is, do the same thing
	 we do when we are counting in decimal notation, but rather than flipping
	 back to 0 when a given digit is a 9, flip to 0 when a given digit is equal
	 to the number of 0’s that precede it.  This is exactly how the numbers are
	 ordered for the 4-digit sequences listed above.";

#ConvertToPartition := function(v)

# Computes good partitions of scheme R
# Returns an output list with the first element being
# the number of good partitions and the second element
# being the list of representatives of each equivalence class
# of good partitions
ComputeGoodPartitions := function(R)
	local n, jvectors, representatives, partition, partitions,representative, automorphisms, numGoodPartitions, previouslyDiscovered;
	Newline();

	n := OrderOfScheme(R);
	numGoodPartitions := 0;

  partitions := PartitionsSet([1..n]);

	# the set of representatives of equivalence classes of good partitions
	representatives := [];
	automorphisms := AutomorphismGroupOfScheme(R);

	# for partition in partitionIterator do
	for partition in partitions do
		jvectors := PartitionToJVectors(partition, n);

		# if a good partition
		if IsGoodPartition(R, jvectors) = true then
			numGoodPartitions := numGoodPartitions + 1;
			# if this is the first good partition, add it to the representatives
			if (Length(representatives) = 0) then
				Add(representatives, jvectors);
				# if not, determine if it was previously discovered
			else
				previouslyDiscovered := false;

				# check if good partition is isomorphic to previously discovered good partitions
				for representative in representatives do
					if IsIsomorphic(automorphisms, representative, jvectors) then
						previouslyDiscovered := true;
						break;
					fi;
				od;

				#if nothing else already in the list is isomorphic, then add jvectors
				if not(previouslyDiscovered) then
					Add(representatives, jvectors);
				fi;
			fi;
		fi;
	od;
	Print("Good Partitions: ");
	Println(numGoodPartitions);
	Print("Equivalence Classes: ");
	Println(Length(representatives));
	return [numGoodPartitions, representatives];
end;
