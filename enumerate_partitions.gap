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
S3 := [[0, 1, 2, 3, 4, 5],
			[1, 0, 5, 4, 3, 2],
			[2, 4, 0, 5, 1, 3],
			[3, 5, 4, 0, 2, 1],
			[5, 3, 1, 2, 0, 4],
			[4, 2, 3, 1, 5, 0]];

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

# # Determines whether or not two good partitions are isomorpic
# # 	autoGroup, the automorphism group of the scheme we are inspecting
# # 		knownPartition, the classified good partition
# # 		unknownPartition, the good partition we are unsure we have aready
# # 			discovred.
# IsIsomorphic := function (autoGroup, knownPartition, unknownPartition)
# 	local orbit, part, part2, auto, isIso, oneEq, isomorphicPart, partSentTo;
# 	if not(Length(knownPartition) = Length(unknownPartition)) then
# 		return false;
# 	fi;
# 	orbit := [];
# 	for auto in autoGroup do
# 		isomorphicPart := [];
# 		#boolean starts as true
# 		isIso := true;
# 		for part in knownPartition do
# 			partSentTo := OnPoints(part, auto);
# 			oneEq x:= false;
# 			for part2 in unknownPartition do
# 				if(part2 = partSentTo) then
# 					oneEq := true;
# 					break;
# 				fi;
# 			od;
# 			if oneEq then
# 				continue;
# 			else
# 				break;
# 			fi;
# 			#Add(isomp)
# 				#for loop that goes through unknownpartition
# 				#compares all values of un part  to partsentto
# 				#if one is equal break;
# 				#if none are equal boolean is false
# 			#if boolean = false
# 			#break;
# 		od;
# 		if isIso then
# 			return true;
# 		fi;
# 		#if boolean still true return true;
# 	od;
#    return false;
#  end;
#   # possibly run through all orbits genrated by knownPartition and compare
#   # w/ unknown partition to see if isomorphic

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

# returns true iff partition1 and partition2 are equivalent under
# the given automorphism
CheckEquivalent := function(automorphism, partition1, partition2)
	local part, permuted;
	for part in partition1 do
		permuted := PermutePart(automorphism, part);
		if not(permuted in partition2) then
			return false;
		fi;
  od;
	return true;
end;

# we are assuming that if the image of partition1 is a subset of
# partition2, then we must have that Im(partition1) = partition2
IsIsomorphic := function(autoGroup, partition1, partition2)
	local auto, isIso;
	if not(Length(partition1) = Length(partition2)) then
		return false;
	fi;
	for auto in autoGroup do
		isIso := CheckEquivalent(auto, partition1, partition2);
		if isIso then
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

ComputeGoodPartitions := function(R)
	local n, partitionIterator, jvectors, representatives, partition, partitions,representative, automorphisms, numGoodPartitions, previouslyDiscovered;
	Newline();
	n := OrderOfScheme(R);
	numGoodPartitions := 0;
	#partitionIterator := IteratorOfPartitions(n);
  partitions := PartitionsSet([1..n]);
	# the set of representatives of equivalence classes of good partitions
	representatives := [];
	automorphisms := AutomorphismGroupOfScheme(R);

	#for partition in partitionIterator do
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
	Println(numGoodPartitions);
	Println(Length(representatives));
	# for rep in representatives do
	# 	PrintMatrix(TransposedMat(rep));
	# od;
	# return representatives;
end;

#repS3 := ComputeGoodPartitions(S3);
# PrintMatrix(TransposedMat(PermutePartition((2,3)(5,6), repS3[2])));
#Println(IsIsomorphic(AutomorphismGroupOfScheme(S3), repS3[2], repS3[4]));

# ComputeGoodPartitions := function(R)
# 	local numGoodPartitions, representatives, automorphisms, partitions, partition, isGoodPartition, representative;
#   numGoodPartitions := 0;
# 	# enumerate potential partitions
# 	partitions := EnumeratePartitionBases(OrderOfScheme(R));
#
# 	# the set of representatives of equivalence classes of good partitions
# 	representatives := [];
# 	automorphisms := AutomorphismGroupOfScheme(R);
#
# 	# for each potential partition
# 	for partition in partitions do
#
# 		isGoodPartition := IsGoodPartition(R, partition);
#
# 		# if a good partition
# 		if isGoodPartition = true then
# 			# have we already discovered an isomorphic version of this?
# 			if (Length(representatives) = 0) then
# 				Add(representatives, partition);
# 			else
# 				for representative in representatives do
# 					if not(IsIsomorphic(automorphisms, representative, partition)) then
# 						Add(representatives, partition);
# 					fi;
# 	      od;
# 			fi;
# 			numGoodPartitions := numGoodPartitions + 1;
# 			#PrintMatrix(TransposedMat(partition));
# 		fi;
# 	od;
#
# 	Print("The number of good partitions is : ");
# 	Println(numGoodPartitions);
# 	Print("The number of good partition equivalence classes is : ");
# 	Println(Length(representatives));
# end;
