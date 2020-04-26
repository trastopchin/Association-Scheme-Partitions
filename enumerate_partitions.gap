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

	#for each component of the partition, create a jvector
	for part in partition do

		# create a j vector of all zeros
		jvector := MyZeroVector(n);

		# if element is in partition set it's corresponding
		# element in the j vector to 1
		for element in part do
			jvector[element] := 1;
		od;

		# add this jvector to our set of jvectors
		Add(jvectors, jvector);
	od;
	return jvectors;
end;

# Creates the auxiliary array for a given
# partition number array. The ith component of the
# array corresponds to how many zeroes preceed
# the ith compnent in the array.
CreateAuxiliary := function(partition)
	local i, aux, val;
	# local variables
	aux := [];
	val := 0;
	# iterate through the partition, creating the auxiliary
	for i in partition do
		Add(aux, val);
		# when we encounter a zero, set the appropriate value
		if i = 0 then
			val := val + 1;
		fi;
	od;
	return aux;
end;

# Given a partition sequence prev, produces the
# next partition sequence in our ourdering of
# the partition sequences
NextPartitionSequence := function(prev)
	local aux, length, new, index;

	aux := CreateAuxiliary(prev);
	length := Length(prev);

	# our new partition sequence
	new  := ShallowCopy(prev);

	# if we can add one, add it
	if aux[length] > prev[length] then
     new[length]:= new[length] + 1;
	# otherwise find last spot in prev where aux is greater
	# and add 1 and make all of the following elements 0
	else
		index := length;

		# moves backwards along array
		while index>0 do
			if aux[index] > prev[index] then
				break;
			else
				# make the current element 0
				new[index] := 0;
			fi;
			index := index - 1;
		od;
		# increment the desired element
		new[index] := new[index] + 1;
	fi;
	return new;
end;

# Turns a partition sequence into a partition
SequenceToPartition := function(seq)
	local partition, index, component;
	# each partition starts like this
	partition := [[1]];
	for index in [2..Length(seq)] do
		component := seq[index];
		# if component is 0, we add a new part
		if component = 0 then
			Add(partition, [index]);
		# if not, we add into the corresponding part
		else
		  Add(partition[component], index);
		fi;
	od;
	return partition;
end;

# Computes good partitions of scheme R
# Returns an output list with the first element being
# the number of good partitions and the second element
# being the list of representatives of each equivalence class
# of good partitions.
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
