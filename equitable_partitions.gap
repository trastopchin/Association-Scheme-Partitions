# Tal Rastopchin and Gabby Masini
# May 13, 2020

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./hanaki_association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes
# as05 is a list of the schemes of order 5
Read("./classified_schemes/schemes_order05.gap");
Read("./classified_schemes/schemes_order06.gap");
Read("./classified_schemes/schemes_order07.gap");
Read("./classified_schemes/schemes_order08.gap");
Read("./classified_schemes/schemes_order09.gap");
Read("./classified_schemes/schemes_order10.gap");
Read("./classified_schemes/schemes_order11.gap");
Read("./classified_schemes/schemes_order12.gap");

# makes a zero vector of dimension n
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

# Creates a basis for the module given a set of jvectors.
# 	jvectors, the set of partition vectors of the all
#			one vector of length n.
VectorSpaceBasisFromPartition := function(jvectors)
	local V, basis;
	V := VectorSpace(Rationals, jvectors);
	basis := Basis(V, jvectors);
	return basis;
end;

# Computes the coefficients of a vector sigma * jvector
# in terms of the basis. If the vector sigma * jvector does
# not lie in the vector space defined by the basis V,
# returns fail.
# 	sigma, an n x n adjacency matrix of a scheme R
#		jvector, a length n vector of the partition basis
#		basis, the partition basis of our vector space (module?)
ComputeCoefficients := function(sigma, jvector, basis)
	return Coefficients(basis, sigma * jvector);
end;

# Determines whether or not a specific set of jvectors is
# an equitable partition of the given scheme. Does this by
# performing an exhaustive search.
# 	R, the relational matrix of our given scheme.
# 	jvectors, the partition we are inspecting.
IsEquitablePartition := function(R, jvectors)
	local basis, jvector, sigma, coefficients;
	# create the vector space over the rationals with our partition basis
	basis := VectorSpaceBasisFromPartition(jvectors);
	for jvector in jvectors do
		for sigma in AdjacencyMatrices(R) do
			# we compute the coefficients of sigma * jvector in terms of basis
			coefficients := ComputeCoefficients(sigma, jvector, basis);
			# if coefficients don't exist, then return false
			if coefficients = fail then
				return false;
			fi;
		od;
	od;
	return true;
end;

# Applies a permutation to a jvector
#		permutation, our input permutation
#		jvector, the j vector we are permuting
PermuteJVector := function(permutation, jvector)
	local result, i;
	# duplicate jvector so we can edit it and create the permuted jvector
	result := ShallowCopy(jvector);
	for i in [1..Length(jvector)] do
		# i gets mapped to permuation applied to i
		result[i^permutation] := jvector[i];
	od;
	return result;
end;

# Returns true if and only if jvectors1 and jvectors2 are
# equivalent under the given automorphism
CheckEquivalent := function(automorphism, jvectors1, jvectors2)
	local jvector, permuted;
	# if every jvector of jvectors1 gets sent to a jvector belonging to
	# jvectors2 by the automorphism, then jvectors1 must be equivalent
	# to jvectors2.
	for jvector in jvectors1 do
		permuted := PermuteJVector(automorphism, jvector);
		if not(permuted in jvectors2) then
			return false;
		fi;
  od;
	return true;
end;

# Returns true if and only if jvectors1 and jvectors2 are
# equivalent under the given automorphism group, autoGroup
IsIsomorphic := function(autoGroup, jvectors1, jvectors2)
	local auto;
	# make sure jvectors have same cardinality
	if not(Length(jvectors1) = Length(jvectors2)) then
		return false;
	fi;
	for auto in autoGroup do
		if CheckEquivalent(auto, jvectors1, jvectors2) then
			return true;
		fi;
	od;
	# if every automorphism does not send jvectors1 to jvectors2,
	# then they are not isomorphic
	return false;
end;

# Applies a permutation to every jvector of a partition
#		permutation, our input permutation
#		jvectors, the set of jvectors we are permuting
PermuteJVectors := function (permutation, jvectors)
  local result, i;
	result := [];
	for i in [1..Length(jvectors)] do
		Add(result, PermuteJVector(permutation, jvectors[i]));
	od;
	return result;
end;

# Converts a partitions of integers [1..n] to
# a set of jvectors of length n of 1s and 0s
# depending on the indices indicated by the partition
# returns a set of jvectors based on partition
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
# partition sequence. The ith component of the
# auxiliary array corresponds to how many zeroes preceed
# the ith component in the array.
CreateAuxiliary := function(partitionSequence)
	local i, aux, val;
	# local variables
	aux := [];
	val := 0;
	# iterate through the sequence, creating the auxiliary
	for i in partitionSequence do
		Add(aux, val);
		# when we encounter a zero, set the appropriate value
		if i = 0 then
			val := val + 1;
		fi;
	od;
	return aux;
end;

# Given a partition sequence prev, produces the
# next partition sequence in our ordering of
# the partition sequences. Returns fail if we
# finished enumerating the partition sequences.
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
		if index > 0 then
			new[index] := new[index] + 1;
		else
			return fail;
		fi;
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

# Turns a set of jvectors into a partition
# The resulting partition is NOT necesarrily sorted
# by the first element of each part.
JVectorsToPartition := function(jvectors)
	local jvector, n, index, partition, currentPart;
	# figure out the size of the universe
	n := Length(jvectors[1]);
	partition := [];
	# iterate over each jvector
	for jvector in jvectors do
		currentPart := [];
		# create the parts of the partition (each jvector becomes a part)
		for index in [1..n] do
			if (jvector[index] = 1) then
				Add(currentPart, index);
			fi;
		od;
		# add the created part to the partition
		Add(partition, currentPart);
	od;
	return partition;
end;

# Turns a partition back into a partition sequence
# For the correct sequence to be output,
# the partition must be sorted by the first element
# of each part.
PartitionToSequence := function(partition, n)
  local sequence, part, index, count;
	sequence := MyZeroVector(n);
	count := 0;
  for part in partition do
		count := count + 1;
		# we skip the first element bc it is a zero
		# already made a zero by MyZeroVector
		for index in [2..Length(part)] do
      #skip first element somehow
			sequence[part[index]] := count;
		od;
	od;
	return sequence;
end;

# Given a set of jvectors, give its corresponding partition sequence
JVectorsToSequence := function(jvectors)
	local equivalentPartition, n;
	n := Length(jvectors[1]);
	# we turn the set of JVectors into a regular partition
	equivalentPartition := JVectorsToPartition(jvectors);
	# we sort the partition by its first element
	Sort(equivalentPartition, function(v,w) return v[1] < w[1]; end );
	# we turn the sorted partition back into the original sequence
	return PartitionToSequence(equivalentPartition, n);
end;

# Computes equitable partitions of scheme R
# Returns an output list with the first element being
# the number of equitable partitions and the second element
# being the list of representatives of each equivalence class
# of equitable partitions.
EquitablePartitions := function(R)
	local n, currentSequence, jvectors, representatives, partition, representative, automorphisms, numEquitablePartitions, previouslyDiscovered;

	n := OrderOfScheme(R);
	numEquitablePartitions := 0;

	# the set of representatives of equivalence classes of equitable partitions
	representatives := [];
	automorphisms := AutomorphismGroupOfScheme(R);

	# initialize our partition sequence
  currentSequence := MyZeroVector(n);

	# iterate through the partition sequences
	while not(currentSequence = fail) do
		# transform the current partition sequence into a set of jvectors
		partition := SequenceToPartition(currentSequence);
		jvectors := PartitionToJVectors(partition, n);

		# if the partition is a equitable partition
		if IsEquitablePartition(R, jvectors) = true then
			numEquitablePartitions := numEquitablePartitions + 1;
			# if this is the first equitable partition, add it to the representatives
			if (Length(representatives) = 0) then
				Add(representatives, jvectors);
				# if not, determine if it was previously discovered
			else
				previouslyDiscovered := false;

				# check if equitable partition is isomorphic to previously discovered equitable partitions
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

		# iterate to the next partition in the sequence
		currentSequence := NextPartitionSequence(currentSequence);
	od;
	# return the number of equitable partitions and the list of representatives
	return [numEquitablePartitions, representatives];
end;

# Computes equitable partitions of scheme R
# Returns an output list with the first element being
# the number of equitable partitions and the second element
# being the list of representatives of each equivalence class
# of equitable partitions.
# This does the same things as EquitablePartitions,
# but it is faster and uses a bit more space. (The accompanying paper
# explains this algorithm more deeply).
EquitablePartitionsFast := function(R)
	local added, equivalentJVectors, partitionSequence, auto, representativeIndex, equitablePartitionTable, n, currentSequence, jvectors, representatives, partition, automorphisms, numEquitablePartitions;

	# determine n
	n := OrderOfScheme(R);
	numEquitablePartitions := 0;

	# the set of representatives of equivalence classes of equitable partitions
	representativeIndex := 0;
	representatives := [];
	automorphisms := AutomorphismGroupOfScheme(R);

	# initialize our partition sequence
  currentSequence := MyZeroVector(n);

	# create a equitable partition table with currentSequence as the sample object
	equitablePartitionTable := HTCreate(currentSequence);

	# iterate through the partition sequences
	while not(currentSequence = fail) do
		partition := SequenceToPartition(currentSequence);
		jvectors := PartitionToJVectors(partition, n);

		# if it's not in the table, is it equitable?
		if HTValue(equitablePartitionTable, currentSequence) = fail then

			# if the partition is a equitable partition
			if IsEquitablePartition(R, jvectors) = true then

				# add to our list of representatives
				Add(representatives, jvectors);
				representativeIndex := representativeIndex + 1;

				# generate the entire equivalence class and add it to the table
				for auto in automorphisms do
					# apply the automorphism to get a set of equivalent jvectors
					equivalentJVectors := PermuteJVectors(auto, jvectors);
					# turn the equivalent jvectors into its corresponding partition sequence
					partitionSequence := JVectorsToSequence(equivalentJVectors);
					# add this partition sequence to the table
					added := HTAdd(equitablePartitionTable, partitionSequence, representativeIndex);
					# if we added a new equitable partition, increase the count
					if not(added = fail) then
						numEquitablePartitions := numEquitablePartitions + 1;
					fi;
				od;
			fi;
		fi;
		# iterate to the next partition in the sequence
		currentSequence := NextPartitionSequence(currentSequence);
	od;
	# return the number of equitable partitions and the list of representatives
	return [numEquitablePartitions, representatives];
end;

# Given a scheme and an equitable partition of that scheme, return
# a list of the entire equivalence class of that equitable partition.
# jvectors is assumed to be an equitable partition of that scheme.
ComputeIsomorphicPartitions := function(scheme, jvectors)
	local auto, automorphisms, equivalentJVectors, exampleSequence, partitionTable, partitionSequence, equivalenceClass, added;

	# compute automorphism group of scheme
	automorphisms := AutomorphismGroupOfScheme(scheme);

	# example sequence to initialize our hash table
  exampleSequence := MyZeroVector(OrderOfScheme(scheme));

	# create a partition table
	partitionTable := HTCreate(exampleSequence);

	# create our empty equivalence class
  equivalenceClass := [];

	# recall the first automorphism is always the identity automorphism
	for auto in automorphisms do
		# apply the automorphism to get a set of equivalent jvectors
		equivalentJVectors := PermuteJVectors(auto, jvectors);
		# turn the equivalent jvectors into its corresponding partition sequence
		partitionSequence := JVectorsToSequence(equivalentJVectors);
		# add this partition sequence to the table
		added := HTAdd(partitionTable, partitionSequence, auto);
		# if we added a new equitable partition, increase the count
		if not(added = fail) then
			Add(equivalenceClass, equivalentJVectors);
		fi;
	od;
	return equivalenceClass;
end;

# prints out a matrix A with better formatting
PrintMatrixToFile := function(file, A)
	local i;
	AppendTo(file, "[ ");
	for i in [1..Length(A)] do
		if not (i = 1) then
			AppendTo(file, "  ");
		fi;
		AppendTo(file, A[i]);
		if not(i = Length(A)) then
			AppendTo(file, ",\n");
	  else
			AppendTo(file, " ]");
		fi;
	od;
end;

# Given the output from EquitablePartitions(), prints the
# computed number of equitable partitions, the computed number
# of equivalence classes of equitable partitions, and prints
# a representative of each equivalence class.
# Each j vector is printed as a row vector, so accessing them
# is easier when they are read back in to GAP.
# schemeName is the name of the variable containing the list of
# good partitions.
# file represents the file name as a string that is either created if
# it doesn't exist or appended to if it does
# Using "*stdout*" for file prints to the terminal
# Mainly used as a helper function, but you can use it for
# more specific file output
PrintEquitablePartitionsToFile := function(file, output, schemeName)
	local i;
	AppendTo(file, "#----------------------------------------\n");
	AppendTo(file, "#Scheme name: \n");
	AppendTo(file, schemeName);
	AppendTo(file, " := ");
	AppendTo(file, "\n\n");
	AppendTo(file, "# Number of Equitable Partitions: ");
	AppendTo(file, output[1]);
	AppendTo(file, "\n");
	AppendTo(file, "# Number of Equivalence Classes: ");
	AppendTo(file, Length(output[2]));
	AppendTo(file, "\n");
  AppendTo(file, "[ \n");
	for i in [1..Length(output[2])] do
		AppendTo(file, "\n# No. ");
		AppendTo(file, i);
		AppendTo(file, "\n");
		PrintMatrixToFile(file, output[2][i]);
		if not(i = Length(output[2])) then
			AppendTo(file, ",\n");
		fi;
	od;
	AppendTo(file, "\n];\n");
end;

# Logs the equitable partitions of a specific scheme R to file
# schemeName is the name of the variable containing the list of
# good partitions. If file does not yet exist, it is created.
LogPartitionsScheme := function(file, R, schemeName)
	PrintEquitablePartitionsToFile(file, EquitablePartitionsFast(R), schemeName);
end;

# Logs the equitable partitions of list of schemes to file
# schemessName is the name of the prefix for each variable containing
# the list of good partitions. If file does not yet exist, it is created.
LogPartitionsSchemeList := function(file, schemes, schemesName)
	local scheme, count, schemeName;
	count := 0;
	for scheme in schemes do
		count := count+1;
		schemeName := Concatenation(schemesName, "No", String(count));
		LogPartitionsScheme(file, scheme, schemeName);
	od;
end;

# Given a scheme R, print to the terminal the output of running
# the EquitablePartitions algorithm.
PrintPartitionsScheme := function(R, schemeName)
	LogPartitionsScheme("*stdout*", R, schemeName);
end;

# Given a list of schemes, print to the terminal the output of running
# the EquitablePartitions algorithm on each scheme
# schemesName is the name of the prefix for each variable containing
# the list of good partitions.
PrintPartitionsSchemeList := function(schemes, schemesName)
	LogPartitionsSchemeList("*stdout*", schemes, schemesName);
end;
