# Tal Rastopchin and Gabby Masini
# April 7, 2020

# Read("enumerate_partitions.gap");

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("./association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes of order 05-07
# as05 is a list of the schemes of order 5
Read("./classified_schemes/schemes_order05.gap");

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

  # Println(partitions);
  # Newline();
  # Println(partitionBases);

  return partitionBases;
end;

# https://en.wikipedia.org/wiki/Free_module

# might be a better way to do this using an algebraic
# construction from within GAP, but this is an exhaustive
# approach

# choose the second association schemes of order 5
R := as05[2];

# partition bases of order 5
partitionBases5 := EnumeratePartitionBases(5);

for partitionBasis in partitionBases5 do
  # create a module from our basis over Z
  basis := MutableBasis(Integers, partitionBasis);
  # check if current partition basis is closed
  closed := true;
  for j in partitionBasis do
    for sigma in AdjacencyMatrices(R) do
      # sigma_p * j_i
      result := sigma * j;
      # if it is in the span, even though the zero vector
      # is in the basis, that should not affect the computation
      if not(IsContainedInSpan(basis, result)) then
        closed := false;
      fi;
    od;
  od;

  # if all of the b numbers exist
  if closed = true then
    Println("The partition basis ");
    PrintMatrix(TransposedMat(partitionBasis));
    Println("Is a good partitioning of ");
    PrintMatrix(R);
  fi;
od;
