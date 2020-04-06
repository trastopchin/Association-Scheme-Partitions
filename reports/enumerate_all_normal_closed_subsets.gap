# Tal Rastopchin and Gabby Masini
# March 8, 2020

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("research_project/association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes of order 05-07
# as05 is a list of the schemes of order 5, and so on
Read("research_project/schemes_order05.gap");
Read("research_project/schemes_order06.gap");
Read("research_project/schemes_order07.gap");

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

# print new line
Print("\n");

# we choose which order schemes we will work with
relationMatrices := as06;

# we will store our normal closed subset as a list of sets
normalClosedSubsets := [];
for i in [1..Length(relationMatrices)] do
	normalClosedSubsets[i] := [];
od;

# let us find every normal closed subset of the schemes given in the relationMatrices
for i in [1..Length(relationMatrices)] do
	# get the current relation matrix and class of the scheme
	R := relationMatrices[i];
	class := ClassOfAssociationScheme(R);

	# generate all possible subsets of the relations
	powerSet := Combinations([1..class]);

	# keep track of the closed subsets we already generated
	generatedClosedSubsets := [];

	for subset in powerSet do
		# generate a closed subset
		t := GeneratedClosedSubset(R, subset);

		# if t is a new closed subset, check if its normal
		if not(t in generatedClosedSubsets) then
			# add t to the generated set of closed subsets
			AddSet(generatedClosedSubsets, t);
			# if t is normal, print t and its associated scheme
			if IsNormalClosedSubset(R, t) then
				AddSet(normalClosedSubsets[i], t);
			fi;
		fi;
	od;
od;

# print out the computed normal closed subsets
for i in [1..Length(relationMatrices)] do
	Print("The ");
	Print(i);
	Println("th scheme has the following normal closed subsets:");
	for ncs in normalClosedSubsets[i] do
		Print(ncs);
	od;
	Print("\n");
od;
Print("\n");
