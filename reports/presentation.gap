# Tal Rastopchin and Gabby Masini
# March 8, 2020

# load Miyamoto and Hanaki's elementary functions for association schemes
Read("research_project/association_scheme.gap");

# load Miyamoto and Hanaki's classified association schemes of order 05-07
# as05 is a list of the schemes of order 5
Read("research_project/schemes_order05.gap");

# prints with line at end
Println := function(object)
	Print(object);
	Print("\n");
end;

M2 := as05[2];

Print("Is M2 an association scheme ? : ");
Println(IsAssociationScheme(M2));
Println("What is the automorphism group of M2 ? :");
Println(AutomorphismGroupOfScheme(M2));
Print("What is the structure constant a_qrs ? ");
Println(IntersectionNumber(M2, 1,2,2));
