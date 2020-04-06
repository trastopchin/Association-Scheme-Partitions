
#simple association schemes w/ five points taken from Miyamoto and Hanaki's list
M1:=[ [ 0,1,1,1,1 ],
  [ 1,0,1,1,1 ],
  [ 1,1,0,1,1 ],
  [ 1,1,1,0,1 ],
  [ 1,1,1,1,0 ] ];

M2:=[ [ 0,1,1,2,2 ],
  [ 1,0,2,1,2 ],
  [ 1,2,0,2,1 ],
  [ 2,1,2,0,1 ],
  [ 2,2,1,1,0 ] ];

#Print(IsAssociationScheme(M1), "\n");
#M3:=CompleteGraphScheme(5);
#Print(M3, "\n");
#Print(OrderOfScheme(M1), "\n");

#Print some information about M2
Print("M2: \n");
Print(M2, "\n");
Print("Is association scheme?\n");
Print(IsAssociationScheme(M2), "\n");
Print("Order of scheme:\n");
Print(OrderOfScheme(M2), "\n");
Print("Automorphism group:\n");
Print(AutomorphismGroupOfScheme(M2), "\n");
Print("Algebraic automorphism group:\n");
Print(AlgebraicAutomorphismGroupOfScheme(M2), "\n");
Print("Transposed relations:\n");
Print(TransposedRelations(M2), "\n");
Print("Symmetric relations:\n");
Print(SymmetricRelations(M2), "\n");
Print("Non-symmetric relations:\n");
Print(NonSymmetricRelations(M2), "\n");
Print("Valencies:\n");
Print(Valencies(M2), "\n");
Print("Relation (1, 1):\n");
Print(Relation(M2, 1, 1), "\n");
Print("Relation (1, 2):\n");
Print(Relation(M2, 1, 2), "\n");
Print("Relation (2, 3):\n");
Print(Relation(M2, 2, 3), "\n");
Print("Involutions:\n");
Print(Involutions(M2), "\n");
Print("IntersectionNumber(1,1,1): \n");
Print(IntersectionNumber(M2, 1,1,1), "\n");
Print("IntersectionNumber(0,1,1): \n");
Print(IntersectionNumber(M2, 0,1,1), "\n");
Print("IntersectionNumber(0,1,2): \n");
Print(IntersectionNumber(M2, 0,1,2), "\n");
Print("IntersectionNumber(1,2,2): \n");
Print(IntersectionNumber(M2, 1,2,2), "\n");
Print("IntersectionNumber(1,1,2): \n");
Print(IntersectionNumber(M2, 1,1,2), "\n");
Print("Intersection matrices:\n");
Print(IntersectionMatrices(M2), "\n");
Print("Adjacency algebra over rationals:\n");
Print(AdjacencyAlgebra(M2), "\n");
Print("Some properties:\n");
Print("Is primitive scheme?\n");
Print(IsPrimitiveScheme(M2), "\n");
Print("Is commutative scheme?\n");
Print(IsCommutativeScheme(M2), "\n");
Print("Is symmetric scheme?\n");
Print(IsSymmetricScheme(M2), "\n");
Print("Is group-like scheme?\n");
Print(IsGroupLikeScheme(M2), "\n");
