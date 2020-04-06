###################################################
#
# Elementary functions for association schemes on GAP
#
# for GAP 4.x (http://www.gap-system.org/)
# Akihide Hanaki (hanaki@math.shinshu-u.ac.jp)
# 2013/10/08
#
###################################################

_IsGrapeLoaded := (LoadPackage("grape") <> fail);

###########################
###########################
# Elementary functions
###########################
###########################

AdjacencyMatrices := function(M)
    local A, i, j, k, n, d;

    n := Length(M);
    d := Maximum(List(M, Maximum));
    A := [];
    for i in [1..d+1] do
        A[i] := NullMat(n, n);
    od;
    for i in [1..n] do
        for j in [1..n] do
            A[M[i][j]+1][i][j] := 1;
        od;
    od;

    return A;
end;

###########################
IntersectionNumber := function(M, i, j, k)
    local a, b, n, x;

    n := Length(M);
    a := Position(M[1], k);
    b := Number([1..n], x -> (M[1][x]=i) and (M[x][a]=j));
    return b;
end;

###########################
IsAssociationScheme := function(M)
    local i, j, k, A, V, d, n, l, Mat2Vec;

    Mat2Vec := function(M)
        local i, V;

        V := M[1] * 1;
        for i in [2..Length(M)] do
            Append(V, M[i]);
        od;

        return V;
    end;

    n := Length(M);
    A := AdjacencyMatrices(M);
    d := Length(A);
    V := [];
    for i in [1..d] do
        Add(V, Mat2Vec(A[i]));
    od;

    if (IdentityMat(n) in A) = false then
        return false;
    fi;

    if ForAll([1..d], i -> (TransposedMat(A[i]) in A)) = false then
        return false;
    fi;

    for i in [1..d] do
        for j in [1..d] do
            V[d+1] := Mat2Vec(A[i]*A[j]);
            if RankMat(V) <> d then
                return false;
            fi;
        od;
    od;

    return true;
end;

###########################
IntersectionMatrices := function(M)
    local AM, i, j, k, IM, d, AM1, V, VB, P;

    AM := AdjacencyMatrices(M);
    d := Length(AM);
    AM1 := List([1..Length(AM)], x -> AM[x][1]);
    V := VectorSpace(Rationals, AM1);
    VB := Basis(V, AM1);
    IM := [];
    for i in [1..d] do
        IM[i] := [];
        for j in [1..d] do
            P := AM[i] * AM[j];
            IM[i][j] := Coefficients(VB, P[1]);
        od;
    od;

    return IM;
end;

###########################
###########################
# Construction
###########################
###########################

DirectProductScheme := function(M, N)
    local L, i, j, MM, NM, A;

    MM := AdjacencyMatrices(M);
    NM := AdjacencyMatrices(N);
    A := [];
    for i in MM do
        for j in NM do
            Add(A, KroneckerProduct(i, j));
        od;
    od;

    L := A[1] * 0;
    for i in [1..Length(A)] do
        L := L + A[i] * (i-1);
    od;

    return L;
end;

############################
WreathProductScheme := function(M, N)
    local L, i, j, MM, NM, A, J;

    MM := AdjacencyMatrices(M);
    NM := AdjacencyMatrices(N);
    A := [];
    for i in MM do
        Add(A, KroneckerProduct(NM[1], i));
    od;

    J := NullMat(Length(M), Length(M));
    for i in [1..Length(J)] do
        for j in [1..Length(J)] do
            J[i][j] := 1;
        od;
    od;

    for i in [2..Length(NM)] do
        Add(A, KroneckerProduct(NM[i], J));
    od;

    L := A[1] * 0;
    for i in [1..Length(A)] do
        L := L + A[i] * (i-1);
    od;

    return L;
end;

##############################################################
# G : group
# H : subgroup of G
# make the transitive permutation group scheme
# ---
# if G is given as a transitive permutation group,
# use "TransitivePermutationGroupScheme(G, Stabilizer(G,1))"
# ---
# to make the regular group scheme,
# use "TransitivePermutationGroupScheme(G, Subgroup(G,[()]))"
##############################################################
TransitivePermutationGroupScheme := function(G, H)
    local M, D, x, y, n, C, i, d;

    C := RightCosets(G, H);
    C := List(C, x -> Representative(x));
    n := Length(C);
    D := DoubleCosets(G, H, H);
    D := List(D, Elements);
    d := Length(D);
    M := NullMat(n, n);
    for x in [1..n] do
        for y in [1..n] do
            for i in [1..d] do
                if C[y]*C[x]^(-1) in D[i] then
                    M[x][y] := i - 1;
                fi;
            od;
        od;
    od;

    return M;
end;

SchurianScheme := TransitivePermutationGroupScheme;

##############################################################
# group association scheme
##############################################################
GroupAssociationScheme := function(G)
    local C, x, y, i, n, c, GE, M;

    C := ConjugacyClasses(G);
    C := List(C, Elements);
    c := Length(C);
    GE := Elements(G);
    n := Length(GE);
    M := NullMat(n, n);
    for x in [1..n] do
        for y in [1..n] do
            for i in [1..c] do
                if GE[x]^(-1)*GE[y] in C[i] then
                    M[x][y] := i - 1;
                fi;
            od;
        od;
    od;

    return M;
end;

##############################################################
NormalSubgroupScheme := function(G, N)
    local C, x, y, i, n, c, GE, M, CG;

    CG := ConjugacyClasses(G);
    CG := List(CG, Elements);
    C := [];
    for i in [1..Length(CG)] do
        if CG[i][1] in N then
            Add(C, ShallowCopy(CG[i]));
        fi;
    od;
    CG := [];
    c := Length(C);
    GE := Elements(N);
    n := Length(GE);
    M := NullMat(n, n);
    for x in [1..n] do
        for y in [1..n] do
            for i in [1..c] do
                if GE[x]^(-1)*GE[y] in C[i] then
                    M[x][y] := i - 1;
                fi;
            od;
        od;
    od;

    return M;
end;

##############################################################
# cyclotomic scheme
##############################################################
CyclotomicScheme := function(n, d)
    local i, m, FEr, x, y, M, r, F, FE;

    m := (n - 1) / d;
    F := GF(n);
    FE := Elements(F);
    FEr := [];
    for i in [1..d] do
        Add(FEr, Set(List([0..m], x -> Z(n)^(d*x+i))));
    od;
    Add(FEr, [0*Z(n)]);
    FEr := Set(FEr);
    M := NullMat(n, n);
    for x in [1..n] do
        for y in [1..n] do
            for i in [1..d+1] do
                if FE[x] - FE[y] in FEr[i] then
                    M[x][y] := i - 1;
                fi;
            od;
        od;
    od;

    return M;
end;

##############################################################
# Hamming scheme
##############################################################
HammingScheme := function(n, q)
    local S, x, y, i, d, M;

    S := Tuples([1..q], n);
    d := Length(S);
    M := NullMat(d, d);

    for x in [1..d] do
        for y in [1..d] do
            M[x][y] := Number([1..n], i-> S[x][i] <> S[y][i]);
        od;
    od;

    return M;
end;

##############################################################
# Johnson scheme
##############################################################
JohnsonScheme := function(v, k)
    local S, x, y, i, d, M;

    S := Combinations([1..v], k);
    d := Length(S);
    M := NullMat(d, d);

    for x in [1..d] do
        for y in [1..d] do
            M[x][y] := k - Length(Intersection(S[x], S[y]));
        od;
    od;

    return M;
end;

###########################
FusionScheme := function(M, L)
    local i, j, k, N, s, l, L2;

    if L[1] = [0] then
        L2 := L{[2..Length(L)]};
        L := L2;
    fi;
    s := Length(M);
    l := Length(L);
    N := NullMat(s, s);
    for i in [1..s] do
        for j in [1..s] do
            if i <> j then
                for k in [1..l] do
                    if M[i][j] in L[k] then
                        N[i][j] := k;
                    fi;
                od;
            fi;
        od;
    od;

    return N;
end;

###########################
SymmetrizationOfScheme := function(M)
    local i, L, n, N;

    n := Length(M);
    L := [];
    for i in [1..n] do
        Add(L, Set([M[i][1], M[1][i]]));
    od;

    L := Set(L);
    N := FusionScheme(M, L);

    return N;
end;

###########################
CompleteGraphScheme := function(n)
    local i, j, M;

    M := NullMat(n, n);
    for i in [1..n] do
        for j in [1..n] do
            M[i][j] := 1;
        od;
        M[i][i] := 0;
    od;

    return M;
end;

###########################
CompleteMultipartiteGraphScheme := function(n, m)
    local i, j, k, M, d;

    d := n / m;
    M := NullMat(n, n);
    for i in [1..n] do
        for j in [1..n] do
            M[i][j] := 2;
        od;
    od;

    for i in [1..d] do
        for j in [(m * (i - 1) + 1)..(m * i)] do
            for k in [(m * (i - 1) + 1)..(m * i)] do
                M[j][k] := 1;
            od;
        od;
    od;

    for i in [1..n] do
        M[i][i] := 0;
    od;

    return M;
end;

###########################
###########################
# Properties
###########################
###########################

IsSymmetricScheme := function(M)
    return M = TransposedMat(M);
end;

###########################
IsCommutativeScheme := function(M)
    local adj, i, j, d;

    if IsSymmetricScheme(M) then
        return true;
    fi;

    adj := AdjacencyMatrices(M);
    d := Length(adj);
    for i in [2..d] do
        if ForAll([2..d], j -> (adj[i] * adj[j] = adj[j] * adj[i])) = false then
            return false;
        fi;
    od;

    return true;
end;

###########################
NrCharacters := function(M)
    local i, j, k, l, m, n, AM, v, d, AM1, V, VB, U, P, C;

    n := Length(M);
    AM := AdjacencyMatrices(M);
    d := Length(AM);
    if IsCommutativeScheme(M) then
        return d;
    fi;
    v := List(AM, x -> Sum(x[1]));
    AM1 := List([1..Length(AM)], x -> AM[x][1]);
    V := VectorSpace(Rationals, AM1);
    VB := Basis(V, AM1);
    P := NullMat(d, d);
    for i in [1..d] do
        for j in [1..d] do
            for k in [1..d] do
                U := TransposedMat(AM[k]) * AM[i] * AM[k];
                C := Coefficients(VB, U[1]);
                P[i][j] := P[i][j] + C[j]/v[k];
            od;
        od;
    od;

    return Rank(P);
end;

NumCharacters := NrCharacters;

###########################
###########################
# Closed subsets
###########################
###########################

Valency := function(R, a)
    local x;

    return Number(R[1], x-> (x=a));
end;

###########################
Valencies := function(R)
    local L, i, x, d, n;

    d := Maximum(R[1]) + 1;
    n := Length(R);
    L := [];
    for i in [1..d] do
        L[i] := 0;
    od;
    for i in [1..n] do
        L[R[1][i] + 1] := L[R[1][i] + 1] + 1;
    od;

    return L;
end;

###########################
OrderOfScheme := function(R)
    return Length(R);
end;

###########################
SumOfValencies := function(R, L)
    local ans, i, V;

    V := Valencies(R);
    ans := 0;
    for i in L do
        ans := ans + V[i + 1];
    od;

    return ans;
end;

###########################
SumOfAdjacencyMatrices := function(R, L)
    local e, adj, i;

    adj := AdjacencyMatrices(R);
    e := adj[1] * 0;
    for i in L do
        e := e + adj[i + 1];
    od;

    return e;
end;

###########################
# R relation matrix, L closed subset
IsClosedSubset := function(R, L)
    local e, m;

    # sum of the adjacency matrices that have the relations in L
    # what will this look like ?
    e := SumOfAdjacencyMatrices(R, L);

    #  sum of entries of first row of e (row valency)
    m := Sum(e[1]);

    # we know the diagonal is all 1s bc L must contain the identity relation
    return (e^2 = m * e);
end;

###########################
IsNormalClosedSubset := function(R, L)
    local e, i, m, adj;

    e := SumOfAdjacencyMatrices(R, L);
    m := Sum(e[1]);

    if (e^2 <> m * e) then
        return false;
    fi;

    adj := AdjacencyMatrices(R);
    return ForAll(adj, i -> (e * i = i * e));
end;

###############################################
# local function, but some functions need this
###############################################
_ToOne := function(M)
    local i, j, n, M2;

    M2 := M * 0;
    n := Length(M);
    for i in [1..n] do
        for j in [1..n] do
            if M[i][j] <> 0 then
                M2[i][j] := 1;
            else
                M2[i][j] := 0;
            fi;
        od;
    od;

    return M2;
end;

###############################################
# local function, but some functions need this
###############################################
_Mat2List := function(R, M)
    local e, ans, n, i;

    ans := [];
    n := Length(M);
    e := _ToOne(M);
    for i in [1..n] do
        if e[1][i] <> 0 then
            AddSet(ans, R[1][i]);
        fi;
    od;

    return ans;
end;

###########################
GeneratedClosedSubset := function(R, L)
    local i, M, e, adj, m;

    if not (0 in L) then
        AddSet(L, 0);
    fi;
    e := SumOfAdjacencyMatrices(R, L);
    adj := AdjacencyMatrices(R);
    m := Sum(e[1]);

    while (e^2 <> m * e) do
        e := _ToOne(e^2);
        m := Sum(e[1]);
    od;

    return _Mat2List(R, e);
end;

###########################
ThinRadical := function(R)
    local L, V, i;

    V := Valencies(R);
    L := Filtered([1..Length(V)], i -> (V[i] = 1));
    for i in [1..Length(L)] do
        L[i] := L[i] - 1;
    od;

    return L;
end;

###########################
ThinResidue := function(R)
    local i, M, e, adj, m;

    e := R * 0;
    adj := AdjacencyMatrices(R);
    for i in adj do
        e := e + i * TransposedMat(i);
    od;
    e := _ToOne(e^2);
    m := Sum(e[1]);

    while (e^2 <> m * e) do
        e := _ToOne(e^2);
        m := Sum(e[1]);
    od;

    return _Mat2List(R, e);
end;

###########################
ComplexProduct := function(R, L1, L2)
    local e1, e2;

    e1 := SumOfAdjacencyMatrices(R, L1);
    e2 := SumOfAdjacencyMatrices(R, L2);
    return _Mat2List(R, e1 * e2);
end;

###########################
Neighbors := function(R, p, L)
    local i, ans, n;

    ans := [];
    n := Length(R);
    for i in [1..n] do
        if R[p][i] in L then
            AddSet(ans, i);
        fi;
    od;

    return ans;
end;

###########################
PartitionOfPointsByClosedSubset := function(R, L)
    local ans, n, i;

    ans := [];
    n := Length(R);
    for i in [1..n] do
        AddSet(ans, Neighbors(R, i, L));
    od;

    return ans;
end;

###########################
###########################
# Others
###########################
###########################

Mat2TeX := function(R)
    local i, j, m, n;

    m := Length(R);
    n := Length(R[1]);
    for i in [1..m] do
        for j in [1..n-1] do
            Print(R[i][j], " & ");
        od;
        Print(R[i][n], " \\\\\n");
    od;

    return;
end;

###########################
PermutePoints := function(R, P)
    local p;

    p := PermutationMat(P, Length(R));
    return p * R * TransposedMat(p);
end;

###########################
RenumberRelations := function(R, N)
    local i, j, m, n, R2;

    R2 := R * 1;
    m := Length(R);
    n := Length(R[1]);
    for i in [1..m] do
        for j in [1..n] do
            R2[i][j] := R[i][j]^N;
        od;
    od;

    return R2;
end;

###########################
RelationMatrixSortedByClosedSubset := function(R, L)
    local PP, Li, i;

    PP := PartitionOfPointsByClosedSubset(R, L);
    Li := [];
    for i in PP do
        Append(Li, i);
    od;
    PP := PermList(Li);

    return PermutePoints(R, PP);
end;

###########################
Subscheme := function(arg)
    local points, i, j, R2, m, n, S, R, L, p;

    R := arg[1];
    L := arg[2];
    if Length(arg) = 2 then
        p := 1;
    else
        p := arg[3];
    fi;

    points := Neighbors(R, p, L);
    n := Length(R);
    m := Length(points);
    S := List([1..m], i -> R[p][points[i]]);
    S := Set(S);
    R2 := NullMat(m, m);
    for i in [1..m] do
        for j in [1..m] do
            R2[i][j] := Position(S, R[points[i]][points[j]]) - 1;
        od;
    od;

    return R2;
end;

###########################
FactorScheme := function(R, L)
    local PP, i, j, m, n, q, R2, x, y, rel;

    n := Length(R);
    m := SumOfValencies(R, L);
    q := n / m;
    PP := PartitionOfPointsByClosedSubset(R, L);
    R2 := NullMat(q, q);
    for i in [1..q] do
        for j in [1..q] do
            R2[i][j] := [];
            for x in [1..m] do
                for y in [1..m] do
                    AddSet(R2[i][j], R[PP[i][x]][PP[j][y]]);
                od;
            od;
        od;
    od;

    rel := Set(R2[1]);

    for i in [1..q] do
        for j in [1..q] do
            R2[i][j] := Position(rel, R2[i][j]) - 1;
        od;
    od;

    return R2;
end;

###########################
QuotientScheme := function(R, L)
    return FactorScheme(R,L);
end;

###########################
IsThin := function(R)
    return Length(Set(Valencies(R))) = 1;
end;

###########################
IsQuasiThin := function(R)
    return Set(Valencies(R)) = [1, 2];
end;

###########################
IsPrimitiveScheme := function(M)
    local i, d;

    d := Maximum(M[1]);
    return ForAll([1..d], i -> (GeneratedClosedSubset(M, [i]) = [0..d]));
end;

###########################
IsSchurian := function(R)
    local G, adj, gp, gr, n, x, y, i;

    if _IsGrapeLoaded = false then
        return fail;
    fi;

    adj := AdjacencyMatrices(R);
    n := Length(R);
    G := SymmetricGroup(n);
    for i in [2..(Length(adj) - 1)] do
        gr := Graph(Group(()), [1..n], OnPoints, function(x,y) return adj[i][x][y]=1; end);
        gp := AutomorphismGroup(gr);
        G := Intersection(G, gp);
        if IsTransitive(G, [1..n]) = false then
            return false;
        fi;
    od;

    return (Length(Orbits(Stabilizer(G, 1), [1..n])) = Length(adj));
end;

###########################
AutomorphismGroupOfScheme:=function(R)
    local G, adj, gp, gr, n, x, y, i;

    if _IsGrapeLoaded = false then
        return fail;
    fi;

    adj := AdjacencyMatrices(R);
    n := Length(R);
    G := SymmetricGroup(n);
    for i in [2..(Length(adj) - 1)] do
        gr := Graph(Group(()), [1..n], OnPoints, function(x,y) return adj[i][x][y]=1; end);
        gp := AutomorphismGroup(gr);
        G := Intersection(G, gp);
    od;

    return G;
end;

###########################
IsPPolynomialScheme := function(R)
    local i, adj, d, gr, n, x, y;

    if _IsGrapeLoaded = false then
        return fail;
    fi;

    n := Length(R);
    adj := AdjacencyMatrices(R);
    d := Length(adj) - 1;
    for i in [2..d+1] do
        gr := Graph(Group(()), [1..n], OnPoints, function(x,y) return adj[i][x][y]=1; end);
        if IsDistanceRegular(gr) then
            if Diameter(gr) = d then
                return true;
            fi;
        fi;
    od;

    return false;
end;

###########################
PPolynomialOrdering := function(R, a)
    local ord, i, j, d;

    d := Maximum(R[1]);
    ord := [a];
    for j in [1..d-1] do
        for i in [1..d] do
            if not (i in ord) then
                if IntersectionNumber(R, a, ord[Length(ord)], i) <> 0 then
                    Add(ord, i);
                fi;
            fi;
        od;
    od;

    return ord;
end;

###########################
AllPPolynomialOrdering := function(R)
    local i, adj, d, gr, n, x, y, ans;

    if _IsGrapeLoaded = false then
        return fail;
    fi;

    ans := [];
    n := Length(R);
    adj := AdjacencyMatrices(R);
    d := Length(adj) - 1;
    for i in [2..d+1] do
        gr := Graph(Group(()), [1..n], OnPoints, function(x,y) return adj[i][x][y]=1; end);
        if IsDistanceRegular(gr) then
            if Diameter(gr) = d then
                Add(ans, PPolynomialOrdering(R, i - 1));
            fi;
        fi;
    od;

    return ans;
end;

###########################
IntersectionArray := function(R, a)
    local i, j, k, ord, d, ar, IM;

    d := Maximum(R[1]);
    IM := IntersectionMatrices(R);
    IM := IM[a + 1];
    ord := PPolynomialOrdering(R, a);
    ar := NullMat(3, d + 1);
    ar[1][1] := "*";
    ar[2][1] := 0;
    ar[3][1] := IM[ord[1] + 1][1         ];
    ar[1][2] := IM[1         ][ord[1] + 1];
    ar[2][2] := IM[ord[1] + 1][ord[1] + 1];
    ar[3][2] := IM[ord[2] + 1][ord[1] + 1];
    for i in [3..d] do
        ar[1][i] := IM[ord[i - 2] + 1][ord[i - 1] + 1];
        ar[2][i] := IM[ord[i - 1] + 1][ord[i - 1] + 1];
        ar[3][i] := IM[ord[i    ] + 1][ord[i - 1] + 1];
    od;
    ar[1][d + 1] := IM[ord[d - 1] + 1][ord[d] + 1];
    ar[2][d + 1] := IM[ord[d    ] + 1][ord[d] + 1];
    ar[3][d + 1] := "*";

    return ar;
end;

###########################
SymmetricRelations := function(R)
    local i, adj;

    adj := AdjacencyMatrices(R);
    return Filtered([0..(Length(adj) - 1)], i -> IsZero(adj[i+1] - TransposedMat(adj[i+1])));
end;

###########################
NonSymmetricRelations := function(R)
    local i, adj;

    adj := AdjacencyMatrices(R);
    return Filtered([1..(Length(adj) - 1)], i -> not (IsZero(adj[i+1] - TransposedMat(adj[i+1]))));
end;

###########################
Involutions := function(R)
    local i;

    return Filtered([1..Maximum(R[1])], i -> IsClosedSubset(R, [0, i]));
end;

##################################################
CharacterTableOfJohnsonScheme := function(arg)
    local i, j, T, a, v, k, flag;

    v := arg[1];
    k := arg[2];
    flag := -1;
    if Length(arg) = 3 then
        flag := arg[3];
    fi;

    if flag = -1 then
        T := NullMat(k+1, k+2);
    else
        T := NullMat(k+1, k+1);
    fi;

    for i in [0..k] do
        for j in [0..k] do
            for a in [0..i] do
                T[j+1][i+1] := T[j+1][i+1] + (-1)^(i-a) * Binomial(k-a, i-a)
                               * Binomial(k-j, a) * Binomial(v-k+a-j, a);
            od;
        od;
    od;

    if flag = -1 then
        for i in [0..k] do
            T[i+1][k+2] := (v-2*i+1) / (v-i+1) * Binomial(v, i);
        od;
    fi;

    return T;
end;

##################################################
CharacterTableOfHammingScheme := function(arg)
    local i, j, T, k, n, q, flag;

    n := arg[1];
    q := arg[2];
    flag := -1;
    if Length(arg) = 3 then
        flag := arg[3];
    fi;

    if flag = -1 then
        T := NullMat(n+1, n+2);
    else
        T := NullMat(n+1, n+1);
    fi;

    for k in [0..n] do
        for i in [0..n] do
            for j in [0..k] do
                T[i+1][k+1] := T[i+1][k+1] + (-1)^j * (q-1)^(k-j)
                               * Binomial(i, j) * Binomial(n-i, k-j);
            od;
        od;
    od;

    if flag = -1 then
        for i in [0..n] do
            T[i+1][n+2] := (q-1)^i * Binomial(n, i);
        od;
    fi;

    return T;
end;

##################################################
DRG2PPolynomialScheme := function(Gra)
    local M, d, n, i, j;

    d := Diameter(Gra);
    n := Length(Vertices(Gra));
    M := NullMat(n, n);
    for i in [1..n] do
        for j in [1..n] do
            M[i][j] := Distance(Gra, i, j);
        od;
    od;

    return M;
end;

##################################################
AMofDRG2PPolynomialScheme := function(A)
    local M, d, n, i, j, Gra;

    n := Length(A);
    Gra := Graph(Group(()), [1..n], OnPoints, function(i,j) return A[i][j]=1; end);
    d := Diameter(Gra);
    M := NullMat(n, n);
    for i in [1..n] do
        for j in [1..n] do
            M[i][j] := Distance(Gra, i, j);
        od;
    od;

    return M;
end;

#############################################################
ClassOfAssociationScheme := function(M)
    return Maximum(M[1]);
end;

#############################################################
TransposedRelations := function(M)
    local L, i, x;

    L := [];
    for i in [1..ClassOfAssociationScheme(M)] do
        L[i] := M[Position(M[1], i)][1];
    od;

    return L;
end;

#############################################################
CharacterTableOfAssociationScheme := function(arg)
    local nc, ct, d, i, j, k, am, A, n, Val, ct2, im, alg, Idem, M, F, flag;

    F := Rationals;
    flag := -1;
    M := arg[1];
    if Length(arg) > 1 then
        F := arg[2];
        if Length(arg) > 2 then
            flag := arg[3];
        fi;
    fi;

    nc := NrCharacters(M);
    im := IntersectionMatrices(M);
    alg := Algebra(F, im);
    Idem := CentralIdempotentsOfAlgebra(alg);

    if nc <> Length(Idem) then
        return false;
    fi;

    d := ClassOfAssociationScheme(M) + 1;
    ct := NullMat(nc, d + 1);
    am := AdjacencyMatrices(M);
    n := Length(M);
    Val := Valencies(M);
    for i in [1..nc] do
        ct[i][1] := Sqrt(Trace(Idem[i])); # degree
        A := NullMat(n, n);
        for j in [1..d] do
            A := A + am[j] * Idem[i][1][j];
        od;
        ct[i][d + 1] := Trace(A) / ct[i][1]; # multiplicity
        for j in [2..d] do
            ct[i][j] := Idem[i][1][j] * n * Val[j] / ct[i][d + 1];
        od;
    od;

    ct2 := NullMat(nc, d + 1);
    for j in [1..d] do
        ct2[1][j] := Val[j];
    od;
    ct2[1][d + 1] := 1;
    k := 2;
    for i in [1..nc] do
        if ct[i] <> ct2[1] then
            for j in [1..d+1] do
                ct2[k][j] := ct[i][j];
            od;
            k := k + 1;
        fi;
    od;

    for i in [2..nc] do
        for j in [nc, nc-1..i+1] do
            if ct2[j][d+1] / ct2[j][1] < ct2[j-1][d+1] / ct2[j-1][1] then
                k := ct2[j]; ct2[j] := ct2[j-1]; ct2[j-1] := k;
            fi;
        od;
    od;

    if flag <> -1 then
        for i in [1..nc] do
            Unbind(ct2[i][d + 1]);
        od;
    fi;

    return ct2;
end;

#############################################################
RegularCharacterOfAssociationScheme := function(R)
    local L, d, i, j, k, intnum;

    intnum := IntersectionMatrices(R);
    d := ClassOfAssociationScheme(R);
    L := [];
    for i in [0..d] do
        L[i + 1] := 0;
        for j in [1..(d+1)] do
            L[i + 1] := L[i + 1] + intnum[i + 1][j][j];
        od;
    od;

    return L;
end;

#############################################################
FrameNumber := function(T)
    local d, n, i, k, ans;

    ans := 1;
    d := Length(T[1]) - 2;
    k := Length(T);
    n := 0;
    for i in [1..(d+1)] do
        n := n + T[1][i];
    od;
    for i in [1..(d+1)] do
        ans := ans * T[1][i] * n;
    od;
    for i in [1..k] do
        ans := ans / (T[i][d + 2]^(T[i][1]^2));
    od;

    return ans;
end;

#############################################################
FrameQuotient := function(T)
    local d, n, i, k, ans;

    ans := 1;
    d := Length(T[1]) - 2;
    k := Length(T);
    n := 0;
    for i in [1..(d+1)] do
        n := n + T[1][i];
    od;
    for i in [1..(d+1)] do
        ans := ans * T[1][i] * n;
    od;
    for i in [1..k] do
        ans := ans / (T[i][d + 2]^(T[i][1]^2));
    od;

    return ans/(n^2);
end;

#############################################################
IsResiduallyThin := function(R)
    local tr, tr2;

    repeat
        if IsThin(R) then
            return true;
        fi;

        tr := ThinResidue(R);
        if tr = [0..Maximum(R[1])] then
            return false;
        fi;
        R := Subscheme(R, tr);
    until false;
end;

#############################################################
KreinParameters := function(T)
    local d, i, j, k, F, V, B, T2, Q, vec;

    d := Length(T);
    if Length(T[1]) <> d + 1 then
        return false;
    fi;

    T2 := NullMat(d, d);
    for i in [1..d] do
        for j in [1..d] do
            T2[i][j] := T[i][j];
        od;
    od;

    Q := [];
    for i in [1..d] do
        Q[i] := NullMat(d, d);
    od;

    F := DefaultFieldOfMatrix(T2);
    V := VectorSpace(F, T2);
    B := Basis(V, T2);

    for i in [1..d] do
        for j in [1..d] do
            vec := [];
            for k in [1..d] do
                vec[k] := T2[i][k] * T2[j][k] / T2[1][k];
            od;
            Q[i][j] := Coefficients(B, vec);
            for k in [1..d] do
                Q[i][j][k] := Q[i][j][k] * T[i][d+1] * T[j][d+1] / T[k][d+1];
            od;
        od;
    od;

    return Q;
end;

######################################
AllQPolynomialOrdering := function(T)
    local L, i, d, K, ord, TridiagonalOrdering;

    TridiagonalOrdering := function(M)
        local rest, pres, x, L, LT, d;

        d := Length(M);
        rest := [2..d];
        LT := [1];
        pres := 1;
        while rest <> [] do
            L := Filtered(rest, x-> M[pres][x] <> 0);
            if Length(L) <> 1 then
                return false;
            fi;
            pres := L[1];
            Add(LT, pres);
            rest := Difference(rest, [pres]);
        od;

        return LT;
    end;

    K := KreinParameters(T);
    if K = false then
        return [];
    fi;

    d := Length(T);
    L := [];
    for i in [2..d] do
        ord := TridiagonalOrdering(K[i]);
        if ord <> false then
            Add(L, ord);
        fi;
    od;

    return L;
end;

######################################
IsQPolynomialScheme := function(T)
    local ans;

    ans := AllQPolynomialOrdering(T);
    return (ans <> false) and (ans <> []);
end;

######################################
SemidirectProductScheme := function(R, G, hom)
    local nR, nG, dR, eG, RG, a, b, x, y, theta;

    nR := Length(R);
    dR := Maximum(R[1]) + 1;
    nG := Size(G);
    eG := Elements(G);
    RG := NullMat(nR*nG, nR*nG);
    for a in [1..nR] do
        for b in [1..nR] do
            for x in [1..nG] do
                for y in [1..nG] do
                    theta := eG[x]^(-1)*eG[y];
                    RG[a+nR*(x-1)][b+nR*(y-1)] :=
                      R[a^Image(hom,theta)][b]+dR*(Position(eG, theta)-1);
                od;
            od;
        od;
    od;

    return RG;
end;

######################################
TerwilligerAlgebra := function(arg)
    local am, F, R, i, j, p, d, nei, m, n;

    R := arg[1];
    if Length(arg) <= 1 then
        F := Rationals;
    else
        F := arg[2];
    fi;
    if Length(arg) <= 2 then
        p := 1;
    else
        p := arg[3];
    fi;

    am := AdjacencyMatrices(R);
    d := Maximum(R[1]);
    n := Length(R);
    for i in [0.. d] do
        nei := Neighbors(R, p, [i]);
        m := [];
        for j in [1..n] do
            if j in nei then
                m[j] := 1;
            else
                m[j] := 0;
            fi;
        od;
        m := DiagonalMat(m);
        Add(am, m);
    od;

    am := List(am, i -> i * One(F));
    return Algebra(F, am);
end;

######################################
AdjacencyAlgebra := function(arg)
    local am, F, R, i;

    R := arg[1];
    if Length(arg) = 1 then
        F := Rationals;
    else
        F := arg[2];
    fi;

    am := List(AdjacencyMatrices(R), i -> i * One(F));
    return Algebra(F, am);
end;

######################################
IntersectionAlgebra := function(arg)
    local im, F, R, i;

    R := arg[1];
    if Length(arg) = 1 then
        F := Rationals;
    else
        F := arg[2];
    fi;

    im := List(IntersectionMatrices(R), i -> i * One(F));
    return Algebra(F, im);
end;

#############################################################
IsCharacterTableAS := function(M, T)
    local e, i, j, k, AM, r, s, d, v, x, m, n, NT;

    e := []; m := []; n := Length(M);
    AM := AdjacencyMatrices(M);
    d := Length(AM);
    v := List(AM, x -> Sum(x[1]));
    r := NrCharacters(M);

    if (r <> Length(T)) then
        return false;
    fi;

    for i in [1..r] do
        s := 0;
        for j in [1..d] do
            s := s + T[i][j] * ComplexConjugate(T[i][j]) / v[j];
        od;
        m[i] := T[i][1] * n / s;
        e[i] := NullMat(n, n);
        for j in [1..d] do
            e[i] := e[i] + m[i] / n / v[j] * ComplexConjugate(T[i][j]) * AM[j];
        od;
    od;

    for i in [1..r] do
        if ForAll([1..d], j -> (e[i] * AM[j] = AM[j] * e[i])) = false then
            return false;
        fi;
    od;

    if ForAll([1..r], i -> (e[i] * e[i] = e[i])) = false then
        return false;
    fi;

    for i in [1..r] do
        if ForAll([(i + 1)..r], j -> (IsZero(e[i] * e[j]))) = false then
            return false;
        fi;
    od;

    if IsOne(Sum(e)) = false then
        return false;
    fi;

    if Length(T[1]) = d then
        return true;
    else
        return  ForAll([1..r], i -> (m[i] = T[i][d + 1]));
    fi;
end;

######################################
RenameRelations := function(M)
    local R, L, i, j, n;

    i := M[1][1];
    L := Set(M[1]);
    RemoveSet(L, i);
    L := Concatenation([i], L);
    n := Length(M);
    R := NullMat(n, n);
    for i in [1..n] do
        for j in [1..n] do
            R[i][j] := Position(L, M[i][j]) - 1;
        od;
    od;

    return R;
end;

######################################
Relation := function(R, x, y)
    return R[x][y];
end;

######################################
HadamardProduct := function(M, N)
    local i, j, L, m, n;

    m := Length(M);
    n := Length(M[1]);

    L := NullMat(m, n);
    for i in [1..m] do
        for j in [1..n] do
            L[i][j] := M[i][j] * N[i][j];
        od;
    od;

    return L;
end;

######################################
FrobeniusSchurIndicatorAS := function(R, T, a)
    local IM, i, M, d, n;

    d := Maximum(R[1]);
    n := Length(R);
    IM := IntersectionMatrices(R);
    M := IM[1];
    for i in [2..d+1] do
        M := M + IM[i]^2 / T[1][i];
    od;
    M := M[1];
    M[d + 2] := 0;
    M := M * T[a];
    M := M * T[a][d+2]/n/T[a][1];

    return M;
end;

######################################
CenterOfAdjacencyAlgebra := function(arg)
    local am, F, R, i;

    R := arg[1];
    if Length(arg) = 1 then
        F := Rationals;
    else
        F := arg[2];
    fi;

    return Center(AdjacencyAlgebra(R, F));
end;

######################################
IsGroupLikeScheme := function(R)
    local B, i, j, C, d;

    C := CenterOfAdjacencyAlgebra(R, Rationals);
    B := Elements(Basis(C));
    d := Length(B);
    for i in [1..d] do
        for j in [i..d] do
            if ((HadamardProduct(B[i], B[j]) in C) = false) then
                return false;
            fi;
        od;
    od;

    return true;
end;

######################################
CommutatorOfScheme := function(R, i, j)
    local e, adj;

    adj := AdjacencyMatrices(R);
    e := TransposedMat(adj[i+1]) * TransposedMat(adj[j+1])
         * adj[i+1] * adj[j+1];
    e := _ToOne(e);

    return _Mat2List(R, e);
end;

###########################
CommutatorOfSubsets := function(R, S1, S2)
    local e, adj, i, j, L;

    adj := AdjacencyMatrices(R);
    L := [];
    for i in S1 do
        for j in S2 do
            e := TransposedMat(adj[i+1]) * TransposedMat(adj[j+1])
                 * adj[i+1] * adj[j+1];
            e := _ToOne(e);
            L := Union(L, _Mat2List(R, e));
        od;
    od;

    return GeneratedClosedSubset(R, L);
end;

###########################
AllMultiplicityFreeSubgroups := function(G)
    local C, x, C2, C1, C3, H, H2, prevlen;

    C := ConjugacyClassesMaximalSubgroups(G);
    C := Filtered(C, x -> IsCommutativeScheme(SchurianScheme(G, Representative(x))));
    C1 := ShallowCopy(C);
    C2 := ShallowCopy(C);
    C3 := [];
    prevlen := 0;
    while true do
        for H in C2 do
            H2 := Representative(H);
            C := ConjugacyClassesMaximalSubgroups(H2);
            C := Filtered(C, x ->
                         IsCommutativeScheme(SchurianScheme(G, Representative(x))));
            Append(C3, C);
            C := [];
        od;
        C3 := List(C3, x ->  ConjugacyClassSubgroups(G, Representative(x)));
        C3 := Set(C3);
        Append(C1, C3);
        C1 := Set(C1);
        if prevlen = Length(C1) then
            return List(C1, Representative);
        fi;
        prevlen := Length(C1);
        C2 := ShallowCopy(C3);
        SubtractSet(C2, C1);
    od;

    return C1;
end;

###########################
AllMultiplicityFreeSelfPairedSubgroups := function(G)
    local C, x, C2, C1, C3, H, H2, prevlen;

    C := ConjugacyClassesMaximalSubgroups(G);
    C := Filtered(C, x -> IsSymmetricScheme(SchurianScheme(G, Representative(x))));
    C1 := ShallowCopy(C);
    C2 := ShallowCopy(C);
    C3 := [];
    prevlen := 0;
    while true do
        for H in C2 do
            H2 := Representative(H);
            C := ConjugacyClassesMaximalSubgroups(H2);
            C := Filtered(C, x ->
                         IsSymmetricScheme(SchurianScheme(G, Representative(x))));
            Append(C3, C);
            C := [];
        od;
        C3 := List(C3, x ->  ConjugacyClassSubgroups(G, Representative(x)));
        C3 := Set(C3);
        Append(C1, C3);
        C1 := Set(C1);
        if prevlen = Length(C1) then
            return List(C1, Representative);
        fi;
        prevlen := Length(C1);
        C2 := ShallowCopy(C3);
        SubtractSet(C2, C1);
    od;

    return C1;
end;

###########################
P_Valuation := function(n, p)
    local i;

return Number(FactorsInt(n), i -> i=p);
end;

###########################
P_ValuationRat := function(n, p)
    return P_Valuation(NumeratorRat(n), p) - P_Valuation(DenominatorRat(n), p);
end;

###########################
PValuationFiltration := function(R, p)
    local AM, val, L, x, y, d, L2, n;

    n := Length(R);
    val := Valencies(R);
    d := Length(val);
    L := [];
    for x in [1..n] do
        L[x] := [];
    od;
    for x in [1..d] do
        Add(L[PValuation(val[x], p)+1], x);
    od;
    y := 1;
    for x in [1..d] do
        if Length(L[x]) > 0 then
            y := x;
        fi;
    od;
    L2 := [];
    for x in [1..y] do
        L2[x] := L[x];
    od;

    return L2 - 1;
end;

###########################
IsPrimePower := function(n, p)
    if n = 1 then
        return true;
    fi;

    return Set(FactorsInt(n)) = [p];
end;

###########################
IsPValencedScheme := function(R, p)
    local Val, i;

    Val := Set(Valencies(R));
    for i in Val do
        if IsPrimePower(i, p) = false then
            return false;
        fi;
    od;

    return true;
end;

###########################
IsPPrimeValencedScheme := function(R, p)
    local Val, i;

    Val := Set(Valencies(R));
    for i in Val do
        if (i mod p) = 0 then
            return false;
        fi;
    od;

    return true;
end;

###########################
IsPScheme := function(R, p)
    local n;

    if IsPValencedScheme(R, p) = false then
        return false;
    fi;

    if IsPrimePower(OrderOfScheme(R), p) = false then
        return false;
    fi;

    return true;
end;

###########################
IsPiValencedScheme := function(R, pi)
    local Val, i;

    Val := Set(Valencies(R));
    for i in Val do
        if Intersection(FactorsInt(i), pi) <> [] then
            return false;
        fi;
    od;

    return true;
end;

###########################
# 2009/09/11
###########################
CosetDecompositionOfScheme := function(R, L)
    local PP, i, j, m, n, q, R2, x, y, rel;

    n := Length(R);
    m := SumOfValencies(R, L);
    q := n / m;
    PP := PartitionOfPointsByClosedSubset(R, L);
    R2 := NullMat(q, q);
    for i in [1..q] do
        for j in [1..q] do
            R2[i][j] := [];
            for x in [1..m] do
                for y in [1..m] do
                    AddSet(R2[i][j], R[PP[i][x]][PP[j][y]]);
                od;
            od;
        od;
    od;

    rel := Set(R2[1]);

    return rel;
end;

###########################
CosetRepresentativesOfScheme := function(R, L)
    local i;

    return List(CosetDecompositionOfScheme(R, L), i -> i[1]);
end;

###########################
CoefficientsOfAdjacencyMatrices := function(arg)
    local adj, b, vs, M, R, F;

    R := arg[1];
    M := arg[2];
    if Length(arg) < 3 then
        F := Rationals;
    else
        F := arg[3];
    fi;

    adj := AdjacencyMatrices(R) * One(F);
    vs := VectorSpace(F, adj);
    b := Basis(vs, adj);

    return Coefficients(b, M);
end;
###########################
IsCoherentConfiguration := function(M)
    local am, alg;

    am := AdjacencyMatrices(M);
    alg := Algebra(Rationals, am);

    return (Length(am) = Dimension(alg));
end;
###########################
AssociationScheme2Configuration := function(M, P)
    local em, r, i, j, k, n, am, am2, rm, tm;

    n := Length(M);
    r := Length(P);
    em := [];
    for i in [1..r] do
        em[i] := NullMat(n, n);
        for j in [1..Length(P[i])] do
            em[i][P[i][j]][P[i][j]] := 1;
        od;
    od;

    am := AdjacencyMatrices(M);
    am2 := [];
    for i in [1..r] do
        for j in [1..Length(am)] do
            for k in [1..r] do
                tm := em[i]*am[j]*em[k];
                if (IsZero(tm) = false) then
                    Add(am2, tm);
                fi;
            od;
        od;
    od;

    am2 := Set(am2);
    SubtractSet(am2, NullMat(n, n));

    rm := NullMat(n, n);
    for i in [1..Length(am2)] do
        rm := rm + (Length(am2) - i) * am2[i];
    od;

    return rm;
end;
###########################
IsTriplyRegularAssociationScheme := function(R)
    local am, d, n, i, nei, m, CC;

    am := AdjacencyMatrices(R);
    d := Maximum(R[1]);
    n := Length(R);
    nei := [[1]];
    for i in [1.. d] do
        nei[i+1] := Neighbors(R, 1, [i]);
    od;
    CC := AssociationScheme2Configuration(R, nei);

    return IsCoherentConfiguration(CC);
end;

###########################
AlgebraicAutomorphismGroupOfScheme := function(M)
    local IM, d, G, flag, T, newaut, a, genlist, isxinaut, b, c, i;

    isxinaut := function(y)
        local i, j, k;

        for i in [1..d] do
            for j in [1..d] do
                for k in [1..d] do
                    if IM[i^y][j^y][k^y] <> IM[i][j][k] then
                        return false;
                    fi;
                od;
            od;
        od;
        return true;
    end;

    newaut := function()
        local x;

        for x in T do
            if isxinaut(x) = true then
                return x;
            fi;
        od;
        return ();
    end;

    d := ClassOfAssociationScheme(M) + 1;
    IM := IntersectionMatrices(M);
    G := Subgroup(SymmetricGroup(d), [()]);
    genlist := [];

    while true do
        T := RightTransversal(SymmetricGroup(d), G);
        T := Set(T);
        RemoveSet(T, ());

        a := newaut();

        if a = () then
            genlist := List(genlist, ListPerm);
            for i in [1..Length(genlist)] do
                b := genlist[i];
                b := b{[2..Length(b)]};
                for c in [1..Length(b)] do
                    b[c] := b[c] - 1;
                od;
                b := PermList(b);
                genlist[i] := b;
            od;
            return Subgroup(SymmetricGroup(d), genlist);
        fi;

        Append(genlist, [a]);
        G := Subgroup(SymmetricGroup(d), genlist);
    od;
end;
###########################
AlgebraicAutomorphismGroupOfAssociationScheme := function(M)
    return AlgebraicAutomorphismGroupOfScheme(M);
end;

###########################
AlgebraicFusionOfScheme := function(M, G)
    local MF, orb, d;

    d := ClassOfAssociationScheme(M);
    orb := Orbits(G, [1..d]);

    return FusionScheme(M, orb);
end;

###########################
# 2012/04/08
###########################
CanonicalDualBasisAS := function(R)
    local AM, DAM, d, i, TR, val;

    d := ClassOfAssociationScheme(R) + 1;
    TR := TransposedRelations(R);
    val := Valencies(R);
    AM := AdjacencyMatrices(R);
    DAM := [];
    DAM[1] := AM[1];
    for i in [2..d] do
        DAM[i] := AM[TR[i-1]+1] / val[i];
    od;

    return [AM, DAM];
end;

###########################
CasimirOperatorAS := function(a, R)
    local CDB, i, ans;

    CDB := CanonicalDualBasisAS(R);
    ans := 0 * R;
    for i in [1..Length(CDB[1])] do
        ans := ans + CDB[1][i] * a * CDB[2][i];
    od;

    return ans;
end;

###########################
CasimirElementAS := function(R)
    return CasimirOperatorAS(R^0, R);
end;

###########################
VolumeAS := function(R)
    return CasimirOperatorAS(R^0, R);
end;

###########################
# 2012/04/11
###########################
CentralPrimitiveIdempotentsByCharacterTable := function(R, T)
    local idem, i, AM, d, n, j, val;

    AM := AdjacencyMatrices(R);
    val := Valencies(R);
    d := Length(AM);
    n := Length(R);
    idem := [];
    for i in [1..Length(T)] do
        idem[i] := 0 * R;
        for j in [1..d] do
            idem[i] := idem[i] + AM[j] * ComplexConjugate(T[i][j]) / val[j];
        od;
        idem[i] := idem[i] * T[i][d + 1] / n;
    od;

    return idem;
end;

###########################
CentralPrimitiveIdempotentByCharacterTable := function(R, T, i)
    local idem, AM, d, n, j, val;

    AM := AdjacencyMatrices(R);
    val := Valencies(R);
    d := Length(AM);
    n := Length(R);
    idem := 0 * R;
    for j in [1..d] do
        idem := idem + AM[j] * ComplexConjugate(T[i][j]) / val[j];
    od;
    idem := idem * T[i][d + 1] / n;

    return idem;
end;

###########################
ShrinkhandeGraphScheme := function()
    return
      [ [ 0,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2 ],
        [ 1,0,1,1,2,2,2,1,1,1,2,2,2,2,2,2 ],
        [ 1,1,0,2,1,2,2,1,2,2,1,1,2,2,2,2 ],
        [ 1,1,2,0,2,1,2,2,1,2,2,2,1,1,2,2 ],
        [ 1,2,1,2,0,2,1,2,2,2,1,2,1,2,1,2 ],
        [ 1,2,2,1,2,0,1,2,2,2,2,1,2,1,2,1 ],
        [ 1,2,2,2,1,1,0,2,2,1,2,2,2,2,1,1 ],
        [ 2,1,1,2,2,2,2,0,2,1,2,1,2,1,1,2 ],
        [ 2,1,2,1,2,2,2,2,0,1,1,2,1,2,2,1 ],
        [ 2,1,2,2,2,2,1,1,1,0,2,2,2,2,1,1 ],
        [ 2,2,1,2,1,2,2,2,1,2,0,1,1,2,2,1 ],
        [ 2,2,1,2,2,1,2,1,2,2,1,0,2,1,2,1 ],
        [ 2,2,2,1,1,2,2,2,1,2,1,2,0,1,1,2 ],
        [ 2,2,2,1,2,1,2,1,2,2,2,1,1,0,1,2 ],
        [ 2,2,2,2,1,2,1,1,2,1,2,2,1,1,0,2 ],
        [ 2,2,2,2,2,1,1,2,1,1,1,1,2,2,2,0 ] ];
end;

###########################
DoobGraphScheme := function(m, n)
    local R, i, DP;

    DP := function(A, B)
        local m, n, i1, i2, j1, j2, C;

        m := Length(A);
        n := Length(B);
        C := NullMat(m*n, m*n);
        for i1 in [1..m] do
            for i2 in [1..m] do
                for j1 in [1..n] do
                    for j2 in [1..n] do
                        C[(i1-1)*n+j1][(i2-1)*n+j2] := A[i1][i2] + B[j1][j2];
                    od;
                od;
            od;
        od;

        return C;
    end;

    R := ShrinkhandeGraphScheme();
    for i in [1..(m-1)] do
        R := DP(R, ShrinkhandeGraphScheme());
    od;
    for i in [1..n] do
        R := DP(R, CompleteGraphScheme(4));
    od;

    return R;
end;

###########################
CommutatorOfRingElements := function(a, b)
    return a*b-b*a;
end;

###########################
CommutatorOfAlgebra := function(arg)
    local I1, I2, B1, B2, i, j, ans, field;

    field := arg[1];
    I1 := arg[2]; B1 := Elements(Basis(I1));
    if Length(arg) = 2 then
        I2 := I1; B2 := B1;
    else
        I2 := arg[3]; B2 := Elements(Basis(I2));
    fi;

    ans := [];
    for i in B1 do
        for j in B2 do
            AddSet(ans, CommutatorOfRingElements(i, j));
        od;
    od;

    return VectorSpace(field, ans);
end;

###########################
NrModularSimpleModules := function(R, p)
    local field, A, K, J;

    if p = 0 then
        field := Rationals;
    else
        field := GF(p);
    fi;

    A := AdjacencyAlgebra(R, field);
    K := CommutatorOfAlgebra(field, A, A);
    J := RadicalOfAlgebra(A);

    return Dimension(A)-Dimension(K+J);
end;

###########################
STS2CC := function(M)
    local CC, i, j, k, v, b, sum;

    v := Length(M);
    b := Length(M[1]);

    CC := NullMat(v+b, v+b);

    for i in [1..v] do
        for j in [1..v] do
            if (i <> j) then
                CC[i][j] := 1;
            fi;
        od;
    od;

    for i in [1..v] do
        for j in [1..b] do
            CC[i][v+j] := 3 - M[i][j];
            CC[v+j][i] := 5 - M[i][j];
        od;
    od;

    for i in [1..b] do
        for j in [1..b] do
            if (i = j) then
                CC[v+i][v+j] := 6;
            else
                CC[v+i][v+j] := 8;
                for k in [1..v] do
                    CC[v+i][v+j] := CC[v+i][v+j] - M[k][i] * M[k][j];
                od;
            fi;
        od;
    od;

    return CC;
end;
