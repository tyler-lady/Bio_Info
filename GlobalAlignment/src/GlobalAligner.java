public class GlobalAligner {

    //Compute Matrix "F"
    d ← Gap penalty score
for i = 0 to length(A)
    F(i,0) ← d * i
for j = 0 to length(B)
    F(0,j) ← d * j
for i = 1 to length(A)
    for j = 1 to length(B)
    {
        Match ← F(i−1, j−1) + S(Ai, Bj)
        Delete ← F(i−1, j) + d
        Insert ← F(i, j−1) + d
        F(i,j) ← max(Match, Insert, Delete)
    }

    //"Backtrack"
    AlignmentA ← ""
    AlignmentB ← ""
    i ← length(A)
    j ← length(B)
while (i > 0 or j > 0)
    {
        if (i > 0 and j > 0 and F(i, j) == F(i−1, j−1) + S(Ai, Bj))
        {
            AlignmentA ← Ai + AlignmentA
            AlignmentB ← Bj + AlignmentB
            i ← i − 1
            j ← j − 1
        }
    else if (i > 0 and F(i, j) == F(i−1, j) + d)
        {
            AlignmentA ← Ai + AlignmentA
            AlignmentB ← "−" + AlignmentB
            i ← i − 1
        }
    else
        {
            AlignmentA ← "−" + AlignmentA
            AlignmentB ← Bj + AlignmentB
            j ← j − 1
        }
    }
}
