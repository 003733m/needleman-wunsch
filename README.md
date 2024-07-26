# Needleman-Wunsch Algorithm
Needleman-Wunsh Algorithm Example (Global Search)
This is a pairwise sequence alignment of amino acid (protein) sequences via dynamic programming.
It includes: a scoring matrix (should be able to take any scoring scheme in the format of square a
matrix),
- a gap opening penalty value (a negative integer), and
- a gap extension penalty value (a negative integer).
# Background Information
The Needleman–Wunsch algorithm is an algorithm used in bioinformatics to align protein or nucleotide sequences. It was one of the first applications of dynamic programming to compare biological sequences. The algorithm was developed by Saul B. Needleman and Christian D. Wunsch and published in 1970.[1] The algorithm essentially divides a large problem (e.g. the full sequence) into a series of smaller problems, and it uses the solutions to the smaller problems to find an optimal solution to the larger problem.[2] It is also sometimes referred to as the optimal matching algorithm and the global alignment technique. The Needleman–Wunsch algorithm is still widely used for optimal global alignment, particularly when the quality of the global alignment is of the utmost importance. The algorithm assigns a score to every possible alignment, and the purpose of the algorithm is to find all possible alignments having the highest score. ![image](https://github.com/user-attachments/assets/ee567a6f-154b-4464-9829-b81b07358d88)


# Details in the comments...
