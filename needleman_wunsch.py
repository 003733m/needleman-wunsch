import sys

def read_scoring_matrix(blosum62_file):
    with open(blosum62_file, 'r') as file:
        lines = file.readlines()

        amino_acids = lines[0].split() # First lines are amino acids



        matrix = {amino_acid: {} for amino_acid in amino_acids}
        for i in range(1, len(lines)):
            values = lines[i].split()
            amino_acid1 = values[0]
            for j in range(1, len(values)):
                amino_acid2 = amino_acids[j - 1]
                matrix[amino_acid1][amino_acid2] = int(values[j])



    return amino_acids, matrix


def needleman_wunsch(seq1, seq2, amino_acids, scoring_matrix, gap_open, gap_extension):

     # Match, Up, Diagonal, Right and Leftcharacters here for traceback
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score_matrix = [[0] * cols for _ in range(rows)]
    traceback_matrix = [[''] * cols for _ in range(rows)]

    for i in range(1, rows):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_extension #Fill the first rows with gap extensions here.
        traceback_matrix[i][0] = 'U'

    for j in range(1, cols):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_extension #Fill the first cols with gap extensions.
        traceback_matrix[0][j] = 'L'

    for i in range(1, rows):
        for j in range(1, cols):
            amino_acid1 = seq1[i-1]
            amino_acid2 = seq2[j-1]

            match = score_matrix[i-1][j-1] + scoring_matrix[amino_acid1][amino_acid2]
            delete = score_matrix[i-1][j] + gap_extension if traceback_matrix[i-1][j] != 'D' else score_matrix[i-1][j] + gap_open
            insert = score_matrix[i][j-1] + gap_extension if traceback_matrix[i][j-1] != 'R' else score_matrix[i][j-1] + gap_open

            max_score = max(match, delete, insert)
            if max_score == match:
                traceback_matrix[i][j] = 'M'  # Create rotations for tracebacking algo here.
            elif max_score == delete:
                traceback_matrix[i][j] = 'U'
            else:
                traceback_matrix[i][j] = 'L'

            score_matrix[i][j] = max_score

    aligned_seq1, aligned_seq2 = '', ''

    i, j = rows - 1, cols - 1


    while i > 0 or j > 0:   #Tracebacking
        if i > 0 and j > 0 and traceback_matrix[i][j] == 'M':
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1

            j -= 1

        elif i > 0 and traceback_matrix[i][j] == 'U':
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1

        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1


    return aligned_seq1, aligned_seq2, score_matrix[rows-1][cols-1]

def calculate_identity(aligned_seq1, aligned_seq2):
    match_count = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2))
    total_length = len(aligned_seq1)
    identity_percent = (match_count / total_length) * 100 # Percentile score
    return identity_percent

def print_alignment(seq1, seq2, aligned_seq1, aligned_seq2, score, identity_percent):
    print("Aligned Sequences:")
    print(aligned_seq1)
    print('|' * len(aligned_seq1))
    print(aligned_seq2)

    print("\nAlignment Score:", score)
    print("Percent Identity:", identity_percent, "%")

def main():

    input_file = sys.argv[1]

    with open(input_file, 'r') as file:
        sequences = file.read().splitlines()

    seq1, seq2 = sequences[0], sequences[1]

    # Reading the blosum62 matrix here
    amino_acids, scoring_matrix = read_scoring_matrix("blosum62.txt")

    # Gap scores
    gap_open = -10
    gap_extension = -5

    aligned_seq1, aligned_seq2, alignment_score = needleman_wunsch(seq1, seq2, amino_acids, scoring_matrix, gap_open, gap_extension)

    identity_percent = calculate_identity(aligned_seq1, aligned_seq2)

    print_alignment(seq1, seq2, aligned_seq1, aligned_seq2, alignment_score, identity_percent)


if __name__ == "__main__":
    main()

