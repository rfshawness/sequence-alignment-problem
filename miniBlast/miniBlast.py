def generate_words(sequence, word_length):
    words = [sequence[i:i+word_length] for i in range(len(sequence) - word_length + 1)]
    return words

def generate_high_scoreing_neighbors(words, matrix):
    # for the given word, generate all possible neighboring words
    # if the score of the new word is greater than threshold, add it to the list of neighboring words  

    word_len = len(words[0])

    # Common values to use with BLOSUM62 matrix:
        # K-mer length = 2, threshold = 8
        # K-mer length = 3, threshold = 11
        # K-mer length = 4, threshold = 15
    if word_len == 2:
        threshold = 8
    elif word_len == 3:
        threshold = 11
    elif word_len == 4:
        threshold = 15
    else:
        threshold = 0

    # generate all 20^len(word) possible neighboring words
    # TODO: remove duplicates from the list and only consider unique words
    neighboring_words = set()
    for word in words:
        for i in range(word_len):
            for j in matrix.keys():
                new_word = word[:i] + j + word[i+1:]
                score = 0
                for k in range(word_len):
                    score += matrix[word[k]][new_word[k]]
                if score >= threshold:
                    neighboring_words.add(new_word)

    return list(neighboring_words)

def queryProcessing(matrix):
    # read the input file
    # the input file will be in a directory called "input"
    with open("./input/input.txt", "r") as file:
        lines = file.readlines()
        word_length = int(lines[0])
        sequence = lines[1].strip()

    # generate the words
    words = generate_words(sequence, word_length)

    return words



def find_matches_in_database(high_scoring_words, database_sequences, k):
    matches = {}
    for seq_id, db_seq in enumerate(database_sequences):
        for i in range(len(db_seq) - k + 1):
            k_mer = db_seq[i:i + k]
            if k_mer in high_scoring_words:
                if k_mer not in matches:
                    matches[k_mer] = []
                matches[k_mer].append((seq_id, i))
    return matches


def extend_match(match, query_sequence, database_sequence, scoring_matrix):
    left_extension = ""
    right_extension = ""
    current_score = scoring_matrix[match]
    max_score = current_score

    # Extend to the left
    i = 1
    while current_score >= max_score and i <= len(match):
        left_extension = query_sequence[-i] + left_extension
        current_score += scoring_matrix[query_sequence[-i]]
        if current_score > max_score:
            max_score = current_score
        i += 1

    # Reset current score
    current_score = scoring_matrix[match]

    # Extend to the right
    i = 1
    while current_score >= max_score and i <= len(match):
        right_extension = right_extension + database_sequence[i]
        current_score += scoring_matrix[database_sequence[i]]
        if current_score > max_score:
            max_score = current_score
        i += 1

    return left_extension + match + right_extension

def main():
   # BLOSUM 62 Scoring Matrix
   matrix = {
      "A": {"A": 4, "R": -1, "N": -2, "D": -2, "C": 0, "Q": -1, "E": -1, "G": 0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 0, "W": -3, "Y": -2, "V": 0},
      "R": {"A": -1, "R": 5, "N": 0, "D": -2, "C": -3, "Q": 1, "E": 0, "G": -2, "H": 0, "I": -3, "L": -2, "K": 2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
      "N": {"A": -2, "R": 0, "N": 6, "D": 1, "C": -3, "Q": 0, "E": 0, "G": 0, "H": 1, "I": -3, "L": -3, "K": 0, "M": -2, "F": -3, "P": -2, "S": 1, "T": 0, "W": -4, "Y": -2, "V": -3},
      "D": {"A": -2, "R": -2, "N": 1, "D": 6, "C": -3, "Q": 0, "E": 2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S": 0, "T": -1, "W": -4, "Y": -3, "V": -3},
      "C": {"A": 0, "R": -3, "N": -3, "D": -3, "C": 9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
      "Q": {"A": -1, "R": 1, "N": 0, "D": 0, "C": -3, "Q": 5, "E": 2, "G": -2, "H": 0, "I": -3, "L": -2, "K": 1, "M": 0, "F": -3, "P": -1, "S": 0, "T": -1, "W": -2, "Y": -1, "V": -2},
      "E": {"A": -1, "R": 0, "N": 0, "D": 2, "C": -4, "Q": 2, "E": 5, "G": -2, "H": 0, "I": -3, "L": -3, "K": 1, "M": -2, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
      "G": {"A": 0, "R": -2, "N": 0, "D": -1, "C": -3, "Q": -2, "E": -2, "G": 6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0, "T": -2, "W": -2, "Y": -3, "V": -3},
      "H": {"A": -2, "R": 0, "N": 1, "D": -1, "C": -3, "Q": 0, "E": 0, "G": -2, "H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2, "V": -3},
      "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
      "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I": 2, "L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
      "K": {"A": -1, "R": 2, "N": 0, "D": -1, "C": -3, "Q": 1, "E": 1, "G": -2, "H": -1, "I": -3, "L": -2, "K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
      "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q": 0, "E": -2, "G": -3, "H": -2, "I": 1, "L": 2, "K": -1, "M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
      "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I": 0, "L": 0, "K": -3, "M": 0, "F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1},
      "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
      "S": {"A": 1, "R": -1, "N": 1, "D": 0, "C": -1, "Q": 0, "E": 0, "G": 0, "H": -1, "I": -2, "L": -2, "K": 0, "M": -1, "F": -2, "P": -1, "S": 4, "T": 1, "W": -3, "Y": -2, "V": -2},
      "T": {"A": 0, "R": -1, "N": 0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 5, "W": -2, "Y": -2, "V": 0},
      "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F": 1, "P": -4, "S": -3, "T": -2, "W": 11, "Y": 2, "V": -3},
      "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H": 2, "I": -1, "L": -1, "K": -2, "M": -1, "F": 3, "P": -3, "S": -2, "T": -2, "W": 2, "Y": 7, "V": -1},
      "V": {"A": 0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I": 3, "L": 1, "K": -2, "M": 1, "F": -1, "P": -2, "S": -2, "T": 0, "W": -3, "Y": -1, "V": 4}
      }
   
   # read the input file and generate words
   words = queryProcessing(matrix)
   # generate neighboring words for each word
   word_neighbors = generate_high_scoreing_neighbors(words, matrix)
   print(word_neighbors)

   # check for matches in the database
   database_sequences = ["ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY"]
   k = 3
   matches = find_matches_in_database(word_neighbors, database_sequences, k)
   
   # extend the matches
   for match in matches:
        for seq_id, i in matches[match]:
            print(extend_match(match, words[seq_id], database_sequences[seq_id], matrix))

if __name__ == "__main__":
    main()