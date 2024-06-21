import os

# Take query sequence and return all possible words of length word_length
def generate_words(sequence, word_length):

    words = [(sequence[i:i+word_length], i) for i in range(len(sequence) - word_length + 1)]
    return words



# Take a list of words and generate all possible neighboring words with a score greater than threshold
def generate_high_scoreing_neighbors(words, scoring_matrix): 

    word_len = len(words[0][0]) # words are tuples of (word, index)

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
    for word, index in words:  # unpack word and index from each tuple
        for i in range(word_len):
            for j in scoring_matrix.keys():
                new_word = word[:i] + j + word[i+1:]
                score = 0
                for k in range(word_len):
                    score += scoring_matrix[word[k]][new_word[k]]
                if score >= threshold:
                    neighboring_words.add((new_word, index))  

    return list(neighboring_words)



# Read the input file and return the words and the query sequence
def queryProcessing():

    # the input file will be in a directory called "input"
    with open("./input/input.txt", "r") as file:
        lines = file.readlines()
        word_length = int(lines[0])
        query_sequence = lines[1].strip()

    # generate the words
    words = generate_words(query_sequence, word_length)

    return words, query_sequence



# Find all matches in the database and their positions
def find_matches_in_database(high_scoring_words, directory, file_names, k):

    # Initialize a dictionary to store the matches
    matches = {}

    # Process each file one at a time
    for file_name in file_names:
        # Read the sequence from the file
        with open(os.path.join(directory, file_name), "r") as file:
            database_sequence = file.readline().strip()  # Since each file has only one line
            print(f"Database sequence: {database_sequence}")

        # Search for the high scoring words in the database sequence
        for word_tuple in high_scoring_words:
            word, query_position = word_tuple
            print(f"Searching for word {word} in database sequence")
            for i in range(len(database_sequence) - k + 1):
                if database_sequence[i:i+k] == word:
                    # Store the match in the dictionary along with the file name and position
                    if word not in matches:
                        matches[word] = []
                    matches[word].append((file_name, i, query_position))
                    print(f"Match found at position {i} in database sequence")

    return matches



# Extend the found matches to the left and right
def extend_match(match, query_sequence, database_sequence, scoring_matrix, match_start_db, match_start_query):

    left_extension = ""
    right_extension = ""
    current_score = sum(scoring_matrix[query_sequence[match_start_query + i]][database_sequence[match_start_db + i]] for i in range(len(match))) # calculate initial score of match 
    total_score = 0
    match_length = len(match)

    # Exdend to the left only if next character alignment is greater or equal to 0
    i = 1
    while match_start_query - i >= 0 and match_start_db - i >= 0:
        if scoring_matrix[query_sequence[match_start_query - i]][database_sequence[match_start_db - i]] > 0:
            left_extension = query_sequence[match_start_query - i] + left_extension
            current_score += scoring_matrix[query_sequence[match_start_query - i]][database_sequence[match_start_db - i]]
        i += 1
    total_score = current_score
    
    # Reset current score
    current_score = 0

    # Extend to the right only if next character alignment is greater or equal to 0
    i = 0
    while match_start_query + match_length + i < len(query_sequence) and match_start_db + match_length + i < len(database_sequence):
        if scoring_matrix[query_sequence[match_start_query + match_length]][database_sequence[match_start_db + match_length]] > 0:
            right_extension = right_extension + query_sequence[match_start_query + match_length]
            current_score += scoring_matrix[query_sequence[match_start_query + match_length]][database_sequence[match_start_db + match_length]]
        i += 1
    total_score += current_score

    return left_extension + match + right_extension, total_score




def main():
   # BLOSUM 62 Scoring Matrix
   scoring_matrix = {
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
   
   # Read the input file and return the words and the query sequence
   words, query_sequence = queryProcessing()
   print(f"Words: {words}")
   k = len(words[0][0])  # k-mer length

   # Generate neighboring words for each word
   word_neighbors = generate_high_scoreing_neighbors(words, scoring_matrix)
   print(f"Query sequence: {query_sequence}")
   print(f"Word neighbors: {word_neighbors}")

   # Get the directory that the script is in, and then the directory where the sequence files are stored
   script_directory = os.path.dirname(os.path.abspath(__file__))
   directory = os.path.join(script_directory, "sequenceDatabase")
   print(f"Database directory: {directory}")
   file_names = os.listdir(directory) # get all the file names in the directory
   
   print(f"Database files: {file_names}")

   # Find all matches in the database
   matches = find_matches_in_database(word_neighbors, directory, file_names, k)
   print(matches)
   
   # extend the matches
   for match, positions in matches.items():
       for position in positions:
        file_name, match_start_db, match_start_query = position
        with open(os.path.join(directory, file_name), "r") as file:
            database_sequence = file.readline().strip()
        extended_match, score = extend_match(match, query_sequence, database_sequence, scoring_matrix, match_start_db, match_start_query)
        print(f"Extended match for {match} at db position {match_start_db} in file {file_name}: {extended_match} with score {score}")

if __name__ == "__main__":
    main()