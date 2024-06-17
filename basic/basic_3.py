# # CSCI 570 Final Project 
# # Sequence Alignment Problem 


import sys
from resource import * 
import time
import psutil

def process_memory():
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_consumed = int(memory_info.rss/1024)
    return memory_consumed

def time_wrapper(s, t, delta, alphas):
    start_time = time.time()
    cost, s_alignment, t_alignment = basic_algorithm(s, t, delta, alphas)
    end_time = time.time()
    time_taken = (end_time - start_time)*1000
    return time_taken, cost, s_alignment, t_alignment




def generate_strings(input_file):
    with open(input_file) as f:
        lines = f.readlines()
        s = lines[0].strip()
        for i in range(1, len(lines)):
            # if line is an integer
            if lines[i].strip().isdigit():
                n = int(lines[i].strip())
                s = s[:n+1] + s + s[n+1:]
            # if line is an string
            elif lines[i].strip().isalpha():
                t_start = i
                break
            else:
                # this is an error
                print("Invalid input")
        t = lines[t_start].strip()
        for i in range(t_start + 1, len(lines)):
            if lines[i].strip().isdigit():
                n = int(lines[i].strip())
                t = t[:n+1] + t + t[n+1:]
            else:
                # this is an error
                print("Invalid input")
    return s, t


delta = 30
alphas = {
    'A': {'A': 0, 'C': 110, 'G': 48, 'T': 94},
    'C': {'A': 110, 'C': 0, 'G': 118, 'T': 48},
    'G': {'A': 48, 'C': 118, 'G': 0, 'T': 110},
    'T': {'A': 94, 'C': 48, 'G': 110, 'T': 0}
}


def basic_algorithm(s, t, delta, alphas):
    m = len(s)
    n = len(t)
    dp = [[0 for _ in range(n+1)] for _ in range(m+1)]

    # initalize the zeroth row and column with delta*i and delta*j
    for i in range(1, m+1):
        dp[i][0] = dp[i-1][0] + delta
    for j in range(1, n+1):
        dp[0][j] = dp[0][j-1] + delta
    for i in range(1, m+1):
        for j in range(1, n+1):
            dp[i][j] = min(dp[i-1][j-1] + alphas[s[i-1]][t[j-1]], dp[i-1][j] + delta, dp[i][j-1] + delta)
    cost = dp[m][n]
    
    i, j = m, n
    s_alignment = ""
    t_alignment = ""
    # start from the bottom right corner and move towards the top left corner
    while i > 0 and j > 0:
        # case 1: diagonal move
        if dp[i][j] == dp[i-1][j-1] + alphas[s[i-1]][t[j-1]]:
            s_alignment = s[i-1] + s_alignment
            t_alignment = t[j-1] + t_alignment
            i -= 1
            j -= 1
        # case 2: vertical move
        elif dp[i][j] == dp[i-1][j] + delta:
            s_alignment = s[i-1] + s_alignment
            t_alignment = "_" + t_alignment
            i -= 1
        # case 3: horizontal move
        else:
            s_alignment = "_" + s_alignment
            t_alignment = t[j-1] + t_alignment
            j -= 1
    # if there are still elements left in the first string
    while i > 0:
        s_alignment = s[i-1] + s_alignment
        t_alignment = "_" + t_alignment
        i -= 1
    # if there are still elements left in the second string
    while j > 0:
        s_alignment = "_" + s_alignment
        t_alignment = t[j-1] + t_alignment
        j -= 1
    return cost, s_alignment, t_alignment



def main():
    # generate the strings from the input file
    s, t = generate_strings(sys.argv[1])
    # calculate the time taken by the basic algorithm
    time_taken, cost, s_alignment, t_alignment = time_wrapper(s, t, delta, alphas)
    # calculate the memory used by the basic algorithm
    memory_used = process_memory()
    # write the output to the output file
    with open(sys.argv[2], "w") as f:
        f.write(str(cost) + "\n")
        f.write(s_alignment + "\n")
        f.write(t_alignment + "\n")
        f.write(str(time_taken) + "\n")
        f.write(str(memory_used) + "\n")

if __name__ == "__main__":
    main()
