# Sequence Alignment Problem Memory Efficient Algorithm
# Hirshberg Alogirithm

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
    cost, s_alignment, t_alignment = memory_efficient_algorithm(s, t, delta, alphas)
    end_time = time.time()
    time_taken = (end_time - start_time)*1000
    return time_taken, cost, s_alignment, t_alignment

# generate the strings from the input file
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


def dp_calculator(s_, t_, m, n, mid, delta, alphas):
    dp_1 = [0] * (n + 1)
    dp_2 = [0] * (n + 1)

    for j in range(n + 1):
        dp_1[j] = j * delta

    
    # if the length of s_left is 0, then then all zeros is the cost we want
    if len(s_) == 0:
        pass
    else:
        for i in range(1, mid + 1):
            for j in range(n + 1):
                # case where j = 0
                if j == 0:
                    dp_2[j] = i * delta
                else:
                    dp_2[j] = min(dp_2[j - 1] + delta, dp_1[j] + delta, dp_1[j - 1] + alphas[s_[i - 1]][t_[j - 1]])
            
            # set dp_1 to dp_2
            dp_1 = dp_2[:]

    return dp_2



# implementation of the memory efficient algorithm
def memory_efficient_algorithm(s, t, delta, alphas):
    m, n = len(s), len(t)

    # Base cases
    if m == 0:
        gaps = "_" * n
        return n * delta, gaps, t
    if n == 0:
        gaps = "_" * m
        return m * delta, s, gaps
    # both are single characters
    if m == 1 and n == 1:
        # same character
        if s == t:
            return 0, s, t
        # different but mismatch is cheaper
        elif alphas[s][t] < 2 * delta:
            return alphas[s][t], s, t
        # different and mismatch is more expensive / gaps are cheaper
        else:
            return 2 * delta, "_" + s, t + "_"
    # just s is a single character
    if m == 1:
        # cases:
        # same char is within
        if s in t:
            index = t.index(s)
        # there is a cheaper mismatch
        elif s == "A" and "G" in t:
            index = t.index("G")
        elif s == "C" and "T" in t:
            index = t.index("T")
        elif s == "G" and "A" in t:
            index = t.index("A")
        elif s == "T" and "C" in t:
            index = t.index("C")
         # all gaps is cheaper
        else:
            return n * delta, s + ("_" * n), "_" + t
        
        return alphas[s][t[index]] + ((n - 1) * delta), (index * "_") + s + ((n - index - 1) * "_"), t




    # now all normal cases when m > 1 and n > 1
    # Split s into s_left and s_right
    mid = m // 2
    s_left, s_right = s[:mid], s[mid:]


    dp_left = dp_calculator(s_left, t, m, n, mid, delta, alphas)



    # make a copy of s_right and reverse it
    s_right_rev = s_right[::-1]
    # make a copy of t and reverse it
    t_rev = t[::-1]


    # fix mid for m--1
    dp_right = dp_calculator(s_right_rev, t_rev, m, n, m - mid, delta, alphas)



    # for now make a new array
    final_cost = [0] * (n + 1)

    # Find the optimal place to split t
    min_cost = float('inf')
    for j in range(n + 1):
        cost = dp_left[j] + dp_right[n - j]
        final_cost[j] = cost

        if cost < min_cost:
            min_cost = cost
            t_split = j


    # Split t into t_left and t_right
    t_left, t_right = t[:t_split], t[t_split:]

    # Recursively apply the same process to s_left and t_left, and s_right and t_right
    cost_left, s_left, t_left = memory_efficient_algorithm(s_left, t_left, delta, alphas)
    cost_right, s_right, t_right = memory_efficient_algorithm(s_right, t_right, delta, alphas)

    # return: cost, s_alignment, t_alignment
    return cost_left + cost_right, s_left + s_right, t_left + t_right



def main():
    # generate the strings from the input file
    s, t = generate_strings(sys.argv[1])

    # for testing purposes
    # s, t = generate_strings("SampleTestCases/" + "input5.txt")

    # calculate the time taken by the basic algorithm
    time_taken, cost, s_alignment, t_alignment = time_wrapper(s, t, delta, alphas)
    # calculate the memory used by the basic algorithm
    memory_used = process_memory()


    # print the final results
    # print("\n-----Final Results-----\n")
    # print("cost: ", cost)
    # print("s_alignment: ", s_alignment)
    # print("t_alignment: ", t_alignment)

    # write the output to the output file
    with open(sys.argv[2], "w") as f:
        f.write(str(cost) + "\n")
        f.write(s_alignment + "\n")
        f.write(t_alignment + "\n")
        f.write(str(time_taken) + "\n")
        f.write(str(memory_used) + "\n")

if __name__ == "__main__":
    main()

