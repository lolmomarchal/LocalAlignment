# Helper functions
import argparse
import os.path
from typing import List, Dict, Iterable, Tuple
import matplotlib.pyplot as plt

def read_reads(file_path):
    reads = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if ">" not in line and line:
                reads.append(line)
    return reads
def write_sequences(output_file, alignment, aligned):
    if os.path.isfile(output_file):
        option = "w"
    else:
        option ="x"

    with open(output_file, option) as file:
        for keys in alignment.keys():
            if keys != "reverse":
                if aligned:
                    file.write(f">seq {keys}|score: {alignment[keys]["score"]}| length: {len(alignment[keys]["a"])}\n{alignment[keys]["a"]}\n{alignment[keys]["b"]}\n")
                else:
                    file.write(f">seq {keys}|score: {alignment[keys]["score"]}| length: {len(alignment[keys]["a"])}\n")


# Q1 local alignment
def locAL(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                   s: str, t: str) -> Tuple[int, str, str]:
    score= [[0 for j in range(len(t)+1)] for i in range(len(s)+1)]
    backtrack = [["" for j in range(len(t)+1)] for i in range(len(s)+1)]
    for i in range(len(s)+1):
        backtrack[i][0] = "0"
    for j in range(len(t)+1):
        backtrack[0][j] = "0"

    best_score = 0
    location_best_score = (0,0)
    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            if s[i-1] == t[j-1]:
                f = match_reward
            else:
                f = mismatch_penalty
            # obtain score
            score[i][j] = max(0,score[i-1][j-1] + f, score[i-1][j] +indel_penalty, score[i][j-1] + indel_penalty)
            if score[i][j] >= best_score:
                best_score = score[i][j]
                location_best_score = (i,j)

            # do backtrack
            if score[i][j] == score[i-1][j-1] + f:
                backtrack[i][j] = "↘"
            elif score[i][j] == score[i-1][j] + indel_penalty:
                backtrack[i][j] ="↓"
            elif score[i][j] == 0:
                backtrack[i][j] = "0"
            else:
                backtrack[i][j] = "→"
    i = location_best_score[0]
    j = location_best_score[1]
    seq_s = ""
    seq_t = ""
    # backtrack
    while i*j != 0 and backtrack[i][j] != "0":
        if backtrack[i][j] == "↘":
            i = i -1
            j = j-1
            seq_s += s[i]
            seq_t += t[j]
        elif backtrack[i][j] == "↓":
            i = i-1
            seq_s += s[i]
            seq_t += "-"
        elif backtrack[i][j] == "→":
            j = j-1
            seq_t += t[j]
            seq_s += "-"
        elif backtrack[i][j] == "0":
            break
    print(location_best_score)
    return [best_score, seq_s[::-1], seq_t[::-1]]

def parse_arguments_locAL():
    parser = argparse.ArgumentParser(description='Local alignment program for DNA sequences')
    parser.add_argument('seq_files', metavar='seq_files', type=str, nargs=1,
                        help='FASTA-formatted text file containing two sequences to be aligned')
    parser.add_argument('-m', '--match', type=int, required=True,
                        help='match reward')
    parser.add_argument('-s', '--mismatch', type=int, required=True,
                        help='mismatch penalty')
    parser.add_argument('-d', '--indel', type=int, required=True,
                        help='indel penalty')
    parser.add_argument('-a', '--alignment', action='store_true',
                        help='align sequences')
    parser.add_argument('-o', '--output', type=str, default='./Alignment.txt',
                        help='Output file path')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='plot distribution')
    return parser.parse_args()
def get_distribution(sequences):
    score_values = {}
    for keys in sequences.keys():
        if len(sequences[keys]["a"]) not in score_values.keys():
            score_values[len(sequences[keys]["a"]) ] = 1
        else:
            score_values[len(sequences[keys]["a"]) ] += 1
    return score_values

def main():
    args = parse_arguments_locAL()
    match_score = args.match
    mismatch_penalty = args.mismatch
    indel_penalty = args.indel
    alignment = args.alignment
    output_file = args.output
    plot = args.plot
    
    all_sequences = {}
    seq_count = 1
    for j in range(len(args.seq_files)):
        sequences = read_reads(args.seq_files[j])
        i = 0
        while i < len(sequences):
                s = sequences[i]
                t = sequences[i+1]

                # Align sequences
                score, aligned_s, aligned_t = locAL(match_score, mismatch_penalty, indel_penalty, s, t)
                all_sequences[seq_count] = {"a": aligned_s, "b": aligned_t, "score":score}
                i = i +2
                seq_count +=1
                # Output alignment results
        #first sort by score value
        sorted_aligned = dict(sorted(all_sequences.items(), key=lambda item: item[1]['score'], reverse=True))
        write_sequences(output_file, sorted_aligned, alignment)
    if plot:
        dist = get_distribution((all_sequences))
        x = list(dist.keys())
        y = list(dist.values())
        plt.bar(x, y)
        plt.xlabel("Length of alignment")
        plt.ylabel("Number of Pair Sequences")
        plt.title(f"Parameters -m {match_score} -s {mismatch_penalty} -d {indel_penalty}")
        plt.show()



if __name__ == "__main__":
    main()
