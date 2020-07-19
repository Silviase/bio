import sys
from typing import Sequence
import copy
import tkinter as tk

sys.setrecursionlimit(10000)

__gap_penalty__ = -1
__match_score__ = 1
__unmatch_score__ = -3
memory_score = {0: 0}
memory_trace = {}

__max_length__ = 10


# 何文字目まで見たかをHashにする
# List[int]だとhash化できないのでDictionaryのKeyに出来ない
def hash_code(inspect: Sequence[int]) -> int:
    res = 0
    for i in range(len(inspect)):
        res = res * (__max_length__ + 1) + inspect[i]
    return res


def print_result(frame: tk.Frame):
    global test_ins
    test_trace_hashed = trace_back(hash_code(test_ins))
    # print(test_trace_hashed)
    test_trace_unhashed = unhash_trace(test_trace_hashed, len(test_ins))
    # print(test_trace_unhashed)
    for i in range(len(test_trace_unhashed[0])):
        aligned_dna = ""
        for j in range(1, len(test_trace_unhashed)):
            if test_trace_unhashed[j][i] > test_trace_unhashed[j - 1][i]:
                aligned_dna += test_nuc[i][test_trace_unhashed[j - 1][i]]
                # print(test_nuc[i][test_trace_unhashed[j - 1][i]], end="")
            else:
                aligned_dna += "_"
                # print("_", end="")
        aligned_dna_label = tk.Label(frame, text=aligned_dna, font="Consolas")
        aligned_dna_label.grid(row=i + 3)
        # print()
    frame.grid(row=2, sticky=tk.W + tk.E)
    return


def do_alignment():
    global test_nuc, test_ins, memory_trace
    memory_trace = {}
    max_score = dynamic_alignment(test_nuc, test_ins)

    # 結果を示すフレーム
    result_frame = tk.Frame(root, width=400, height=100, bg="green")
    # ラベルの作成
    label = tk.Label(result_frame, text="Score->" + str(max_score))
    label.grid(row=0)
    alignment_result = tk.Label(result_frame, text="Result below")
    alignment_result.grid(row=1, sticky=tk.N)

    dynamic_alignment(test_nuc, test_ins)
    # print(score)

    print_result(result_frame)

    return


def dynamic_alignment(nucleotide_sequence: Sequence[str], inspect: Sequence[int]) -> int:
    # メモ化再帰して計算量を削減する
    hashed_inspect = hash_code(inspect)
    if hashed_inspect in memory_score:
        return memory_score[hashed_inspect]
    else:
        pass
        # print(str(inspect) + " not found")

    max_score = -1e9
    # 各bitについて計算を行うbit全探索をする. O(2^N)
    for bit in range(1, 1 << len(inspect)):
        f = False
        trace_back = copy.deepcopy(inspect)
        for i in range(len(inspect)):
            if trace_back[i] > 0 and (bit & (1 << i)) > 0:
                trace_back[i] -= 1
            elif trace_back[i] == 0 and (bit & (1 << i)) > 0:
                f = True
                break
        if f:
            continue

        basic_score = dynamic_alignment(nucleotide_sequence, trace_back)
        for i in range(0, len(inspect)):
            for j in range(i + 1, len(inspect)):
                if (bit & (1 << i)) > 0 and (bit & (1 << j)) > 0:
                    if nucleotide_sequence[i][inspect[i] - 1] == nucleotide_sequence[j][inspect[j] - 1]:
                        basic_score += __match_score__
                    else:
                        basic_score += __unmatch_score__
                elif (bit & (1 << i)) == 0 and (bit & (1 << j)) == 0:
                    pass
                else:
                    basic_score += __gap_penalty__

        if basic_score > max_score:
            # print(str(trace_back) + " to " + str(inspect) + " has updated, score is " + str(basic_score))
            max_score = basic_score
            memory_trace[hashed_inspect] = hash_code(trace_back)

    memory_score[hashed_inspect] = max_score

    with open('output.txt', 'a', encoding='utf-8') as f:
        f.write("put " + str(inspect) + ", score = " + str(max_score) + "\n")

    return memory_score[hashed_inspect]


# memory_traceによって示された文字列のマッチと対応するHashの配列を返す
def trace_back(to_trace_back: int) -> Sequence[int]:
    res = [to_trace_back]
    key = to_trace_back
    while key in memory_trace:
        next_key = memory_trace.get(key)
        res.append(next_key)
        key = next_key
    return res


def unhash_trace(hashed_trace: Sequence[int], length: int) -> Sequence[Sequence[int]]:
    res = []
    for i in range(len(hashed_trace)):
        hashed = hashed_trace[i]
        push = []
        for j in range(length):
            push.insert(0, hashed % (__max_length__ + 1))
            hashed //= (__max_length__ + 1)
        res.insert(0, push)
    return res


def add_dna():
    new_dna = add_text_box.get()
    global test_nuc, test_ins
    test_nuc.append(new_dna)
    test_ins.append(len(new_dna))
    current_nuc_label["text"] = str(test_nuc)
    return


test_nuc = ["CTAGGAG", "CTGGAAG", "CGAGGAT", "ATAGGAG", ]
test_ins = [len(test_nuc[i]) for i in range(len(test_nuc))]
score = 0

# root の雛形を作成
root = tk.Tk()
root.title("Alignment")
root.geometry("400x400")

add_frame = tk.Frame(root, width=400, height=100, bg="cyan")

add_dna_button = tk.Button(add_frame, text="Add", command=add_dna)
add_text_box = tk.Entry(add_frame)
add_text_box.grid(row=0, column=0)
add_dna_button.grid(row=0, column=1)
add_frame.grid(row=0, column=0, sticky=tk.W + tk.E)

# 二段目 : 現在の様子を示す
current_frame = tk.Frame(root, width=400, height=100, bg="navy")
current_nuc_label = tk.Label(current_frame, text=str(test_nuc))
current_nuc_label.grid()
button = tk.Button(current_frame, text="DP Alignment!", command=do_alignment)
button.grid()
current_frame.grid(row=1, column=0)

# ----------------
if __name__ == '__main__':
    root.mainloop()
