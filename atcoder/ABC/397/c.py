import io
import sys
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
10
2 5 6 5 2 1 7 9 7 2

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N = int(input())
A = list(map(int, input().split()))

# split_A = list(chunked(A, 3))
# print(split_A)
# n = 2
# result = [A[i:i + n] for i in range(0, len(A), n)]
# print(result)  # [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10]]
# print(A)

def split_list(lst, index):
    first_half = lst[:index]
    second_half = lst[index:]
    return first_half, second_half

ans = 0
for i in range(N-1):
    first = A[:i+1]
    second = A[i+1:]
    set_first = set(first)
    set_second = set(second)
    # print(len(set_first)+len(set_second))
    ans = max(ans,len(set_first)+len(set_second))
    if ans == N:
        break
print(ans)











