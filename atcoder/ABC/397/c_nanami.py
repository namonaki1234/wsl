import io
import sys
from collections import Counter,OrderedDict


# 下記に標準入力を記載
_INPUT = """\
5
3 1 4 1 5

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N = int(input())
A = list(map(int, input().split()))


ans = 0
for i in range(N-1):
    # first, second = split_list(A, i+1)
    first = A[:i+1]
    second = A[i+1:]
    counter_first = Counter(first)
    counter_second = Counter(second)
    # print(len(counter_first)+len(counter_second))
    ans = max(ans,len(counter_first)+len(counter_second))
    if ans == N:
        break
print(ans)














