import io
import sys
from collections import Counter,OrderedDict


# 下記に標準入力を記載
_INPUT = """\
7 3
.o???o.


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N, K = map(int, input().split())
S = input().rstrip()

S = S.replace("o", "1").replace("?", "0").replace(".", "2")
S = list(map(int, S))

Counter(S)

for i in range(N):
    if S[i] == 0:
        for j in range(1<<len(S)):
            if S[j] == 0:
                S[j] = 1
                break
            if S[j] == 0:
                S[j] = 1
                break
        
    












