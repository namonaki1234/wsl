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

S = S.replace("o", "1").replace("?", "0").replace(".", "-1")
S = list(map(int, S))

Counter(S)

for i in S:
    if i == 0:
        if S.count(0) == K:
            S.replace(0, 1)
        else:













