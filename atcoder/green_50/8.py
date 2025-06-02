import io
import sys

# 下記に標準入力を記載
_INPUT = """\
314 2

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# 目的は関数実装のマスター
def g1(x):
    x = str(x)
    x = list(x)
    x.sort(reverse=True)
    x = ''.join(x)
    return int(x)

def g2(x):
    x = str(x)
    x = list(x)
    x.sort()
    x = ''.join(x)
    return int(x)

def f(x):
    return g1(x) - g2(x)

N,K = map(int, input().split())
a = N
for i in range(K):
    a = f(a)

print(a)