import io
import sys

# 下記に標準入力を記載
_INPUT = """\
8

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# 目的はsetを使えるようになる
N=int(input())

able_num=set()

for a in range(2,10**5+10):
    for b in range(2,100):
        if a**b<=N:
            able_num.add(a**b)
        else:
            break

print(N-len(able_num))
