import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8
3 3
2 7 1
4 5
2 9 16 14
4 1
1 2 3 4
5 2
5 6 7 8 9
2 10
7 9
1 1000000000
1
1 371
1000000000
3 6
1 3 5
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

t = int(input())
output = []
for _ in range(t):
    n,k = map(int, input().split())
    a = list(map(int, input().split()))
    # print(n,k,a)
    m = k + 1

    # 最終値: ai + k * (ai % (k+1))
    result = [ai + k * (ai % m) for ai in a]
    output.append(" ".join(map(str, result)))

print("\n".join(output))