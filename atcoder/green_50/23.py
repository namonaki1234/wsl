import io
import sys

# 下記に標準入力を記載
_INPUT = """\
6
2 4 4 9 4 9

"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

# 制約を見てpypyなら間に合うか？を判断できるようになる
N=int(input())
A=list(map(int, input().split()))

ans=-10**15
for l in range(N):
    A_min=A[l]
    for r in range(l,N):
        A_min=min(A_min,A[r]) #みかんの個数x
        ans=max(ans,A_min*(r-l+1))

print(ans)
