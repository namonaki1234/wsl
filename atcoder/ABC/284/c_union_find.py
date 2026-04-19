import io
import sys
from atcoder.dsu import DSU

# 下記に標準入力を記載
_INPUT = """\
5 3
1 2
1 3
4 5
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# Union-Find は、「同じグループかどうか」を管理するデータ構造
#1-indexでやる
n,m = map(int,input().split())

uf = DSU(n+1)

for _ in range(m):
    a,b = map(int,input().split())
    uf.merge(a,b)

ans = 0
for i in range(1,n+1):
    #頂点 i が属しているグループの代表元を返す。
    #たとえば0,1,2 が同じグループ、3,4 が同じグループなら、0,1,2 の leader は同じ値になり、3,4 の leader も同じ値になる。
    if uf.leader(i) == i:
        ans += 1

print(ans)