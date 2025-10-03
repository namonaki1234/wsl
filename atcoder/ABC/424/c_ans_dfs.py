import io
import sys

# 下記に標準入力を記載 今回のエラーの原因がわかるような例を考えた
_InPUT = """\
6
6 1
2 3
3 1
0 0
4 6
0 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

sys.setrecursionlimit(10**9)

N=int(input())
G=[[] for _ in range(N+1)]

for i in range(1,N+1):
  a,b=map(int,input().split())
  G[a].append(i)
  G[b].append(i)

print(G)
ok=[0]*(N+1)
ok[0]=1

def dfs(v):
  ok[v]=1
  for vv in G[v]:
    if not ok[vv]:
      dfs(vv)

dfs(0)

print(sum(ok)-1)
