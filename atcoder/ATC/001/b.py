import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8 9
0 1 2
0 3 2
1 1 3
1 1 4
0 2 4
1 4 1
0 4 2
0 0 0
1 0 0
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n,q = map(int,input().split())
par = [i for i in range(n+1)]
def find(x):
    if par[x] == x:
        return x
    else:
        par[x] = find(par[x]) #経路圧縮
        return par[x]
def same(x,y):
    return find(x) == find(y)
def unite(x,y):
    x = find(x)
    y = find(y)
    if x == y:
        return 0
    par[x] = y

for i in range(q):
    p,a,b = map(int,input().split())
    if p == 0:
        unite(a,b)
    else:
        if same(a,b):
            print('Yes')
        else:
            print('No')
