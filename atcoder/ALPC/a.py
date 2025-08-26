import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
4 7
1 0 1
0 0 1
0 2 3
1 0 1
1 1 2
0 0 2
1 1 3
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

class UnionFind:
    def __init__(self, n):
        self.n=n
        self.parent_size=[-1]*n
    def merge(self, a, b):
        x, y=self.leader(a), self.leader(b)
        if x == y: return
        if abs(self.parent_size[x])<abs(self.parent_size[y]): x, y=y, x
        self.parent_size[x] += self.parent_size[y]
        self.parent_size[y]=x
        return
    def same(self, a, b):
        return self.leader(a) == self.leader(b)
    def leader(self, a):
        if self.parent_size[a]<0: return a
        self.parent_size[a]=self.leader(self.parent_size[a])
        return self.parent_size[a]
    def size(self, a):
        return abs(self.parent_size[self.leader(a)])
    def groups(self):
        result=[[] for _ in range(self.n)]
        for i in range(self.n):
            result[self.leader(i)].append(i)
        return [r for r in result if r != []]

n,q = map(int,input().split())
uf=UnionFind(n)

# for i in range(q):
#     u,v=map(int, input().split())
#     u-=1
#     v-=1
#     Uni.merge(u,v)

# friends_group=Uni.groups()
# # print(friends_group)
# friends_size=[]
# for fri in friends_group:
#     friends_size.append(len(fri))

# print(max(friends_size))
# print(len(friends_size))

for i in range(q):
    t,a,b = map(int,input().split())
    if t == 0:
        uf.merge(a,b)
    else:
        if uf.same(a,b):
            print(1)
        else:
            print(0)


# def find(x):
#     if par[x] == x:
#         return x
#     else:
#         par[x] = find(par[x]) #経路圧縮
#         return par[x]
# def same(x,y):
#     return find(x) == find(y)
# def unite(x,y):
#     x = find(x)
#     y = find(y)
#     if x == y:
#         return 0
#     par[x] = y