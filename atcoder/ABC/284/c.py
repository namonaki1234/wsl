import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
4 6
1 2
1 3
1 4
2 3
2 4
3 4

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

n,m=map(int, input().split())
Uni=UnionFind(n)

for i in range(m):
    u,v=map(int, input().split())
    u-=1
    v-=1
    Uni.merge(u,v)

friends_group=Uni.groups()
# print(friends_group)
friends_size=[]
for fri in friends_group:
    friends_size.append(len(fri))

# print(max(friends_size))
print(len(friends_size))