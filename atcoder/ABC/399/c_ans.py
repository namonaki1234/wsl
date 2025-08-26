import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
10 10
7 9
4 6
6 10
2 5
5 6
5 9
6 8
4 8
1 5
1 4
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

class UnionFind:
    def __init__(self, n):
        self.n=n
        self.parent_size=[-1]*n
    def merge(self, a, b):
        x, y = self.leader(a), self.leader(b)
        if x == y:
            return False  # 併合なし
        # union by size（サイズ大に小をつける）
        if abs(self.parent_size[x]) < abs(self.parent_size[y]):
            x, y = y, x
        self.parent_size[x] += self.parent_size[y]  # 負のサイズを加算
        self.parent_size[y] = x                    # y の親を x に
        return True  # 併合が起きた
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

n,m = map(int, input().split())
# edges = [tuple(map(int, input().split())) for _ in range(m)]

# print(edges)

uf = UnionFind(n)
k = n  # 最初は N 個の連結成分

for _ in range(m):
    u,v=map(int, input().split())
    u-=1
    v-=1
    if uf.merge(u,v):
        k -= 1

# 森にするために削除すべき最小本数 = M - (N - K)
print(m - (n - k))