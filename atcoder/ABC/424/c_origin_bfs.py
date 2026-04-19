import io
import sys
from collections import deque

# 下記に標準入力を記載 今回のエラーの原因がわかるような例を考えた
_InPUT = """\
6
0 0
1 3
3 2
5 5
4 6
6 4
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n = int(input())

# ab = [int(x) for x in input().split()]

"""
以下はgraph[x] に
「スキル x を持っていると取れるかもしれないスキル」
を入れるためのリスト
いつもは1 3と合ったら1と3の間に辺があったが、今回は一行目ならスキル1という対応なので、
iをgragh[x]の矢印の先として設定する
"""
gragh = [[] for _ in range(n+1)]
visited = [False] * (n+1)
q = deque()

for i in range(1,n+1):
    a,b = map(int,input().split())

    if a == 0 and b == 0:
        visited[i] = True
        q.append(i)
    else:
        gragh[a].append(i)
        gragh[b].append(i)

while q:
    c = q.popleft()

    for d in gragh[c]:
        if visited[d]:
            continue
        visited[d] = True
        q.append(d)

print(visited.count(True))

# for _ in range(n):
#     a,b = map(int,input().split())
#     gragh[a].append(b)
#     gragh[b].append(a)

# q = deque()
# q.append(1)
# visited = [False] * (n+1)
# visited[1] = True
# while q:
#     c = q.popleft()

#     for d in gragh[c]:
#         if visited[d]:
#             continue
#         visited[d] = True
#         q.append(d)

# print(visited.count(True))