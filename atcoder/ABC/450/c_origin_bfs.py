import io
import sys
from collections import Counter,deque
from atcoder.dsu import DSU
# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)


# 下記に標準入力を記載
_InPUT = """\
5 15
##########..###
#...#######.###
####....###..##
######.########
########....###
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
"""
処理の流れを短く言うと
1.全マスを順番に見る
2.未訪問の白マスを見つけたら BFS
3.その BFS で同じ白領域を全部訪問済みにする
4.その領域が外周に触れていなければ ans += 1
"""

h,w = map(int,input().split())

grid = [input().strip() for _ in range(h)]

dx = [0,1,0,-1]
dy = [1,0,-1,0]

visited = [[False]*w for _ in range(h)]
# visited[0][0] = True

def bfs(sx,sy):
    #この問題では、未訪問の白マス (i, j) のたびに、BFS を始める必要があるので、q は bfs の中で毎回作るべき
    q = deque()
    q.append((sx,sy))
    visited[sx][sy] = True
    #返り値となるoutは関数内で宣言＆代入をしておく必要がある？
    out = False
    while q:
        x,y = q.popleft()

        #外周にあるマスを含んでいてたらout = True
        if x == 0 or x == h - 1 or y == 0 or y == w -1:
            out = True
        for i in range(4):
            nx = x + dx[i]
            ny = y + dy[i]
        
            #場外
            if nx < 0 or nx >= h or ny < 0 or ny >= w:
                continue

            #壁,#は飛ばす
            if grid[nx][ny] == "#":
                continue
            if visited[nx][ny]:
                continue
            visited[nx][ny] = True
            q.append((nx,ny))
    return out

ans = 0
for i in range(h):
    for j in range(w):
        if grid[i][j] == "." and not visited[i][j]:
            out = bfs(i,j)

            if not out:
                ans += 1

print(ans)