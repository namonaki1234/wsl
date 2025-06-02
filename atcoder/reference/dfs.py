import io
import sys

# 下記に標準入力を記載
_INPUT = """\
4 5
s####
....#
#.###
#...g
"""

sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# https://atcoder.jp/contests/atc001/tasks/dfs_a

def field_print(field):
    print("-------------------------")
    for i in field:
        print(i)

H,W = map(int, input().split())

field = []
for i in range(H):
    row = list(input())
    for j in range(W):
        if row[j]=="s":
            s_idx = (i,j)
    field.append(row)
print(s_idx)
field_print(field)

#N,S,E,W
direc = [(-1,0),(0,1),(1,0),(0,-1)]

def dfs(h,w):
    print(h,w)
    field[h][w] = "P"
    field_print(field)
    for d in direc:
        next_h = h+d[0]
        next_w = w+d[1]
        if next_h >= H or next_w >= W:
            continue
        if next_h < 0 or next_w < 0:
            continue
        next_point = field[next_h][next_w]
        if next_point == "#":
            continue
        if next_point == "-":
            continue
        if next_point == "g":
            print("Yes")
            exit()
        field[h][w] = "-"
        dfs(next_h,next_w)
    field[h][w] = "-"

dfs(s_idx[0], s_idx[1])
print("No")



        
        
        


