import io
import sys
import numpy as np
import pprint

# 下記に標準入力を記載
_INPUT = """\
2

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
# NxNの行列を作成
grid =[]

for i in range(1,N+1):
    if i == 1 or i == N:
        grid.append(['#']*N)
    
    elif N%2 == 0 and i >= int(N/2)+1:
        for k in range(1, int(N/2)+1):
            if  i == int(N/2)+k:
                grid.append(grid[int(N/2)-k])
                continue
    
    elif N%2 != 0 and i >= int(N/2)+2:
        for k in range(1, int(N/2)+1):
            if  i == int(N/2)+k:
                grid.append(grid[int(N/2)-k+1])
                continue

    # 偶数
    elif i%2 == 0 and i < int(N/2)+1:
        grid.append(['.']*N)
        grid[i-1][0] = '#'
        grid[i-1][N-1] = '#'
        if i == 2 or i == N-1:
            continue
        grid[i-1][i-2] = '#'
        grid[i-1][N-(i-1)] = '#'
    
    elif N%2 != 0 and i == int(N/2)+1:
        grid.append(['#']*N)
        for q in range(1, N+1):
            if q%2 == 0:
                grid[i-1][q-1] = '.'
            

    # 奇数
    elif i%2 != 0:
        grid.append(['#']*N)
        if i == 2 or i == N-1:
            continue
        grid[i-1][1] = '.'
        grid[i-1][N-2] = '.'
        grid[i-1][i-2] = '.'
        grid[i-1][N-(i-1)] = '.'


# pprint.pprint(grid,indent = 0, width=40)
for i in range(N):
    print(''.join(grid[i]))
    










