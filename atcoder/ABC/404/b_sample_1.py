import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
3 3
.#.
###
.#.


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

H,W = map(int,input().split())
paint = [['.']*(W+2)]
for i in range(H):
    temp = ['.'] + list(input()) + ['.']    
    paint.append(temp)
paint.append(['.']*(W+2))

for i in range(1,H+1):
    for j in range(1,W+1):
        if paint[i][j] == '#' and paint[i-1][j] == '.' and paint[i+1][j] == '.' and paint[i][j-1] == '.' and paint[i][j+1] == '.':
            print('No')
            exit()
print('Yes')





    


    











