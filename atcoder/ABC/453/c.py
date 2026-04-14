import io
import sys
import itertools
# 下記に標準入力を記載
_InPUT = """\
20
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 初期座標0.5
# i回目に正負どちらかにl_i移動
# 最大0を何回通れるか
# n<=20の全探索、分岐の数は2**n

n = int(input())
l = [int(l_i) for l_i in input().split( )]

# print(l)

start_pos = 0.5
max_crossings = 0

for pattern in itertools.product([1, -1], repeat=len(l)):
    # print(pattern)
    current_pos = start_pos
    crossing = 0 # 初期化
    for i in range(n):
        next_pos = current_pos + l[i]*pattern[i]
        if current_pos*next_pos < 0:
            crossing += 1
        current_pos = next_pos
        
    max_crossings = max(max_crossings,crossing)

print(max_crossings)