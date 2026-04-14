import io
import sys
import itertools
# 下記に標準入力を記載
_InPUT = """\
5
2 5 2 2 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 初期座標0.5
# i回目に正負どちらかにl_i移動
# 最大0を何回通れるか
# n<=20の全探索、分岐の数は2**n

n = int(input())
l = [int(l_i) for l_i in input().split( )]

start_pos = 0.5

# 2分岐するから、関数の中にその分2個の同じ名前の関数を書く、そして最後の分岐の先まで着たらそこでは分岐する前にreturn 0をして終了にする
def get_max_crossings(index, current_pos):
    if index == n:
        return 0
    
    # 左に行くパターン
    pos_l = current_pos + l[index] 
    cross_l = 1 if current_pos * pos_l < 0 else 0
    total_l = cross_l + get_max_crossings(index+1,pos_l)
    # 右に行くパターン
    pos_r = current_pos - l[index] 
    cross_r = 1 if current_pos * pos_r < 0 else 0
    total_r = cross_r + get_max_crossings(index+1,pos_r)

    return max(total_l,total_r)
# 最初の自分をスタートさせる
print(get_max_crossings(0, 0.5))