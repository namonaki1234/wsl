import io
import sys

_InPUT = """\
5
1 1
4 2
2 3
5 5
3 4
"""
sys.stdin = io.StringIO(_InPUT)

n = int(input())
xy = []
for i in range(n):
    x,y = map(int,input().split())
    xy.append((x,y))

ans = 0
# すべての点 i について調べる
for i in range(n):
    x_i, y_i = xy[i]
    is_ok = True # 最初は「条件を満たしている」と仮定する
    
    # 他のすべての点 j と比較する
    for j in range(n):
        if i == j: # 自分自身との比較はスキップ
            continue
            
        x_j, y_j = xy[j]
        
        # 点 j が点 i の「左」かつ「下」にある場合（＝長方形の中に入ってしまった場合）
        if x_j < x_i and y_j < y_i:
            is_ok = False # 条件を満たさないことが確定
            break         # これ以上調べる必要がないので内側のループを抜ける
            
    # 最後まで長方形の中に点が入らなかったらカウント
    if is_ok:
        ans += 1
            
print(ans)