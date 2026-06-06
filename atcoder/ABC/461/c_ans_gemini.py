import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
5 3 3
1 30
1 40
1 50
2 10
3 20
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, k, m = map(int, input().split())
cv = []
for i in range(n):
    c, v = map(int, input().split())
    cv.append((c, v))

# 色(x[0])ではなく、価値(x[1])の降順にソートする
cv.sort(key=lambda x: x[1], reverse=True)

# ===== フェーズ1：とりあえず価値の高い順にK個選ぶ =====
selected = cv[:k]      # 選んだ宝石
unselected = cv[k:]    # 選ばれなかった宝石

# 選んだ宝石の色の種類数と、価値の合計を計算
color_counts = Counter([c for c, v in selected])
ans = sum([v for c, v in selected])
current_m = len(color_counts)

# ===== フェーズ2：色がM種類になるまで入れ替える =====
# 選んだ中から捨てる候補（価値が低い順＝後ろから）
drop_idx = k - 1 
# 選ばなかった中から入れる候補（価値が高い順＝前から）
add_idx = 0      

while current_m < m:
    # 捨てる宝石を探す：色がダブっていない（count == 1）なら捨てられないのでスキップ
    while color_counts[selected[drop_idx][0]] == 1:
        drop_idx -= 1
        
    # 入れる宝石を探す：既に選んでいる色（count > 0）なら新種にならないのでスキップ
    while color_counts[unselected[add_idx][0]] > 0:
        add_idx += 1
        
    # スワップする宝石が確定
    drop_color, drop_value = selected[drop_idx]
    add_color, add_value = unselected[add_idx]
    
    # 価値を更新（捨てた分を引いて、入れた分を足す）
    ans = ans - drop_value + add_value
    
    # 色のカウント状態を更新
    color_counts[drop_color] -= 1
    color_counts[add_color] += 1
    current_m += 1
    
    # 次の探索のためにインデックスをずらす
    drop_idx -= 1
    add_idx += 1

print(ans)