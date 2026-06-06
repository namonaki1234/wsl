import io
import sys
from collections import Counter

# 下記に標準入力を記載
_INPUT = """\
5 3 3
1 30
1 40
1 50
2 10
3 20
"""
sys.stdin = io.StringIO(_INPUT)

# ここからコードを記載

n, k, m = map(int, input().split())

cv = []
for _ in range(n):
    c, v = map(int, input().split())
    cv.append((c, v))

# 解法2では「価値の高い順」に見る
cv.sort(key=lambda x: x[1], reverse=True)

# まず色を無視して，価値が高い順にK個選ぶ
selected = cv[:k]
rest = cv[k:]

# 現在の価値の合計
ans = sum(v for c, v in selected)

# 選んだ宝石の色ごとの個数
color_count = Counter(c for c, v in selected)

# 現在選んでいる色の種類数
kind_cnt = len(color_count)

# すでに条件を満たしているなら，そのまま答え
if kind_cnt >= m:
    print(ans)
    sys.exit()

# 捨てても色数が減らない宝石を集める
# つまり，同じ色が複数ある中の「余り」の宝石
remove_values = []

# 低価値のものから捨てたいので，価値の昇順に見る
tmp_count = color_count.copy()

for c, v in sorted(selected, key=lambda x: x[1]):
    # その色が2個以上あるなら，この宝石を捨てても色数は減らない
    if tmp_count[c] > 1:
        remove_values.append(v)
        tmp_count[c] -= 1

# 追加候補を集める
# まだ選んでいない色の中で，価値が高いものから使う
add_values = []

used_colors = set(color_count.keys())

for c, v in rest:
    # まだ選んでいない色なら，その色を増やせる
    if c not in used_colors:
        add_values.append(v)
        used_colors.add(c)

# 足りない色数だけ交換する
need = m - kind_cnt

for i in range(need):
    ans -= remove_values[i]
    ans += add_values[i]

print(ans)