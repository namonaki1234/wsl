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
N, K, M = map(int, input().split())

# t[color] に，その色の宝石の価値をまとめる
t = [[] for _ in range(N + 1)]

# 以下のように配列をあたかも数値のキーを持つ辞書の様に利用する
# t[1] = 色1の宝石価値リスト
# t[2] = 色2の宝石価値リスト
# t[3] = 色3の宝石価値リスト
for _ in range(N):
    C, V = map(int, input().split())
    t[C].append(V)

# top  : 各色の代表宝石，つまりその色で一番価値が高い宝石
# tail : 代表以外の宝石
top = []
tail = []


# ソートした後の箱（例：[50, 40, 30]）から、
# 一番大きい代表（values[0] つまり 50）を取り除いた残りすべて（values[1:] つまり [40, 30]）を、
# 一発で tail（落選者リスト）に合流させています。

for values in t:
    # その色の宝石が1つもないならスキップ
    if len(values) == 0:
        continue

    # 価値が高い順に並べる
    values.sort(reverse=True)

    # 一番価値が高いものをその色の代表にする
    top.append(values[0])

    # 残りは代表以外の候補に入れる
    tail += values[1:]

# 代表宝石を価値が高い順に並べる
top.sort(reverse=True)

# 代表のうち，上位M個は必ず選ぶ
selected_top = top[:M]

# M個に入らなかった代表は，普通の候補に戻す
tail += top[M:]

# 残り候補を価値が高い順に並べる
tail.sort(reverse=True)

# 代表M個を選んだので，あとK-M個をtailから選べばよい
answer = sum(selected_top) + sum(tail[:K - M])

print(answer)