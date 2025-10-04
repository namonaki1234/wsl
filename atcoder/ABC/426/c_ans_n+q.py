import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_InPUT = """\
8 5
2 6
3 5
1 7
5 7
7 8

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, q = map(int,input().split())

# pc[v] = 「現在バージョン v に存在する PC の台数」
pc = [0] * (n + 1)
for v in range(1, n + 1):
    pc[v] = 1  # 初期状態は v=1..N に 1 台ずつ
pc[0] = 0     # 使わない番地（ガード）

o = 1  # 最古のバージョン（まだ未処理な最小バージョン）= O
out = []

for _ in range(q):
    x, y = map(int,input().split())
    res = 0
    # O ≤ x の範囲だけ、古い順にまとめて y へ「移送」する
    while o <= x:
        res += pc[o]      # 今回アップグレードされる台数に加算
        pc[y] += pc[o]    # そのぶんを y バージョンへ集約
        pc[o] = 0         # 元の古いバージョンは空になる
        o += 1            # 次の最古へ

    out.append(str(res))

print("\n".join(out))
