import io
import sys

# 下記に標準入力を記載
_INPUT = """\
5
3 9 5 3 1
"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
# ワンパス法（last_pos）：最短・最速・省メモリ
"""
さらに軽い：ワンパス “最後に見た位置” 法（おすすめ）

実はこの問題は 右から順に見て、同じ値を最後に見た位置との距離だけ比較すれば十分です。
値 x を位置 i で見たとき、直前の出現位置が j なら、重複を含む最短の区間の一つは [j, i]（長さ i - j + 1）。これを全値で最小化します。

公式の解答コード：値ごとに出現位置を溜めて、隣り合う出現の差 + 1 の最小を取る。O(N)。

尺取り法：重複がある状態を保ちながら左端を詰めて最短化。O(N)。

ワンパス（last_pos）：直前の出現位置だけ覚える超シンプル解。O(N)。
"""

N = int(input())
A = list(map(int, input().split()))

last = {}        # 値 -> 直前に現れた添字
ans = N + 1

for i, x in enumerate(A):
    if x in last:
        ans = min(ans, i - last[x] + 1)
    last[x] = i

print(-1 if ans == N + 1 else ans)
