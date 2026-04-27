import io
import sys
import itertools
from atcoder.dsu import DSU
# 再帰の上限を増やす（Pythonでは必須のおまじない）
sys.setrecursionlimit(10**6)

# 下記に標準入力を記載
_InPUT = """\
3 2
.#
#.
##
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

# h,w = map(int, input().split())
# s = [input() for _ in range(h)]

from functools import lru_cache

H, W = map(int, input().split())
S = [input() for _ in range(H)]


@lru_cache(maxsize=None)
def dfs(h1, h2, w1, w2):
    """
    長方形 (h1, h2, w1, w2) が点対称か判定するDFS。
    h1, h2, w1, w2 はすべて 0-indexed。
    h2, w2 は「含む」。
    """

    # 左上と右下が違うなら点対称ではない
    if S[h1][w1] != S[h2][w2]:
        return False

    # 右上と左下が違うなら点対称ではない
    if S[h1][w2] != S[h2][w1]:
        return False

    # 上下を1行ずつ削った内側を見る
    if h1 + 1 < h2:
        if not dfs(h1 + 1, h2 - 1, w1, w2):
            return False

    # 左右を1列ずつ削った内側を見る
    if w1 + 1 < w2:
        if not dfs(h1, h2, w1 + 1, w2 - 1):
            return False

    return True


ans = 0

for h1 in range(H):
    for h2 in range(h1, H):
        for w1 in range(W):
            for w2 in range(w1, W):
                if dfs(h1, h2, w1, w2):
                    ans += 1

print(ans)