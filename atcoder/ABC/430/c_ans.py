import io
import sys
from bisect import bisect_left

# ---- 入力（デバッグ用の埋め込み）----
_InPUT = """\
11 4 2
abbaaabaaba
"""
sys.stdin = io.StringIO(_InPUT)
# -------------------------------------

n, A, B = map(int, input().split())
s = input().strip()

# 累積和 XA, XB を作る（長さ n+1, 1-indexed 的に使う）
XA = [0]*(n+1)
XB = [0]*(n+1)
for i, ch in enumerate(s, 1):
    XA[i] = XA[i-1] + (ch == 'a')
    XB[i] = XB[i-1] + (ch == 'b')

ans = 0
for l in range(1, n+1):
    # a が A 個以上となる最小 r（存在しなければスキップ）
    if A == 0:
        r_a = l  # A=0 はどの区間でも満たすので r の最小は l
    else:
        targetA = XA[l-1] + A
        # r >= l を保証したいので lo=l を指定
        r_a = bisect_left(XA, targetA, lo=l)
        if r_a > n:
            continue  # これ以降の r では A 個以上にならない

    # b が B 個未満となる最大 r
    # XB[r] < XB[l-1] + B  <=>  bisect_left で targetB の直前が最大 r
    targetB = XB[l-1] + B
    r_b = bisect_left(XB, targetB, lo=l) - 1  # lo=l で r>=l を担保
    if r_b < r_a:
        continue

    ans += (r_b - r_a + 1)

print(ans)

"""
の解説（累積和＋二分探索）に沿って、コードを シンプル＆計算量 O(N log N) に組み直しました。
ポイントは：

XA[i]: 先頭から i 文字目までの 'a' の累積個数

XB[i]: 先頭から i 文字目までの 'b' の累積個数

各開始位置 l について

a_l: 「[l, r] に含まれる a が A 個以上」になる最小の r

累積和で言い換えると XA[r] >= XA[l-1] + A を満たす最小 r

b_l: 「[l, r] に含まれる b が B 個未満」になる最大の r

累積和で言い換えると XB[r] < XB[l-1] + B を満たす最大 r

bisect_left(XB, XB[l-1]+B) のひとつ前の添字

これで各 l の貢献は max(0, b_l - a_l + 1) です（a_l が存在しないときは 0）。
"""