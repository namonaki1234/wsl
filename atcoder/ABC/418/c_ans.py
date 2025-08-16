import sys
import io

# ===== サンプル入力（提出時はこの3行を消す） =====
_InPUT = """\
5 3
13 13 13 13 2
5
12
13
"""
sys.stdin = io.StringIO(_InPUT)
# ==============================================

# 入力
n, q = map(int, input().split())
a = list(map(int, input().split()))
bs = [int(input()) for _ in range(q)]

# 配列サイズは必要最小限（max(a) と max(b)-1 の大きい方まで）
max_a = max(a) if n > 0 else 0
max_bm1 = max((b - 1) for b in bs) if q > 0 else 0
M = max(0, max(max_a, max_bm1))

accum_sum = [0] * (M + 1)
accum_cnt = [0] * (M + 1)

# 値ごとの合計と個数を加算
for v in a:
    # v は 0〜M のはず（b-1 の最大と max(a) に合わせたため）
    accum_sum[v] += v
    accum_cnt[v] += 1

# 累積和
for i in range(1, M + 1):
    accum_sum[i] += accum_sum[i - 1]
    accum_cnt[i] += accum_cnt[i - 1]

# クエリ処理
out = []
for b in bs:
    idx = b - 1
    if idx < 0:
        cnt_bm1 = 0
        sum_bm1 = 0
        bm1 = -1
    else:
        if idx > M:
            # 累積配列の最終値（= 全要素分）を使う
            cnt_bm1 = accum_cnt[M]
            sum_bm1 = accum_sum[M]
        else:
            cnt_bm1 = accum_cnt[idx]
            sum_bm1 = accum_sum[idx]
        bm1 = idx

    if cnt_bm1 == n:
        out.append("-1")
    else:
        ans = 1 + sum_bm1 + (n - cnt_bm1) * bm1
        out.append(str(ans))

print("\n".join(out))
