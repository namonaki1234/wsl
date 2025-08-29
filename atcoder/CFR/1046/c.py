import io
import sys

# 下記に標準入力を記載
_InPUT = """\
6
1
1
2
2 2
4
2 2 1 1
6
1 2 3 3 3 1
8
8 8 8 8 8 8 8 7
10
2 3 3 1 2 3 5 1 1 7

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

def solve():
    t = int(input())
    for _ in range(t):
        n = int(input())
        a = list(map(int, input().split()))

        # 値ごとの出現位置（1-indexed）
        pos = {}
        for idx, val in enumerate(a, start=1):
            pos.setdefault(val, []).append(idx)

        # 各終端 R に対して、そのRで終わるブロックの(始点L, 重みv)を格納
        ends = [[] for _ in range(n + 1)]
        for v, p in pos.items():
            m = len(p)
            if m >= v:
                # 連続する v 個で 1 ブロック（これが最短スパンで最適）
                for i in range(0, m - v + 1):
                    L = p[i]
                    R = p[i + v - 1]
                    ends[R].append((L, v))

        dp = [0] * (n + 1)
        for i in range(1, n + 1):
            # best: 位置 i まで見たときの「最適値（最大長）」を一時的に保存するための変数
            # cand(candidate):「i で終わるブロック [L, i] を使った場合に得られる候補値」
            best = dp[i - 1]  # i を使わない
            for L, w in ends[i]:
                cand = dp[L - 1] + w
                if cand > best:
                    best = cand
            dp[i] = best

        print(dp[n])

if __name__ == "__main__":
    solve()
