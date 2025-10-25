import sys
import io

# 下記に標準入力を記載
_InPUT = """\
5
3 2 5 2 2
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
N = int(input().strip())
A = list(map(int, input().split()))

# 値域が 1..N なので、長さ N+1 のカウント配列（0番は未使用）
B = [0] * (N + 1)
for x in A:
    B[x] += 1

ans = 0
for cnt in B[1:]:  # 1..N だけ見る
    if cnt >= 2:
        ans += cnt * (cnt - 1) // 2 * (N - cnt)
print(ans)

