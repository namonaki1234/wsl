import io
import sys

# 下記に標準入力を記載
_InPUT = """\
17
AAABABABBBABABBABABABABBAAABABABBA
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

N = int(input().strip())
S = input().strip()
# Aの位置（1-indexed）を収集
X = [i for i, ch in enumerate(S, start=1) if ch == 'A']  # 長さ N のはず

# 目標1: ABABAB... → A は奇数位置 (1,3,5,...)
cost_odd = 0
for k, x in enumerate(X):  # k = 0..N-1
    target = 2*k + 1
    cost_odd += abs(x - target)

# 目標2: BABABA... → A は偶数位置 (2,4,6,...)
cost_even = 0
for k, x in enumerate(X):
    target = 2*(k+1)
    cost_even += abs(x - target)

print(min(cost_odd, cost_even))
