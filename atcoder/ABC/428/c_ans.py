import io
import sys

# 下記に標準入力を記載
_InPUT = """\
8
1 (
2
1 (
1 )
2
1 (
1 )
1 )
"""
sys.stdin = io.StringIO(_InPUT)

# ここからコードを記載
q = int(input())
a = [0]  # 現在の累積バランス
b = [0]  # ここまでの最小バランス
out = []

for _ in range(q):
    parts = input().split()
    if parts[0] == '1':
        c = parts[1]
        a.append(a[-1] + (1 if c == '(' else -1))
        b.append(min(b[-1], a[-1]))
    else:  # parts[0] == '2'
        a.pop()
        b.pop()
    out.append("Yes" if a[-1] == 0 and b[-1] == 0 else "No")

print("\n".join(out))
