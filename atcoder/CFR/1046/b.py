import io
import sys

# 下記に標準入力を記載
_InPUT = """\
6
2 1
00
4 3
0010
5 2
11011
7 5
1111110
8 4
00101011
10 2
1000000010
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

def solve():
    t = int(input())
    out = []
    for _ in range(t):
        n, k = map(int, input().split())
        s = input().strip()

        # 1) 「長さkで全て1」の区間があるかをスライドウィンドウで判定
        ones = sum(1 for i in range(min(k, n)) if s[i] == '1')
        possible = True
        if ones == k:
            possible = False
        else:
            for i in range(k, n):
                if s[i] == '1': ones += 1
                if s[i - k] == '1': ones -= 1
                if ones == k:
                    possible = False
                    break

        if not possible:
            out.append("NO")
            continue

        # 2) 構成：0 に大きい値、1 に小さい値
        p = [0] * n
        cur = n
        # まず 0 の場所に大きい値
        for i, ch in enumerate(s):
            if ch == '0':
                p[i] = cur
                cur -= 1
        # 次に 1 の場所に残りの値
        for i, ch in enumerate(s):
            if ch == '1':
                p[i] = cur
                cur -= 1

        out.append("YES")
        out.append(" ".join(map(str, p)))
    print("\n".join(out))

if __name__ == "__main__":
    solve()
