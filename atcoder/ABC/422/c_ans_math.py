import io
import sys

# 下記に標準入力を記載
_InPUT = """\
5
3 1 2
100 0 0
1000000 1000000 1000000
31 41 59
1000000000 10000 1
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

def main():
    T = int(input())
    for _ in range(T):
        n_A, n_B, n_C = map(int, input().split())
        ans = min(n_A, n_C, (n_A + n_B + n_C) // 3)
        print(ans)

if __name__ == "__main__":
    main()
