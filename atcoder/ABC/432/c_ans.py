import io
import sys

# 下記に標準入力を記載
_InPUT = """\
3 6 8
11 10 13
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
# 小さい飴をがxグラム、大きい飴がyグラム、ただしx < y
# n人の子供に飴を配るとき、各子供がもらう飴の重さの合計がすべて同じになるように配ることができるか？
import sys

def solve(n, x, y, a):
    # A を昇順に並べ替えても答えは変わらないので、念のためソートしておく
    a.sort()
    
    mina = a[0]  # A_1（最小値）

    total = 0
    vt = y - x   # (Y - X)

    for ai in a:
        # vs = (A_k - A_1) * Y
        vs = (ai - mina) * y

        # D_k = vs / vt が整数でなければ条件を満たせない
        if vs % vt != 0:
            return -1

        # minor = D_k = (A_k - A_1) * Y / (Y - X)
        minor = vs // vt

        # D_k > A_k になってしまうと x_k が負になってしまうので NG
        if minor > ai:
            return -1

        # x_k = A_k - D_k
        total += (ai - minor)

    return total


def main():
    input = sys.stdin.readline
    n, x, y = map(int, input().split())
    a = list(map(int, input().split()))
    ans = solve(n, x, y, a)
    print(ans)

if __name__ == "__main__":
    main()
