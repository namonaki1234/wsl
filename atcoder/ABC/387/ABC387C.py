import io
import sys

# 下記に標準入力を記載
_INPUT = """\
97 210

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載


def count_snake_numbers(n):
    s = str(n)
    length = len(s)
    total = 0

    # 桁数が少ない場合のヘビ数をカウント
    for i in range(1, length):
        for j in range(1, 10):
            total += j ** (i - 1)

    # 最上位桁がs[0]の数をカウント
    first_digit = int(s[0])
    for i in range(1, first_digit):
        total += i ** (length - 1)

    # 最上位桁がs[0]で、残りの桁を考慮
    for i in range(1, length):
        if int(s[i]) < first_digit:
            total += first_digit ** (length - i - 1)
        elif int(s[i]) == first_digit:
            if i == length - 1:
                total += 1  # 最後の桁が同じ場合は1を加える
            break

    return total

L, R = map(int, input().split())
result = count_snake_numbers(R) - count_snake_numbers(L - 1)
print(result)
