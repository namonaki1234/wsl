
import io
import sys

# 下記に標準入力を記載
_INPUT = """\
1
leap
read


"""
sys.stdin = io.StringIO(_INPUT)

# 入力の読み込み
K = int(input().strip())
S = input().strip()
T = input().strip()

"""
# S と T が等しい場合は "Yes" を出力
if S == T:
    result = "Yes"
else:
    len_s = len(S)
    len_t = len(T)

    if len_s == len_t:
        # 文字列の長さが同じ場合、異なる文字が1つだけか確認
        diff_count = sum(1 for a, b in zip(S, T) if a != b)
        result = "Yes" if diff_count == 1 else "No"
    elif len_s == len_t + 1:
        # S の方が1文字長い場合、1文字削除して一致するか確認
        result = "No"  # 初期値を "No" に設定
        for i in range(len_s):
            if S[:i] + S[i+1:] == T:
                result = "Yes"
                break  # 一致が見つかったらループを終了
    elif len_s + 1 == len_t:
        # S の方が1文字短い場合、1文字挿入して一致するか確認
        result = "No"  # 初期値を "No" に設定
        for i in range(len_t):
            if T[:i] + T[i+1:] == S:
                result = "Yes"
                break  # 一致が見つかったらループを終了
    else:
        result = "No"

# 結果の出力
print(result)
"""



# S と T が等しい場合は "Yes" を出力
if S == T:
    print("Yes")
else:
    len_s = len(S)
    len_t = len(T)

    if len_s == len_t:
        # 変更の場合
        diff_count = sum(1 for a, b in zip(S, T) if a != b)
        if diff_count <= 1:
            print("Yes")
        else:
            print("No")
    elif len_s + 1 == len_t:
        # 挿入の場合
        pc = 0  # 先頭一致のカウント
        sc = 0  # 末尾一致のカウント

        while pc < len_s and S[pc] == T[pc]:
            pc += 1
        while sc < len_s and S[len_s - 1 - sc] == T[len_t - 1 - sc]:
            sc += 1

        if pc + sc >= len_s:
            print("Yes")
        else:
            print("No")
    elif len_s - 1 == len_t:
        # 削除の場合
        # S と T を入れ替える
        S, T = T, S
        len_s, len_t = len(T), len(S)

        pc = 0  # 先頭一致のカウント
        sc = 0  # 末尾一致のカウント

        while pc < len_s and S[pc] == T[pc]:
            pc += 1
        while sc < len_s and S[len_s - 1 - sc] == T[len_t - 1 - sc]:
            sc += 1

        if pc + sc >= len_s:
            print("Yes")
        else:
            print("No")
    else:
        print("No")

