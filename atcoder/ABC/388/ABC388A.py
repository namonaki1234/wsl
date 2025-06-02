import io
import sys

# 下記に標準入力を記載
_INPUT = """\
Tohoku


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input()
S_UPC = S[0] + 'UPC'

print(S_UPC)