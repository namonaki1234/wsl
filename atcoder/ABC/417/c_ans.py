import io
import sys
from collections import Counter

# 下記に標準入力を記載
_InPUT = """\
9
3 1 4 1 5 9 2 6 5

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載
N=int(input())
A=list(map(int,input().split()))

# Ai+i の度数分布の作成
Ai_plus_i = Counter(A[i]+i for i in range(N))
# j-Ajの登場回数の和を計算
print(sum(Ai_plus_i[j-A[j]] for j in range(N)))
