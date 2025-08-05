import io
import sys
from collections import defaultdict,deque,Counter
import numpy as np
from itertools import combinations,permutations
import math

# 下記に標準入力を記載
_InPUT = """\
3 2 6
abc
xxx
abc
"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

n, k, x = map(int, input().split())
S = [input().strip() for _ in range(n)]

ans=[]
# 要はdfsを使って、k 個の文字列を全探索かつ配列で取得している。
def dfs(crr, count):
  # count 個の文字列を結合して crr になった状態
  if count==k:
    ans.append(crr)
    return
  for s in S:
    dfs(crr+s, count+1)
dfs("", 0)
ans.sort()
print(ans[x-1])


'''
dfs("", 0)
├─ dfs("abc", 1)
│  ├─ dfs("abcabc", 2) ← ansに追加
│  ├─ dfs("abcxxx", 2) ← ansに追加
│  └─ dfs("abcabc", 2) ← ansに追加
├─ dfs("xxx", 1)
│  ├─ dfs("xxxabc", 2) ← ansに追加
│  ├─ dfs("xxxxx", 2)  ← ansに追加
│  └─ dfs("xxxabc", 2) ← ansに追加
└─ dfs("abc", 1)
   ├─ dfs("abcabc", 2) ← ansに追加
   ├─ dfs("abcxxx", 2) ← ansに追加
   └─ dfs("abcabc", 2) ← ansに追加

'''