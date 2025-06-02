import io
import sys
from collections import defaultdict,deque,Counter
import math
# from more_itertools import chunked

# 下記に標準入力を記載
_INPUT = """\
5 4
2 1 2
3 3 4 5
3 1 2 5
1 3
1 3 2 5 4


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
N, M = map(int, input().split())
KA = [list(map(int, input().split())) for _ in range(M)]
B =list(map(int, input().split())) 
KA = deque(KA)
count = 0
KA_remove = []
# for _ in range(N):
for b in B:
    for i in range(M):
        # for j in sorted(range(KA[i][0]), reverse=True):
        for j in range(KA[i][0], 0, -1):
            KA_remove.clear()
            
            if len(KA[i]) > j and KA[i][j] == b:
                KA_remove.append(j)
            for u in range(len(KA_remove)):
                # for j in range(len(KA)):
                #     if KA_remove[i] in KA[j]:
                #         KA[j].remove(KA_remove[i])
                KA[i].pop(KA_remove[u])

    for i in range(M):
        if len(KA[i]) == 1:
            KA[i].pop(0)
            count += 1
    print(count)




