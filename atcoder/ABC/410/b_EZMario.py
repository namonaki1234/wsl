import sys
sys.setrecursionlimit(5*10**7)
import math
import bisect
import heapq
from collections import deque, defaultdict
import random
import itertools
#from decimal import Decimal
#from copy import copy, deepcopy
#from sortedcontainers import SortedList

# -------------------------------------------------
pin = sys.stdin.readline
def ST(): return pin().rstrip()
def IN(): return int(pin())
def IM(): return map(int, pin().split())
def IL(): return list(IM())
def SR(n:int)->list: return [pin().rstrip() for _ in range(n)]
def IMatrix(n:int)->list: return [IL() for _ in range(n)]

INF = 2*10**18+1
mod = 998244353
# -------------------------------------------------
#import pypyjit
#pypyjit.set_param("max_unroll_recursion=-1")

N = IN()
A = IL()
A.sort()

def keyfunc(n:int)->bool:
    global A, N
    if N-bisect.bisect_left(A, n) < n:
        return True
    return False

def binary_search(key = keyfunc, MX = 2**62)->int:
    ng = -1
    ok = MX
    while abs(ok - ng) > 1:
        mid = (ok+ng)//2
        if key(mid):
            ok = mid
        else:
            ng = mid
    return ok

ans = binary_search()-1
print(ans)

