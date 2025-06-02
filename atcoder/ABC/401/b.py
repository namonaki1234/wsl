import io
import sys
from collections import defaultdict,deque,Counter

# 下記に標準入力を記載
_INPUT = """\
20
private
login
private
logout
public
logout
logout
logout
logout
private
login
login
private
login
private
login
public
private
logout
private


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
S = [input().rstrip() for _ in range(N)]

log = False

err_count = 0
for i in range(N):
    if S[i] == "login":
        log = True
    elif S[i] == "logout":
        log = False
    
    if log == False:
        if S[i] == "private":
            err_count += 1

print(err_count)

    


    











