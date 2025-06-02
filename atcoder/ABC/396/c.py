import io
import sys

# 下記に標準入力を記載
_INPUT = """\
3 5
-36 -33 -31
12 12 28 24 27

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N,M = map(int,input().split())
B = list(map(int,input().split()))
if N > M:
    W = list(map(int,input().split()))+ [0] * (N-M)
else:
    W = list(map(int,input().split()))

B.sort(reverse=True)
W.sort(reverse=True)

ans = 0
for i in zip(B,W):
    print(i) 
    
    if 0 < i[0]:
        ans += i[0]
        if 0 < i[1]:
            ans += i[1]
        else:
            if i[0]+i[1] < 0:
                if i[0] < i[1] and 0 < i[0]:
                    # このパターンはありえない
                    continue
                # (3,-4)の3はすでにたしてある
                
                print(ans)
                exit()
        
            continue
    else:
        
        if i[0]+i[1] > 0:
            # (-1,5)        
            if i[0] < i[1] and 0 > i[0]:
                ans += i[0]+i[1]
            print(ans)
            exit()
        else:
            # (-5,-1)
            if i[0] < i[1] and 0 > i[0]:
                continue
            # (-3,1)
            elif i[0] > i[1] and 0 > i[0]:
                continue
            print(ans)
            exit()
        
print(ans)













