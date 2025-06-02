import io
import sys

# 下記に標準入力を記載
_INPUT = """\
32
1 2 4 5 8 10 12 16 19 25 33 40 50 64 87 101 149 175 202 211 278 314 355 405 412 420 442 481 512 582 600 641


"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

N = int(input())
A = list(map(int, input().split()))

kind_count = 0
if A[N-1]/2<A[0]:
    print(0)

else:
    for j in range(N):
        for i in range(N):
            if A[N-1-j]/2>=A[i]:
                kind_count += 1
                continue
            
            else:
                
                break
            
                
            
        
    print(kind_count)




