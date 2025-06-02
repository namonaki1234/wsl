import io
import sys

# 下記に標準入力を記載
_INPUT = """\
WWWWA

"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

S = input()
S = list(S)

'''
range(start, stop, step)の形式
start: ループの開始値。
stop: ループの終了値（この値は含まれません）。ループは-1の手前、つまり0まで繰り返されます。
step: ループの増分（または減分）。-1ならループがstartからstopに向かって1ずつ減少することを意味します。
'''

for i in range(len(S) - 2, -1, -1): 
    
    if S[i] == 'W' and S[i+1] == 'A':
             S[i] = 'A'
             S[i+1] = 'C' 
        
 

print(''.join(S))
    
  










