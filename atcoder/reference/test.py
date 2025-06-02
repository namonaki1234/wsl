import io
import sys

# 下記に標準入力を記載
_INPUT = """\
ABC000



"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載
S = input()

if S[3] == "0":
    if S[4] == "0":
        for i in range(0,10):
            if i != 0:
                if S[5] == str(i):
                    print("Yes")
                    break
            if S[5] == "0":
                print("No")
                break
    else:
        for i in range(10,100):
            if S[4:6] == str(i):
                print("Yes")
                break
                
elif S[3] != "0":
    if int(S[3:6]) < 350 and int(S[3:6]) > 1:
        for i in range(100,350):
            if i != 316:
                if S[3:6] == str(i):
                    print("Yes")
                    break
            if int(S[3:6]) == 316:
                print("No")
                break
                
    else:
        print("No")   
        

    
    
