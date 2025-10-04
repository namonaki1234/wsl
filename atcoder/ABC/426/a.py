import io
import sys


# 下記に標準入力を記載
_InPUT = """\
Ocelot Ocelot

"""
sys.stdin = io.StringIO(_InPUT)
# ここからコードを記載

x,y = input().split()

# os = {0:"Serval",1:"Ocelot",2:"Lynx"}

if y == "Serval":
    if x == "Serval" or x == "Lynx":
        print("Yes")
        exit()
elif y == "Ocelot":
    if x == "Ocelot" or x == "Serval" or x == "Lynx":
        print("Yes")
        exit()
elif y == "Lynx":
    if x == "Lynx":
        print("Yes")
        exit()

print("No")
