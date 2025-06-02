import io
import sys

# 下記に標準入力を記載
_INPUT = """\
12 35 7 10
###################################
#.................................#
#..........@......................#
#......@................@.........#
#.............##............@.....#
#...##........##....##............#
#...##........##....##.......##...#
#....##......##......##....##.....#
#....##......##......##..##.......#
#.....#######.........###.........#
#.................................#
###################################
LRURRRUUDDULUDUUDLRLRDRRLULRRUDLDRU



"""
sys.stdin = io.StringIO(_INPUT)
# ここからコードを記載

H,W,X,Y = map(int,input().split())
S = [input().strip() for _ in range(H)]
T = input()
X-=1
Y-=1

found_index = 0

for i in range(len(T)):
    # print(S[X][Y])
    
    if T[i] == "U":
        X -= 1
        if S[X][Y] == "@" or S[X][Y] == ".":
            if S[X][Y] == "@":
                found_index+=1  
                S_list = list(S[X])
                S_list[Y] = "."
                S[X] = "".join(S_list)
        else:
            X += 1
    
    elif T[i] == "D":
        X += 1
        if S[X][Y] == "@" or S[X][Y] == ".":
            if S[X][Y] == "@":
                found_index+=1  
                S_list = list(S[X])
                S_list[Y] = "."
                S[X] = "".join(S_list)
        else:    
            X -= 1
        
    elif T[i] == "L":
        Y -= 1
        if S[X][Y] == "@" or S[X][Y] == ".":
            if S[X][Y] == "@":
                found_index+=1  
                S_list = list(S[X])
                S_list[Y] = "."
                S[X] = "".join(S_list)
        else:    
            Y += 1  
        
    elif T[i] == "R":
        Y += 1
        if S[X][Y] == "@" or S[X][Y] == ".":
            if S[X][Y] == "@":
                found_index+=1  
                S_list = list(S[X])
                S_list[Y] = "."
                S[X] = "".join(S_list)
        else:    
            Y -= 1
    
X,Y = X+1,Y+1
print(X,Y,found_index)