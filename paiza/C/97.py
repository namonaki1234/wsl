N,X,Y=map(int,input().split())

for i in range(1,N+1):
    if i%X==0 and i%Y==0:
        print("AB")

    elif i%Y==0 :
        print("B")
    
    elif i%X==0:
        print("A")
    
    else:
        print("N")