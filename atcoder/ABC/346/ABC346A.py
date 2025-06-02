n= int(input())

a = list(map(int,input().split()))

x = 0
while x < n-1:
    y = a[x]*a[x+1]
    print(y)
    x = x + 1