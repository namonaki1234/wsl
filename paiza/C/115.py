
n,m=map(int,input().split())

carmeter = [int(input()) for _ in range(n-1)]

trafficjam =0


for i in range(n-1):
    if carmeter[i]<= m:
        trafficjam+=carmeter[i]


print(trafficjam)

