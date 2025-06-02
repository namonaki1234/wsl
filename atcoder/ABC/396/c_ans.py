N, M = map(int, input().split())
A = list(map(int, input().split()))
B = list(map(int, input().split()))
A.sort(reverse=True)
B.sort(reverse=True)
S, T, maxT = [0] * (N + 1), [0] * (M + 1), [0] * (M + 1)
for i in range(N):
    S[i + 1] = S[i] + A[i]
for i in range(M):
    T[i + 1] = T[i] + B[i]
    maxT[i + 1] = max(maxT[i], T[i + 1])

ans = 0
for i in range(N + 1):
    ans = max(ans, S[i] + maxT[min(i, M)])
print(ans)
