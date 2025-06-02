def check_ad(S, N, titles):
    results = []
    
    for title in titles:
        if S in title:
            results.append("Yes")
        else:
            results.append("No")
    
    return results

S = input() 
N = int(input())  
titles = [input().strip() for _ in range(N)] 

results = check_ad(S, N, titles)

for result in results:
    print(result)
