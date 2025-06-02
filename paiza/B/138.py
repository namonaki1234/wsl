# 下記に標準入力を記載
_INPUT = """\
5 5
####.
#.###
###.#
#.###
###.#

"""

def count_donuts(H, W, grid):
    
    donut= [
        "###",
        "#.#",
        "###"
    ]
    
    count = 0
    
    
    for i in range(1, H-1):
        for j in range(1, W-1):
            # 3x3の範囲をチェック
            if (grid[i-1][j-1] == "#" and grid[i-1][j] == "#" and grid[i-1][j+1] == "#" and
                grid[i][j-1] == "#" and grid[i][j] == "." and grid[i][j+1] == "#" and
                grid[i+1][j-1] == "#" and grid[i+1][j] == "#" and grid[i+1][j+1] == "#"):
                count += 1
                
    return count


H, W = map(int, input().split())
grid = [input().strip() for _ in range(H)]

print(count_donuts(H, W, grid))
