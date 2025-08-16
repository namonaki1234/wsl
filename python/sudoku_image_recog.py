import tkinter as tk
from tkinter import messagebox
from tkinter import font as tkfont  # ★追加



# ---------- Sudoku logic ----------
def is_valid(grid, r, c, v):
    for i in range(9):
        if grid[r][i] == v: return False
        if grid[i][c] == v: return False
    br = (r // 3) * 3
    bc = (c // 3) * 3
    for i in range(br, br + 3):
        for j in range(bc, bc + 3):
            if grid[i][j] == v:
                return False
    return True

def find_empty(grid):
    for r in range(9):
        for c in range(9):
            if grid[r][c] == 0:
                return r, c
    return None

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    r, c = pos
    for v in range(1, 10):
        if is_valid(grid, r, c, v):
            grid[r][c] = v
            if solve(grid):
                return True
            grid[r][c] = 0
    return False

# ---------- GUI ----------
class SudokuGUI:
    def __init__(self, root):
        self.root = root
        root.title("Sudoku Solver (tkinter)")

    # ★共通フォントとズーム用の状態
        self.zoom_var = tk.IntVar(value=16)            # 初期フォントサイズ
        self.font = tkfont.Font(family="Segoe UI", size=self.zoom_var.get())


        self.entries = [[None]*9 for _ in range(9)]
        self._build_grid()
        self._build_buttons()
        self._build_zoom()   # ★ズームUI

          # ★ここにショートカットのバインドを追加
        root.bind("<Control-plus>", self._zoom_in)       # Ctrl + ＋
        root.bind("<Control-KP_Add>", self._zoom_in)     # Ctrl + テンキー＋
        root.bind("<Control-minus>", self._zoom_out)     # Ctrl + −
        root.bind("<Control-KP_Subtract>", self._zoom_out) # Ctrl + テンキー−
        root.bind("<Control-0>", self._zoom_reset)       # Ctrl + 0

    def _build_grid(self):
        frame = tk.Frame(self.root, padx=6, pady=6)
        frame.pack()
        self.grid_frame = frame  # ★参照保持（必要なら）

        FONT = ("Segoe UI", 16)
        for r in range(9):
            for c in range(9):
                block = (r // 3 + c // 3) % 2
                bg = "#f3f6ff" if block == 0 else "#eef3f0"
                e = tk.Entry(frame, width=2,
                             font=self.font,  # ★共通フォント適用
                             justify="center",
                             bg=bg, relief="solid", bd=1)
                e.grid(row=r, column=c, padx=(0 if c%3 else 2, 2), pady=(0 if r%3 else 2, 2))

                # only empty or single digit 0-9
                vcmd = (self.root.register(self._validate_digit), "%P")
                e.config(validate="key", validatecommand=vcmd)

                # focus: select all
                e.bind("<FocusIn>", lambda ev: ev.widget.select_range(0, tk.END))

                # arrow keys
                e.bind("<Up>",    lambda ev, rr=r, cc=c: self._move_focus(rr, cc, -1,  0))
                e.bind("<Down>",  lambda ev, rr=r, cc=c: self._move_focus(rr, cc,  1,  0))
                e.bind("<Left>",  lambda ev, rr=r, cc=c: self._move_focus(rr, cc,  0, -1))
                e.bind("<Right>", lambda ev, rr=r, cc=c: self._move_focus(rr, cc,  0,  1))

                # Enter to move down / Shift+Enter to move up (also support keypad Enter)
                e.bind("<Return>",       lambda ev, rr=r, cc=c: self._move_focus(rr, cc,  1, 0))
                e.bind("<KP_Enter>",     lambda ev, rr=r, cc=c: self._move_focus(rr, cc,  1, 0))
                e.bind("<Shift-Return>", lambda ev, rr=r, cc=c: self._move_focus(rr, cc, -1, 0))

                self.entries[r][c] = e

    def _build_zoom(self):
        # ★ズーム用のスライダー
        zf = tk.Frame(self.root, pady=4)
        zf.pack()
        tk.Label(zf, text="Zoom (cell size):").pack(side="left", padx=(0,6))
        scale = tk.Scale(zf, from_=12, to=40, orient="horizontal",
                        variable=self.zoom_var, command=self._on_zoom_change)
        scale.pack(side="left")

    def _on_zoom_change(self, _):
        size = int(self.zoom_var.get())
        self.font.configure(size=size)  # ★共通フォントを書き換えるだけで全Entryが追従
        # 必要ならウィンドウサイズも調整したい場合は geometry 再計算などを入れる

    def _zoom_in(self, event=None):
        self.zoom_var.set(min(self.zoom_var.get()+1, 40))
        self._on_zoom_change(None)

    def _zoom_out(self, event=None):
        self.zoom_var.set(max(self.zoom_var.get()-1, 12))
        self._on_zoom_change(None)

    def _zoom_reset(self, event=None):
        self.zoom_var.set(16)
        self._on_zoom_change(None)



    def _move_focus(self, r, c, dr, dc):
        nr = max(0, min(8, r + dr))
        nc = max(0, min(8, c + dc))
        target = self.entries[nr][nc]
        target.focus_set()
        target.icursor(tk.END)
        target.select_range(0, tk.END)
        return "break"  # 既定のキー動作を抑制

    def _build_buttons(self):
        btns = tk.Frame(self.root, pady=8)
        btns.pack()
        tk.Button(btns, text="Solve", width=10, command=self.on_solve).grid(row=0, column=0, padx=4)
        tk.Button(btns, text="Clear", width=10, command=self.on_clear).grid(row=0, column=1, padx=4)
        tk.Button(btns, text="Fill Sample", width=12, command=self.on_sample).grid(row=0, column=2, padx=4)

    # Only allow "", or one digit 0-9
    def _validate_digit(self, proposed: str) -> bool:
        if proposed == "":
            return True
        if len(proposed) > 1:
            return False
        return proposed in "0123456789"

    # Entry -> 2D grid
    def read_grid(self):
        grid = [[0]*9 for _ in range(9)]
        for r in range(9):
            for c in range(9):
                s = self.entries[r][c].get().strip()
                if s == "" or s == "0":
                    grid[r][c] = 0
                else:
                    try:
                        v = int(s)
                        if not (1 <= v <= 9):
                            raise ValueError
                        grid[r][c] = v
                    except ValueError:
                        raise ValueError(f"Please put a digit 1–9 at ({r+1},{c+1}).")
        return grid

    # 2D grid -> Entry
    def write_grid(self, grid):
        for r in range(9):
            for c in range(9):
                self.entries[r][c].delete(0, tk.END)
                v = grid[r][c]
                if v != 0:
                    self.entries[r][c].insert(0, str(v))

    # Check conflicts
    def check_initial_conflict(self, grid):
        for r in range(9):
            seen = set()
            for c in range(9):
                v = grid[r][c]
                if v == 0: continue
                if v in seen:
                    return f"Row {r+1} has duplicates."
                seen.add(v)
        for c in range(9):
            seen = set()
            for r in range(9):
                v = grid[r][c]
                if v == 0: continue
                if v in seen:
                    return f"Column {c+1} has duplicates."
                seen.add(v)
        for br in range(0, 9, 3):
            for bc in range(0, 9, 3):
                seen = set()
                for r in range(br, br+3):
                    for c in range(bc, bc+3):
                        v = grid[r][c]
                        if v == 0: continue
                        if v in seen:
                            return f"Block ({br+1}–{br+3} rows, {bc+1}–{bc+3} cols) has duplicates."
                        seen.add(v)
        return None

    def on_solve(self):
        try:
            grid = self.read_grid()
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        msg = self.check_initial_conflict(grid)
        if msg:
            messagebox.showerror("Initial Grid Error", msg)
            return

        g2 = [row[:] for row in grid]
        if solve(g2):
            self.write_grid(g2)
            messagebox.showinfo("Success", "Solved successfully.")
        else:
            messagebox.showwarning("Unsolved", "This puzzle could not be solved.")

    def on_clear(self):
        for r in range(9):
            for c in range(9):
                self.entries[r][c].delete(0, tk.END)

    def on_sample(self):
        sample = [
            [3,1,0,2,0,8,6,0,0],
            [0,2,0,0,0,6,3,1,0],
            [0,4,0,3,1,5,0,2,0],
            [0,9,0,0,0,0,0,7,6],
            [0,0,0,6,0,0,0,0,0],
            [0,0,1,9,0,0,0,3,2],
            [0,0,0,0,6,0,2,0,9],
            [7,0,2,0,3,9,4,0,1],
            [1,5,9,0,0,0,0,0,0],
        ]
        self.write_grid(sample)

if __name__ == "__main__":
    root = tk.Tk()
    app = SudokuGUI(root)
    root.resizable(True, False)
    root.mainloop()
