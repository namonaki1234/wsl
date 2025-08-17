import cv2
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
DBG  = HERE / "debug"
DBG.mkdir(exist_ok=True)

def show_or_save(name: str, img):
    """WSLでも動く：imshowの代わりに保存"""
    cv2.imwrite(str(DBG / f"{name}.png"), img)
    print(f"[saved] {DBG / f'{name}.png'}")

# 相対 import
from ._detect_corners import (
    convex_poly, select_corners, convex_poly_fitted,
    normalize_corners, fit_size,detect_corners
)
from ._trim_board import trim_board


def show_fitted(img, x):
    cntr = np.int32(x.reshape((4, 2)))
    blank = np.copy(img)
    cv2.drawContours(blank, [cntr], -1, (0,255,0), 2)
    return blank


def split_image(img, x, y, line):
    Path("./raw_img").mkdir(exist_ok=True)
    height, width, channels = img.shape
    h = height // y
    w = width // x
    line_h = round(h * line)
    line_w = round(w * line)
    counter = 0
    for split_y in range(1, y+1):
        for split_x in range(1, x+1):
            counter += 1
            clp = img[(h*(split_y-1))+line_h:(h*split_y)-line_h,
                      (w*(split_x-1))+line_w:(w*split_x)-line_w]
            cv2.imwrite(f"./raw_img/{counter}.png", clp)
    return


# def draw_ruled_line(img, show=True):
#     base_size = 32
#     w = base_size * 15
#     h = base_size * 15
#     img = img.copy()
#     for i in range(10):
#         x = int((w / 9) * i)
#         y = int((h / 9) * i)
#         cv2.line(img, (x, 0), (x, h), (0, 0, 255), 1)
#         cv2.line(img, (0, y), (w, y), (0, 0, 255), 1)
#     if show:
#         show_or_save("lined", img)
#     return img

def draw_ruled_line(img, show=True):
    out = img.copy()
    H, W = out.shape[:2]
    for i in range(10):
        x = int(W * i / 9)
        y = int(H * i / 9)
        cv2.line(out, (x, 0), (x, H), (0, 0, 255), 1)
        cv2.line(out, (0, y), (W, y), (0, 0, 255), 1)
    if show:
        show_or_save("lined", out)
    return out


# この先：数字を「画像認識」する最小ステップ
# 1) 9×9に分割して各マスから数字のROIを抜く
def extract_digit_roi(cell, min_ratio=0.01, max_ratio=0.5):
    # 二値化（照明ムラに強い）
    g = cv2.cvtColor(cell, cv2.COLOR_BGR2GRAY) if cell.ndim==3 else cell
    th = cv2.adaptiveThreshold(g, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                               cv2.THRESH_BINARY_INV, 11, 2)
    # 細線除去
    th = cv2.morphologyEx(th, cv2.MORPH_OPEN, np.ones((2,2), np.uint8), iterations=1)

    cnts, _ = cv2.findContours(th, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts: return None
    h, w = th.shape
    area = h*w
    cnt = max(cnts, key=cv2.contourArea)
    ratio = cv2.contourArea(cnt)/area
    if ratio < min_ratio or ratio > max_ratio:
        return None
    x,y,wc,hc = cv2.boundingRect(cnt)
    roi = th[y:y+hc, x:x+wc]
    # 正方形にパディングして 28x28
    s = max(wc, hc)
    canv = np.zeros((s,s), np.uint8)
    canv[(s-hc)//2:(s-hc)//2+hc, (s-wc)//2:(s-wc)//2+wc] = roi
    return cv2.resize(canv, (28,28), interpolation=cv2.INTER_AREA)

# 2) OCR：テンプレート照合 → Tesseract の順で（テンプレートを用意すると精度が高い。無い場合は Tesseract を使う）
def load_templates(dir_path="templates"):
    T = {}
    p = Path(dir_path)
    for d in range(1,10):
        f = p/f"{d}.png"
        if f.exists():
            img = cv2.imread(str(f), cv2.IMREAD_GRAYSCALE)
            _, img = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
            if np.mean(img) < 127: img = 255 - img
            T[d] = cv2.resize(img, (28,28))
    return T

def match_template(roi28, templates):
    if not templates: return None, 0.0
    r = roi28
    if np.mean(r) < 127: r = 255 - r
    best_d, best_s = None, -1.0
    for d, t in templates.items():
        s = float(cv2.matchTemplate(r, t, cv2.TM_CCOEFF_NORMED).max())
        if s>best_s: best_s, best_d = s, d
    return best_d, best_s

def ocr_tesseract(roi28):
    try:
        import pytesseract
    except Exception:
        return None, 0.0
    r = roi28
    if np.mean(r) < 127: r = 255 - r
    cfg = r'--psm 10 --oem 3 -c tessedit_char_whitelist=123456789'
    text = pytesseract.image_to_string(r, config=cfg).strip()
    return (int(text), 0.5) if len(text)==1 and text.isdigit() and text!='0' else (None, 0.0)

# 3) 盤面から 9×9 配列を作る
def read_grid(board_bgr, use_tesseract=True, templates_dir="templates",
              tm_thresh=0.58, ocr_ok_score=0.55):   # ← 追加: しきい値
    H, W = board_bgr.shape[:2]
    ch, cw = H//9, W//9
    templates = load_templates(templates_dir)
    grid   = np.zeros((9,9), dtype=int)
    scores = np.zeros((9,9), dtype=float)          # ← デバッグ用

    for r in range(9):
        for c in range(9):
            cell = board_bgr[r*ch:(r+1)*ch, c*cw:(c+1)*cw]
            pad = max(2, min(ch,cw)//15)
            cell = cell[pad:-pad, pad:-pad]

            roi = extract_digit_roi(cell)
            if roi is None:
                continue

            # ① テンプレ照合
            d, s = match_template(roi, templates)

            # ② スコアがしきい値未満なら Tesseract へフォールバック
            if (d is None or s < tm_thresh) and use_tesseract:
                d2, s2 = ocr_tesseract(roi)   # s2は0.5を仮に返すよう実装済み
                # OCRが読めた & ある程度自信があるときだけ採用
                if d2 is not None and max(s, s2) >= ocr_ok_score:
                    d, s = d2, max(s, s2)

            # ③ まだ曖昧なら未判定(0)のまま
            if d is not None and s >= tm_thresh:
                grid[r, c] = d
                scores[r, c] = s
            else:
                grid[r, c] = 0
                scores[r, c] = s

    # デバッグ: スコアマップも出しておくと便利
    from pathlib import Path
    here = Path(__file__).resolve().parent
    np.savetxt(here/"debug"/"score_map.csv", scores, fmt="%.3f", delimiter=",")
    return grid

# def read_grid(board_bgr, use_tesseract=True, templates_dir="templates"):
#     H, W = board_bgr.shape[:2]
#     ch, cw = H//9, W//9
#     templates = load_templates(templates_dir)
#     grid = np.zeros((9,9), dtype=int)

#     for r in range(9):
#         for c in range(9):
#             cell = board_bgr[r*ch:(r+1)*ch, c*cw:(c+1)*cw]
#             # 枠の影響を減らすため少し内側を使う
#             pad = max(2, min(ch,cw)//15)
#             cell = cell[pad:-pad, pad:-pad]
#             roi = extract_digit_roi(cell)
#             if roi is None: continue
#             d, s = match_template(roi, templates)
#             if (d is None or s<0.6) and use_tesseract:
#                 d2, s2 = ocr_tesseract(roi)
#                 if d2 is not None:
#                     d, s = d2, max(s, s2)
#             if d is not None:
#                 grid[r,c] = d
#     return grid



if __name__ == "__main__":
    raw_img = cv2.imread("/root/wsl/python/shogi-camera/shogicam/preprocess/sudoku.png")
    if raw_img is None:
        raise RuntimeError("画像が読み込めませんでした")

    # コーナー検出: (points, score) を返す想定
    det = detect_corners(raw_img)
    if isinstance(det, tuple) and len(det) == 2:
        pts, score = det
    else:
        # 実装によっては points だけ返す場合に備えて
        pts, score = det, None

    pts = np.asarray(pts, dtype=np.float32).reshape(-1, 2)
    if pts.shape != (4, 2):
        raise ValueError(f"detect_corners の戻り値が想定外です: shape={pts.shape}")

    # 角の順序を (tl, tr, br, bl) に正規化
    pts_norm = normalize_corners(pts)

    # デバッグ出力
    print("score:", score)
    print("pts_norm:\n", pts_norm)

    # 盤面の射影補正
    fit_img = trim_board(raw_img, pts_norm)   # ← 点だけ渡す（スコアは渡さない）
    show_or_save("fit_img", fit_img)

    # グリッド線描画
    lined = draw_ruled_line(fit_img)
    show_or_save("lined", lined)

    # 9×9に分割（少し内側を使いたいなら line を 0.02〜0.06 に）
    split_image(fit_img, 9, 9, 0.05)

    grid = read_grid(fit_img, use_tesseract=True, templates_dir="templates")
    np.savetxt(HERE/"debug"/"grid.csv", grid, fmt="%d", delimiter=",")
    print(grid)
    print("saved:", HERE/"debug"/"grid.csv")
