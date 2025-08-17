import cv2 as cv
import numpy as np
from pathlib import Path

# --- 1) 盤面抽出（最大四角形を射影変換） ---
def find_board_and_warp(img_bgr, out_size=900):
    img = cv.cvtColor(img_bgr, cv.COLOR_BGR2GRAY)
    blur = cv.GaussianBlur(img, (7,7), 0)
    thr = cv.adaptiveThreshold(blur, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C,
                               cv.THRESH_BINARY_INV, 11, 2)
    # 盤面の太線を強める
    kernel = cv.getStructuringElement(cv.MORPH_RECT, (5,5))
    thr = cv.dilate(thr, kernel, iterations=1)

    # 輪郭から最大四角形
    cnts, _ = cv.findContours(thr, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    if not cnts:
        raise RuntimeError("盤面らしき輪郭が見つかりませんでした")
    cnt = max(cnts, key=cv.contourArea)
    peri = cv.arcLength(cnt, True)
    approx = cv.approxPolyDP(cnt, 0.02*peri, True)
    if len(approx) != 4:
        raise RuntimeError("四角形として近似できませんでした")

    # 頂点の並びを (tl, tr, br, bl) に並べ替え
    pts = approx.reshape(4,2).astype(np.float32)
    s = pts.sum(axis=1); d = np.diff(pts, axis=1).ravel()
    tl = pts[np.argmin(s)]; br = pts[np.argmax(s)]
    tr = pts[np.argmin(d)]; bl = pts[np.argmax(d)]
    src = np.array([tl,tr,br,bl], dtype=np.float32)
    dst = np.array([[0,0],[out_size-1,0],[out_size-1,out_size-1],[0,out_size-1]], dtype=np.float32)

    M = cv.getPerspectiveTransform(src, dst)
    warped = cv.warpPerspective(img_bgr, M, (out_size, out_size))
    return warped

# --- 2) セル切り出し & 数字抽出用の二値画像 ---
def preprocess_board_gray(board_bgr):
    g = cv.cvtColor(board_bgr, cv.COLOR_BGR2GRAY)
    g = cv.GaussianBlur(g, (3,3), 0)
    th = cv.adaptiveThreshold(g, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C,
                              cv.THRESH_BINARY_INV, 11, 2)
    # 細線を少し消して数字を残しやすく
    kernel = cv.getStructuringElement(cv.MORPH_RECT, (2,2))
    th = cv.morphologyEx(th, cv.MORPH_OPEN, kernel, iterations=1)
    return th

def extract_cell(th_img, r, c, N=9):
    h, w = th_img.shape
    ch, cw = h//N, w//N
    y0, y1 = r*ch, (r+1)*ch
    x0, x1 = c*cw, (c+1)*cw
    cell = th_img[y0:y1, x0:x1]

    # 枠の影響を減らすため、少し内側を使う
    pad = max(2, min(ch,cw)//15)
    cell = cell[pad:-pad, pad:-pad]
    return cell

def extract_digit_roi(cell_bin, min_area_ratio=0.01, max_area_ratio=0.5):
    if cell_bin.size == 0:
        return None
    cnts, _ = cv.findContours(cell_bin, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    if not cnts: return None
    h, w = cell_bin.shape
    area = h*w
    # 最大連結成分を採用
    cnt = max(cnts, key=cv.contourArea)
    a = cv.contourArea(cnt) / area
    if a < min_area_ratio or a > max_area_ratio:
        return None
    x,y,wc,hc = cv.boundingRect(cnt)
    roi = cell_bin[y:y+hc, x:x+wc]
    # 余白を付けつつ正方形にパディング
    size = max(wc, hc)
    canvas = np.zeros((size, size), dtype=np.uint8)
    yoff = (size - hc)//2; xoff = (size - wc)//2
    canvas[yoff:yoff+hc, xoff:xoff+wc] = roi
    # 28x28へ正規化
    roi28 = cv.resize(canvas, (28,28), interpolation=cv.INTER_AREA)
    return roi28

# --- 3-A) テンプレート照合（同じアプリ/フォントなら強い） ---
def load_templates(dir_path="templates"):
    T = {}
    for d in range(1,10):
        p = Path(dir_path)/f"{d}.png"
        if p.exists():
            img = cv.imread(str(p), cv.IMREAD_GRAYSCALE)
            _, img = cv.threshold(img, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
            T[d] = cv.resize(img, (28,28))
    return T

def match_with_templates(roi28, templates):
    # 相関で最大のものを採用
    best_d, best_score = None, -1.0
    if not templates: return None, 0.0
    for d, tmpl in templates.items():
        res = cv.matchTemplate(roi28, tmpl, cv.TM_CCOEFF_NORMED)
        score = float(res.max())
        if score > best_score:
            best_score, best_d = score, d
    return best_d, best_score

# --- 3-B) pytesseract フォールバック（任意） ---
def ocr_tesseract_digit(roi28):
    try:
        import pytesseract
    except ImportError:
        return None, 0.0
    cfg = r'--psm 10 --oem 3 -c tessedit_char_whitelist=123456789'
    # 反転して黒字→白字に合わせることも（環境次第）
    inv = cv.bitwise_not(roi28)
    txt = pytesseract.image_to_string(inv, config=cfg)
    txt = txt.strip()
    if len(txt)==1 and txt.isdigit() and txt!='0':
        return int(txt), 0.5  # Tesseractは信頼値APIが重いので暫定
    return None, 0.0

# --- 4) まとめ：画像 → 9x9 配列 ---
def read_sudoku(image_path, use_tesseract=True, template_dir="templates", debug=False):
    bgr = cv.imread(image_path)
    if bgr is None:
        raise FileNotFoundError(image_path)

    board = find_board_and_warp(bgr)
    th = preprocess_board_gray(board)
    templates = load_templates(template_dir)

    grid = np.zeros((9,9), dtype=int)
    scores = np.zeros((9,9), dtype=float)

    for r in range(9):
        for c in range(9):
            cell = extract_cell(th, r, c)
            roi = extract_digit_roi(cell)
            if roi is None:
                continue
            # まずテンプレート照合
            d, s = match_with_templates(roi, templates)
            # スコアが低い場合はTesseractへ
            if (d is None or s < 0.55) and use_tesseract:
                d2, s2 = ocr_tesseract_digit(roi)
                if d2 is not None:
                    d, s = d2, max(s, s2)
            if d is not None:
                grid[r,c] = d
                scores[r,c] = s

    if debug:
        print(grid)
        print("avg score:", scores[grid>0].mean() if (grid>0).any() else 0.0)
    return grid, scores
