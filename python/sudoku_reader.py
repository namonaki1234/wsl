
import cv2 as cv
import numpy as np
from pathlib import Path

# ============================
# Sudoku Reader (OpenCV)
# ============================
# 1) 盤面の射影補正
# 2) 9x9分割 & セル内の数字候補を抽出（黒文字・青文字両対応）
# 3) 認識：テンプレート照合 → Tesseract（任意）
# 4) 結果の9x9グリッドとデバッグ画像を出力
#
# 使い方:
#   python sudoku_reader.py /path/to/sudoku.jpg --out /path/to/outdir --tess
#   python sudoku_reader.py /root/wsl/python/sudoku.jpg --out /root/wsl/python/out --tess
#
# テンプレート:
#   templates/1.png ... templates/9.png を用意するとテンプレ照合が有効になります。
#   無い場合は Tesseract (--tess) を使うか、ROI画像を元に作成してください。

def imwrite(p, img):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    cv.imwrite(str(p), img)

def find_board_and_warp(img_bgr, out_size=900):
    gray = cv.cvtColor(img_bgr, cv.COLOR_BGR2GRAY)
    blur = cv.GaussianBlur(gray, (7,7), 0)
    thr  = cv.adaptiveThreshold(blur, 255, cv.ADAPTIVE_THRESH_GAUSSIAN_C,
                                cv.THRESH_BINARY_INV, 11, 2)
    kernel = cv.getStructuringElement(cv.MORPH_RECT, (5,5))
    thr = cv.dilate(thr, kernel, iterations=1)

    cnts, _ = cv.findContours(thr, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    if not cnts:
        raise RuntimeError("盤面らしき輪郭が見つからない")
    # 面積上位から四角近似を試す
    for cnt in sorted(cnts, key=cv.contourArea, reverse=True)[:20]:
        peri = cv.arcLength(cnt, True)
        approx = cv.approxPolyDP(cnt, 0.02*peri, True)
        if len(approx)==4:
            pts = approx.reshape(4,2).astype(np.float32)
            s = pts.sum(axis=1); d = np.diff(pts, axis=1).ravel()
            tl = pts[np.argmin(s)]; br = pts[np.argmax(s)]
            tr = pts[np.argmin(d)]; bl = pts[np.argmax(d)]
            src = np.array([tl,tr,br,bl], np.float32)
            dst = np.array([[0,0],[out_size-1,0],[out_size-1,out_size-1],[0,out_size-1]], np.float32)
            M = cv.getPerspectiveTransform(src, dst)
            warped = cv.warpPerspective(img_bgr, M, (out_size, out_size))
            return warped
    raise RuntimeError("四角形として近似できない")

def make_digit_mask(board_bgr):
    """黒文字/青文字を同時に拾う2値画像を作る"""
    hsv = cv.cvtColor(board_bgr, cv.COLOR_BGR2HSV)
    # 黒(濃い灰含む)：V低め & S低〜中
    lower_black = np.array([0, 0, 0])
    upper_black = np.array([180, 120, 120])
    mask_black = cv.inRange(hsv, lower_black, upper_black)
    # 青：H ~ [90..140] 近辺（機種によるので広めに）
    lower_blue1 = np.array([90,  50, 50])
    upper_blue1 = np.array([140, 255,255])
    mask_blue = cv.inRange(hsv, lower_blue1, upper_blue1)
    mask = cv.bitwise_or(mask_black, mask_blue)
    # 細線除去・ノイズ除去
    kernel = cv.getStructuringElement(cv.MORPH_RECT, (3,3))
    mask = cv.morphologyEx(mask, cv.MORPH_OPEN, kernel, iterations=1)
    return mask

def extract_cell(img, r, c, N=9):
    h, w = img.shape[:2]
    ch, cw = h//N, w//N
    y0, y1 = r*ch, (r+1)*ch
    x0, x1 = c*cw, (c+1)*cw
    cell = img[y0:y1, x0:x1]
    # 枠の影響を避けるため少し内側を使う
    pad = max(2, min(ch,cw)//15)
    return cell[pad:-pad, pad:-pad]

def extract_digit_roi(cell_mask, min_area_ratio=0.01, max_area_ratio=0.6):
    if cell_mask.size == 0: return None
    cnts, _ = cv.findContours(cell_mask, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    if not cnts: return None
    h,w = cell_mask.shape[:2]
    area = h*w
    cnt = max(cnts, key=cv.contourArea)
    a = cv.contourArea(cnt)/area
    if a < min_area_ratio or a > max_area_ratio:
        return None
    x,y,wc,hc = cv.boundingRect(cnt)
    roi = cell_mask[y:y+hc, x:x+wc]
    # 正方形にパディング → 28x28
    size = max(wc, hc)
    canvas = np.zeros((size,size), dtype=np.uint8)
    yoff = (size-hc)//2; xoff = (size-wc)//2
    canvas[yoff:yoff+hc, xoff:xoff+wc] = roi
    roi28 = cv.resize(canvas, (28,28), interpolation=cv.INTER_AREA)
    return roi28

def load_templates(template_dir="templates"):
    T = {}
    p = Path(template_dir)
    for d in range(1,10):
        f = p/f"{d}.png"
        if f.exists():
            img = cv.imread(str(f), cv.IMREAD_GRAYSCALE)
            # 数字(白) on 黒背景に統一
            _, img = cv.threshold(img, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
            # 反転して数字を白に
            if np.mean(img) < 127:
                img = 255 - img
            img = cv.resize(img, (28,28), interpolation=cv.INTER_AREA)
            T[d] = img
    return T

def match_with_templates(roi28, templates):
    if not templates: return None, 0.0
    # ROIを白字(255)に揃える
    roi = roi28.copy()
    if np.mean(roi) < 127:  # もし黒字なら反転
        roi = 255 - roi
    best_d, best_score = None, -1.0
    for d, tmpl in templates.items():
        res = cv.matchTemplate(roi, tmpl, cv.TM_CCOEFF_NORMED)
        score = float(res.max())
        if score > best_score:
            best_score, best_d = score, d
    return best_d, best_score

def ocr_tesseract(roi28):
    try:
        import pytesseract
    except Exception:
        return None, 0.0
    # 文字は白で与える方が安定
    roi = roi28.copy()
    if np.mean(roi) < 127:
        roi = 255 - roi
    cfg = r'--psm 10 --oem 3 -c tessedit_char_whitelist=123456789'
    txt = pytesseract.image_to_string(roi, config=cfg)
    txt = txt.strip()
    if len(txt)==1 and txt.isdigit() and txt!='0':
        return int(txt), 0.5
    return None, 0.0

def read_sudoku(image_path, out_dir="read_out", template_dir="templates", use_tesseract=True, debug=True):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    bgr = cv.imread(str(image_path))
    if bgr is None: raise FileNotFoundError(image_path)

    board = find_board_and_warp(bgr, out_size=900)
    imwrite(out_dir/"01_warped.png", board)

    mask = make_digit_mask(board)              # 黒+青の2値化
    imwrite(out_dir/"02_mask.png", mask)

    templates = load_templates(template_dir)

    grid = np.zeros((9,9), dtype=int)
    score_map = np.zeros((9,9), dtype=float)

    H, W = mask.shape[:2]
    ch, cw = H//9, W//9
    contact = np.ones((H, W), dtype=np.uint8)*255

    for r in range(9):
        for c in range(9):
            cell = extract_cell(mask, r, c)
            imwrite(out_dir/f"cells/cell_{r}_{c}.png", cell)
            roi28 = extract_digit_roi(cell)
            if roi28 is None:
                continue
            imwrite(out_dir/f"rois/roi_{r}_{c}.png", roi28)

            # 1) テンプレート
            d, s = match_with_templates(roi28, templates)
            # 2) Tesseract フォールバック
            if (d is None or s < 0.6) and use_tesseract:
                d2, s2 = ocr_tesseract(roi28)
                if d2 is not None:
                    d, s = d2, max(s, s2)

            if d is not None:
                grid[r,c] = d
                score_map[r,c] = s

            # コンタクトシート可視化（原画像を反転してセル配置）
            raw_cell = mask[r*ch:(r+1)*ch, c*cw:(c+1)*cw]
            contact[r*ch:(r+1)*ch, c*cw:(c+1)*cw] = 255 - raw_cell

    imwrite(out_dir/"03_contact.png", contact)

    if debug:
        with open(out_dir/"grid.json", "w", encoding="utf-8") as f:
            import json
            json.dump(grid.tolist(), f, ensure_ascii=False, indent=2)

    return grid, score_map

if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser()
    ap.add_argument("image")
    ap.add_argument("--out", default="read_out")
    ap.add_argument("--templates", default="templates")
    ap.add_argument("--tess", action="store_true", help="Use Tesseract fallback")
    args = ap.parse_args()

    grid, scores = read_sudoku(args.image, out_dir=args.out, template_dir=args.templates, use_tesseract=args.tess)
    print(json.dumps(grid.tolist()))
