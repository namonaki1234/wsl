# make_templates_stronger.py
import argparse, uuid
from pathlib import Path
import numpy as np, cv2

def pipelines(gray):
    g1 = cv2.GaussianBlur(gray, (3,3), 0)
    _, th1 = cv2.threshold(g1, 0, 255, cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)

    th2 = cv2.adaptiveThreshold(gray,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                cv2.THRESH_BINARY_INV,11,2)

    g3 = cv2.bilateralFilter(gray, 7, 40, 40)
    _, th3 = cv2.threshold(g3, 0, 255, cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)

    return [th1, th2, th3]

def extract_roi_any(cell, min_ratio=0.004, max_ratio=0.75):
    gray = cv2.cvtColor(cell, cv2.COLOR_BGR2GRAY) if cell.ndim==3 else cell
    for th in pipelines(gray):
        th = cv2.morphologyEx(th, cv2.MORPH_OPEN, np.ones((2,2),np.uint8),1)
        cnts,_ = cv2.findContours(th, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not cnts: continue
        h,w = th.shape; area=h*w
        cnt = max(cnts, key=cv2.contourArea)
        r = cv2.contourArea(cnt)/area
        if r<min_ratio or r>max_ratio: continue
        x,y,wc,hc = cv2.boundingRect(cnt)
        roi = th[y:y+hc, x:x+wc]
        s = max(wc,hc)
        canv = np.zeros((s,s),np.uint8)
        canv[(s-hc)//2:(s-hc)//2+hc, (s-wc)//2:(s-wc)//2+wc] = roi
        roi = cv2.resize(canv,(28,28),cv2.INTER_AREA)
        if np.mean(roi)<127: roi = 255-roi
        return roi
    return None

def ocr_tess(roi28):
    try:
        import pytesseract
    except Exception:
        return None
    big = cv2.resize(roi28,(96,96),cv2.INTER_NEAREST)
    cfg = r'--psm 10 --oem 3 -c tessedit_char_whitelist=123456789'
    s = pytesseract.image_to_string(big, config=cfg).strip()
    return int(s) if (len(s)==1 and s.isdigit() and s!='0') else None

def main(raw_dir, mid_dir, final_dir):
    raw_dir, mid_dir, final_dir = map(Path,[raw_dir, mid_dir, final_dir])
    mid_dir.mkdir(parents=True, exist_ok=True); final_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(raw_dir.glob("*.png"))
    print(f"[info] raw cells: {len(files)} from {raw_dir}")

    per_digit = {d: [] for d in range(1,10)}
    for p in files:
        cell = cv2.imread(str(p))
        if cell is None: continue
        roi = extract_roi_any(cell)
        if roi is None: continue
        d = ocr_tess(roi)
        if d is None: continue
        cv2.imwrite(str(mid_dir/f"{d}_{uuid.uuid4().hex[:6]}.png"), roi)
        per_digit[d].append(roi)

    for d in range(1,10):
        arr = per_digit[d]
        if not arr:
            print(f"[warn] digit {d}: no samples")
            continue
        stack = np.stack(arr).astype(np.float32)
        avg = np.clip(stack.mean(0),0,255).astype(np.uint8)
        cv2.imwrite(str(final_dir/f"{d}.png"), avg)
        print(f"[ok] digit {d}: {len(arr)} samples -> {final_dir}/{d}.png")

if __name__=="__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw",   default="./raw_img")
    ap.add_argument("--mid",   default="./templates_auto")
    ap.add_argument("--final", default="./templates")
    args = ap.parse_args()
    main(args.raw, args.mid, args.final)
