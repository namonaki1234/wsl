# harvest_unlabeled.py
from pathlib import Path
import cv2, numpy as np

RAW = Path("./raw_img"); TO = Path("./to_label")
TO.mkdir(exist_ok=True)

def roi(cell):
    g = cv2.cvtColor(cell, cv2.COLOR_BGR2GRAY)
    _, th = cv2.threshold(g, 0, 255, cv2.THRESH_BINARY_INV+cv2.THRESH_OTSU)
    th = cv2.morphologyEx(th, cv2.MORPH_OPEN, np.ones((2,2),np.uint8),1)
    cnts,_ = cv2.findContours(th, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not cnts: return None
    h,w = th.shape; area=h*w
    cnt = max(cnts, key=cv2.contourArea)
    x,y,wc,hc = cv2.boundingRect(cnt)
    r = th[y:y+hc, x:x+wc]
    s=max(wc,hc); canv=np.zeros((s,s),np.uint8)
    canv[(s-hc)//2:(s-hc)//2+hc, (s-wc)//2:(s-wc)//2+wc]=r
    r=cv2.resize(canv,(28,28),cv2.INTER_AREA)
    if np.mean(r)<127: r=255-r
    return r

for p in sorted(RAW.glob("*.png")):
    img = cv2.imread(str(p))
    r = roi(img)
    if r is None: continue
    cv2.imwrite(str(TO/f"u_{p.stem}.png"), r)  # u_ は unlabeled の意
print("DONE: to_label/ に未ラベルROIを保存しました。")
