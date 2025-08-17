# build_templates_from_labeled.py
from pathlib import Path
import numpy as np, cv2, re

SRC = Path("./to_label")
OUT = Path("./templates")
OUT.mkdir(exist_ok=True)

per_digit = {d: [] for d in range(1,10)}
for p in SRC.glob("*.png"):
    m = re.match(r"([1-9])_", p.name)
    if not m: continue
    d = int(m.group(1))
    im = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
    if im is None: continue
    per_digit[d].append(im)

for d in range(1,10):
    arr = per_digit[d]
    if not arr:
        print(f"[warn] digit {d}: no samples")
        continue
    stack = np.stack(arr).astype(np.float32)
    avg = np.clip(stack.mean(0),0,255).astype(np.uint8)
    cv2.imwrite(str(OUT/f"{d}.png"), avg)
    print(f"[ok] digit {d}: {len(arr)} samples -> templates/{d}.png")
