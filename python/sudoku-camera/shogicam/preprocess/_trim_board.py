import cv2
import numpy as np
from ._detect_corners import detect_corners

BASE_SIZE = 32
# def trim_board(img, corners):
#     w = BASE_SIZE * 15
#     h = BASE_SIZE * 15
#     transform = cv2.getPerspectiveTransform(np.float32(corners), np.float32([[0, 0], [w, 0], [w, h], [0, h]]))
#     normed = cv2.warpPerspective(img, transform, (w, h))
#     return normed

# def trim_board(img, corners, w=288, h=288):
#     # cornersをnumpy配列に強制変換
#     corners = np.array(corners, dtype=np.float32).reshape(4, 2)

#     transform = cv2.getPerspectiveTransform(
#         corners,
#         np.float32([[0, 0], [w, 0], [w, h], [0, h]])
#     )
#     dst = cv2.warpPerspective(img, transform, (w, h))
#     return dst

# _trim_board.py
def trim_board(img, corners, w=288, h=288):
    # (points, score) 形式にも対応
    if isinstance(corners, tuple) and len(corners) == 2:
        corners = corners[0]

    corners = np.asarray(corners, dtype=np.float32).reshape(4, 2)
    M = cv2.getPerspectiveTransform(
        corners,
        np.float32([[0, 0], [w, 0], [w, h], [0, h]])
    )
    return cv2.warpPerspective(img, M, (w, h))

