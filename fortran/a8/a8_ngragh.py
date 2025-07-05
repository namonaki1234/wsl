import numpy as np
import matplotlib.pyplot as plt
import japanize_matplotlib

# データ読み込み
data = np.loadtxt('a8_y_plus_u_plus.csv')
y_plus_data = data[:, 0]
u_plus_data = data[:, 1]

# 理論曲線生成
y_plus_theory = np.logspace(0, 4, 300)
u_plus_theory = np.piecewise(
    y_plus_theory,
    [y_plus_theory <= 11.63, y_plus_theory > 11.63],
    [lambda y: y,
     lambda y: (1/0.41)*np.log(y) + 5.2]
)

# グラフ描画
fig, ax = plt.subplots(figsize=(6,4))
ax.semilogx(y_plus_data, u_plus_data, 'r.', markersize=4, label='計算値')
ax.semilogx(y_plus_theory, u_plus_theory, 'b-', linewidth=1, label='理論値')
ax.set_xlabel(r'$y^+$', fontsize=14)
ax.set_ylabel(r'$u^+$', fontsize=14)
ax.legend()

# x軸の範囲を指定
ax.set_xlim(1, 1e4)

# 目盛線を内向きにする
ax.tick_params(which='both', direction='in')

# plt.grid(which='both', linestyle='--', linewidth=0.5)
ax.grid(False)

plt.tight_layout()

# 画像保存
plt.savefig('yplus_uplus.png', dpi=300)

# 画面表示
plt.show()
