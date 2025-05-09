import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === fldファイルの読み取り関数（修正版） ===
def read_fld_metadata(fld_filename):
    metadata = {}
    with open(fld_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                if key == "label":
                    metadata["labels"] = value.split()
                else:
                    metadata[key] = value
    return metadata

# === メタデータ読み込み ===
fld_file = "mac_method.fld"
dat_file = "mac_method.dat"
meta = read_fld_metadata(fld_file)

# === 必要な情報を取得 ===
nx = int(meta['dim1'])
ny = int(meta['dim2'])
veclen = int(meta['veclen'])
labels = meta['labels']

# === データ読み込みと整形 ===
data = np.loadtxt(dat_file)
variables = {label: data[:, i].reshape((ny, nx)) for i, label in enumerate(labels)}

# === 座標グリッド作成（仮の物理スケール）===
x = np.linspace(0, 11.0, nx)
y = np.linspace(0, 2.0, ny)
X, Y = np.meshgrid(x, y)

# === コンター図の作成と保存（ここでは圧力 'u'） ===
plt.figure(figsize=(10, 4))
contour = plt.contourf(X, Y, variables['u'], levels=50, cmap='viridis')
plt.colorbar(contour, label='Value')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour of Backstep Flow_u')
plt.axis('scaled')
plt.tight_layout()
plt.savefig("contour_u.png", dpi=300)
# plt.show()
plt.close()

# === コンター図の作成と保存（ここでは圧力 'v'） ===
plt.figure(figsize=(10, 4))
contour = plt.contourf(X, Y, variables['v'], levels=50, cmap='viridis')
plt.colorbar(contour, label='Value')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour of Backstep Flow_v')
plt.axis('scaled')
plt.tight_layout()
plt.savefig("contour_v.png", dpi=300)
plt.close()

# === コンター図の作成と保存（ここでは圧力 'p'） ===
plt.figure(figsize=(10, 4))
contour = plt.contourf(X, Y, variables['p'], levels=50, cmap='viridis')
plt.colorbar(contour, label='Value')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour of Backstep Flow_p')
plt.axis('scaled')
plt.tight_layout()
plt.savefig("contour_p.png", dpi=300)
plt.close()
