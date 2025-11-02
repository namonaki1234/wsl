import os
import glob

input_dir = "/root/wsl/fortran/c_grid/abe_GridFlow/sjis_files"
output_dir = "/root/wsl/fortran/c_grid/abe_GridFlow/utf8_files"
os.makedirs(output_dir, exist_ok=True)

for file_path in glob.glob(os.path.join(input_dir, "*.f90")):
    file_name = os.path.basename(file_path)
    output_path = os.path.join(output_dir, file_name)

    with open(file_path, "r", encoding="shift_jis", errors="ignore") as f_in:
        content = f_in.read()
    with open(output_path, "w", encoding="utf-8") as f_out:
        f_out.write(content)

    print(f"✅ {file_name} をUTF-8に変換しました")
