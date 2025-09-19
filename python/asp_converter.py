# Auto Screencapture and PDF Converter, version 1.3
# 電子書籍や資料などを自動で一枚ずつスクリーンショットを撮り，PDF化，PDF圧縮を行うプログラム

##################################################
# PDFファイル名，頭文字
##################################################

pdf_header = "test" # 生成するPDFファイルとフォルダの頭文字に使用する

##################################################
# モード切り替え
##################################################

# キャプチャ範囲の数値：手入力で指定する場合False，マウスカーソルで指定する場合True
coordinate_setup = True
# ページ移動：矢印キーならキーを指定．マウスの左クリックなら'click'
manual_mode = 'right' # 例:'right','down'
# キャプチャするウィンドウを用意するのに何秒時間が欲しいか
wait_for_setup = 3 # 例:1
# PDF出力する場合はTrue, しない場合はFalse
output_pdf = True
# PDFファイル生成後，PDFファイルを自動で開く場合はTrue, しない場合はFalse
open_pdf = True
# PDFファイル生成後，GhostScriptを用いて圧縮する場合はTrue, しない場合はFalse
compress_pdf = True #時間はかかるが，推奨．

##################################################
# 自動スクリーンキャプチャの設定
##################################################

# スクショ取る枚数
page = 1

# 見開きのページを分割して1ページずつPDF化する場合はTrue，しない場合はFalse
split = True # 見開きのページを分割する場合，解像度が落ちるので注意
# 表紙の枚数（見開きの本のとき，最初に中央寄りのページが何枚か）
dont_split = 2

# スクショ間隔(秒)
span = 1

# 出力画像　ファイル頭文字
h_filename = "picture"

# キャプチャ範囲の指定，手入力
if coordinate_setup == False:
    # キャプチャ範囲：左上座標
    position1 = (151.07421875, 142.58984375)
    # キャプチャ範囲：右下座様
    position2 = (1602.8125, 1166.828125)
    # ページ間移動：マウスクリック座標
    position3 = (0, 0)





#############################プログラム開始(以降いじらない)###########################
import time,os,platform,glob,subprocess,datetime,threading #標準ライブラリ
import img2pdf,cv2,pyautogui
from pathlib import Path
from pynput import mouse,keyboard
from natsort import natsorted
from PIL import ImageGrab

print('\nーーーーーASP Converter開始ーーーーー\n')
t0 = time.time()

# 出力フォルダ作成(フォルダ名：頭文字_年月日時分秒)
root = Path.cwd()
folder_name = Path(root, "output", pdf_header + "_" + str(datetime.datetime.now().strftime("%Y%m%d%H%M%S")))
Path(root, "output").mkdir(exist_ok=True)
Path(folder_name).mkdir(exist_ok=True)

if coordinate_setup == True:
    #########################
    # キャプチャ座標計測
    #########################

    #座標を取得
    print("＜画面上の座標を取得します＞")
    if split == True:
        print("分割する場合においても見開き2ページ分の範囲を指定してください")
    print('マウスカーソルを\nキャプチャ範囲の左上座標に合わせて1を')
    print('キャプチャ範囲の右下座標に合わせて2を')
    if manual_mode == 'click':print('マウスクリックでページ間を移動する際に，クリックする座標に合わせて3を')
    print("押してください．")
    print("全て完了したらenterを押してください．")
    print("座標は何度でも上書き可能です．")

    # 押されたキーを確認
    def press0(key):
        global position1, position2, position3
        # print('キー {0} が押されました\n'.format(key))
        if 'char' in dir(key):
            if key.char == '1':
                position1=mouse.Controller().position
                print(' 左上座標：{}'.format(position1))
            elif key.char == '2':
                position2=mouse.Controller().position
                print(' 右下座標：{}'.format(position2))
            elif manual_mode == 'click' and key.char == '3':
                position3=mouse.Controller().position
                print(' クリック座標：{}'.format(position3))
        if key == keyboard.Key.enter:
            if manual_mode == 'click' and 'position1' in globals() and 'position2' in globals() and 'position3' in globals():
                print('全ての座標を取得しました．\n')
                return False
            elif manual_mode != 'click' and 'position1' in globals() and 'position2' in globals():
                print('全ての座標を取得しました．\n')
                return False
            else:
                if 'position1' not in globals():print(" 左上座標を取得していません．")
                if 'position2' not in globals():print(" 右上座標を取得していません．")
                if manual_mode == 'click' and 'position3' not in globals():
                    print(" クリックする座標を取得していません．")
                print("再度お試しください")
    with keyboard.Listener(on_press=press0) as keyboard_listener:keyboard_listener.join()

# この間にスクショするウィンドウをアクティブにする
if wait_for_setup > 0:
    print("{}秒間待機".format(wait_for_setup))
    time.sleep(wait_for_setup)

#########################
# スクリーンショット自動処理
#########################

# ページ数分スクリーンショットをとる
print("＜スクリーンキャプチャ開始＞")

x1=position1[0];y1=position1[1];x2=position2[0];y2=position2[1]
paused = False

def capture(imagesize):
    if platform.system() == "Darwin":
        s=ImageGrab.grab()
        fullsize = pyautogui.size()
        left = round(imagesize[0]*2)
        top = round(imagesize[1]*2)
        right = round(imagesize[2]*2)
        bottom = round(imagesize[3]*2)
        s = s.crop((left, top, right, bottom))
    else:
        s = ImageGrab.grab(imagesize)
    return s

def press1(key): # ESCキーで途中で中断
    global paused
    if key == keyboard.Key.esc:
        paused = True
        return False

def screencapture():
    i=1
    for p in range(int(page)):
        # 出力ファイル名(頭文字_連番.png)
        out_filename1 = h_filename + "_" + str(i).zfill(4) + '.png'
        out_filename2 = h_filename + "_" + str(i+1).zfill(4) + '.png'

        # スクリーンショット取得・保存処理
        if split == False:
            # キャプチャ範囲： 左上のx座標, 左上のy座標, 右下のx座標, 右下のy座標
            cbox = (round(x1),round(y1),round(x2),round(y2))
            screenshot = capture(cbox)
            # 出力パス： 出力フォルダ名 / 出力ファイル名
            screenshot.save(Path(folder_name, out_filename1),optimize=True);i+=1
        elif split == True and i < dont_split+1:
            cbox = (round(x1+(x2-x1)/4),round(y1),round(x1+(x2-x1)*3/4),round(y2))
            screenshot = capture(cbox)
            screenshot.save(Path(folder_name, out_filename1),optimize=True);i+=1
        elif split == True and i > dont_split:
            cbox = (round(x1),round(y1),round(x1+(x2-x1)/2),round(y2))
            screenshot = capture(cbox)
            screenshot.save(Path(folder_name, out_filename1),optimize=True);
            cbox = (round(x1+(x2-x1)/2),round(y1),round(x2),round(y2))
            screenshot = capture(cbox)
            screenshot.save(Path(folder_name, out_filename2),optimize=True);i+=2

        # 次のページへの操作
        if p < page-1:
            if manual_mode == 'click':
                mouse.Controller().position = position3
                mouse.Controller().click(mouse.Button.left,1)
            else:
                line = 'keyboard.Controller().tap(keyboard.Key.'+manual_mode+')'
                exec(line)

        # 次のスクリーンショットまで待機
        time.sleep(span)

        #ESCキーで途中で中断
        if paused == True:
            break

thread = threading.Thread(target=screencapture)
thread.start()
with keyboard.Listener(on_press=press1) as listener:
    thread.join()

print("＜スクリーンキャプチャ終了＞")

#########################
# PDF化の処理
#########################

if output_pdf == True:
    print("＜PDF化開始＞")
    # 画像一覧を取得
    lists = list(glob.glob(str(Path(folder_name, "*.png"))))

    # 透過チャンネルがある場合に削除する処理
    for filename in lists:
        img = cv2.imread(filename,cv2.IMREAD_UNCHANGED)
        if img.shape[2] == 4:
            img2 = cv2.imread(filename,cv2.IMREAD_COLOR)
            cv2.imwrite(filename, img2)

    # PDFファイルを出力
    filename = pdf_header + "_" + str(datetime.datetime.now().strftime("%Y%m%d%H%M%S")) + ".pdf"
    output_file = Path(root, "output", filename)
    with open(output_file,"wb") as f:
        f.write(img2pdf.convert([str(i) for i in natsorted(lists) if ".png" in i]))
    if open_pdf == True and compress_pdf == False:
        if platform.system() == "Windows": subprocess.Popen([output_file], shell=True)
        else: subprocess.run(['open', output_file], check=True)
    print("＜PDF化終了＞")

#########################
# PDF圧縮の処理
#########################

if compress_pdf == True:
    print("＜PDF圧縮開始＞")
    if platform.system() == "Windows":
        gs_command = "gswin64c"
    else:
        gs_command = 'gs'
    output_file_compressed = str(output_file).replace('.pdf','_compressed.pdf')
    arg1= '-sOutputFile=' + output_file_compressed
    p = subprocess.run([gs_command, '-sDEVICE=pdfwrite', '-dCompatibilityLevel=1.4', '-dPDFSETTINGS=/ebook', '-dNOPAUSE', '-dBATCH',  '-dQUIET', arg1, output_file],stdout=subprocess.PIPE)
    if open_pdf == True and compress_pdf == True:
        op1 = Path(str(output_file_compressed))
        if platform.system() == "Windows": subprocess.Popen([op1], shell=True)
        else: subprocess.run(['open', op1], check=True)
    print("＜PDF圧縮終了＞")

#########################
# 出力結果
#########################

original_size = os.path.getsize(output_file)

print(f'\n出力ファイル容量: {original_size/1024/1024:.3f} MB')
if compress_pdf == True:
    compressed_size = os.path.getsize(output_file_compressed)
    print(f'圧縮後: {compressed_size/1024/1024:.3f} MB\n')

t1 = time.time()
seconds = round(t1-t0,2)
minutes, seconds =divmod(seconds,60)
hours, minutes=divmod(minutes,60)
print(f"経過時間:{hours}時間{minutes}分{seconds:.2f}秒")

print('ーーーーーASP Converter終了ーーーーー\n')