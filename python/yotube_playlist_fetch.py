from googleapiclient.discovery import build
import gspread
import os
from dotenv import load_dotenv

load_dotenv()  # .envを読み込む

# YouTube APIの設定
api_key = os.getenv("API_KEY")
youtube = build("youtube", "v3", developerKey=api_key)

playlist_id = "PLbKGrU2uWtjbuakvgNrFz9a0nniTOe9BN"  # 再生リストID (WLは「Watch Later」の略)

# 再生リスト取得
video_items = []
nextPageToken = None

while True:
    res = youtube.playlistItems().list(
        part="snippet",
        playlistId=playlist_id,
        maxResults=50,
        pageToken=nextPageToken
    ).execute()

    video_items += res["items"]

    nextPageToken = res.get("nextPageToken")
    if not nextPageToken:
        break

# 動画タイトル、チャンネル名、チャンネルID、URLを抽出
video_data = []
for item in video_items:
    title = item["snippet"]["title"]
    video_id = item["snippet"]["resourceId"]["videoId"]
    url = f"https://www.youtube.com/watch?v={video_id}"
    channel_title = item["snippet"]["channelTitle"]
    channel_id = item["snippet"]["channelId"]
    video_data.append([title, channel_title, channel_id, url])

# Google Sheets に書き込み
gc = gspread.service_account(filename="/root/wsl/python/credentials.json")
# sh = gc.open("youtube_WL_fetch")
# スプレッドシートIDで開く
spreadsheet_id = "1344XbxD0uQ6gPLGf3Fla1Ur80jpejnc7TaS8jbWCqVA"
sh = gc.open_by_key(spreadsheet_id)
worksheet = sh.sheet1

# ヘッダー行も修正
worksheet.update(
    'A1',
    [["Title", "ChannelName", "ChannelID", "URL"]] + video_data
)

print("書き込み完了！")
