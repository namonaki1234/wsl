import json
import re
import sqlite3
import time
from datetime import datetime
from pathlib import Path
from urllib.parse import quote, urljoin, urlparse

import numpy as np
import cshogi
from cshogi import KIF
from cshogi.cli import usi_info_to_score
from cshogi.usi import Engine
from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
from playwright.sync_api import sync_playwright


# =========================
# 設定
# =========================

# このファイルが置かれているディレクトリ
BASE_DIR = Path(__file__).resolve().parent

# SQLite の保存先
DB_NAME = BASE_DIR / "shogi_wars.db"

# やねうら王の実行ファイルと評価関数ディレクトリ
ENGINE_PATH = BASE_DIR / "YaneuraOu" / "source" / "YaneuraOu-by-gcc"
EVAL_DIR = BASE_DIR / "YaneuraOu" / "source" / "eval"

# 取得対象の将棋ウォーズユーザー
DEFAULT_USER_ID = "namonaki1235813"

# 取得したい検索結果ページ
# Shogi-Extend は page=1 が一番新しく、page=2, page=3... と古くなる想定。
TARGET_PAGES = [1, 2, 3]

# page1〜3 の中から、最新順に何局までを対象にするか。
# 例: 3なら最新3局、10なら最新10局。
RECENT_GAME_LIMIT = 2

# Shogi-Extend 関連URL
BASE_URL = "https://www.shogi-extend.com"
SEARCH_URL = f"{BASE_URL}/swars/search?query={{query}}&page={{page}}"
ADAPTER_URL = f"{BASE_URL}/adapter?body={{body}}"

# Playwright 設定
HEADLESS = True
PLAYWRIGHT_TIMEOUT_MS = 20_000

# adapter ページから KIF を取るときのリトライ設定
KIF_POLL_RETRY = 20
KIF_POLL_INTERVAL_MS = 500

# 1局面あたりの思考時間（ミリ秒）
ENGINE_TIME_MS = 500


# =========================
# ログ
# =========================
def log(message: str) -> None:
    """
    時刻つきでログを表示する。
    """
    now = time.strftime("%H:%M:%S")
    print(f"[{now}] {message}", flush=True)


# =========================
# DB
# =========================
def init_db() -> None:
    """
    SQLite DB を初期化する。
    既存DBでも不足列があれば ALTER TABLE で補う。
    """
    log(f"DB初期化: {DB_NAME}")
    conn = sqlite3.connect(DB_NAME)

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS games (
            id TEXT PRIMARY KEY,
            opponent TEXT,
            start_datetime TEXT,
            sfen_list TEXT,
            scores TEXT,
            worst_move INTEGER,
            moves_usi TEXT,
            black_name TEXT,
            white_name TEXT
        )
        """
    )

    cols = [row[1] for row in conn.execute("PRAGMA table_info(games)").fetchall()]

    if "start_datetime" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN start_datetime TEXT")
    if "moves_usi" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN moves_usi TEXT")
    if "black_name" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN black_name TEXT")
    if "white_name" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN white_name TEXT")

    # Streamlit 側で直近順に表示しやすくするためのインデックス。
    # 既に存在していてもエラーにならない。
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_games_start_datetime ON games(start_datetime DESC)"
    )

    conn.commit()
    conn.close()
    log("DB初期化完了")


def get_existing_game_ids() -> set[str]:
    """
    すでに DB に保存済みの game_id 一覧を取得する。
    """
    conn = sqlite3.connect(DB_NAME)
    ids = {row[0] for row in conn.execute("SELECT id FROM games").fetchall()}
    conn.close()
    return ids


def save_game_to_db(
    game_id: str,
    opponent: str,
    start_datetime: str,
    sfen_list: list[str],
    scores: list[int],
    worst_move_idx: int,
    usi_moves: list[str],
    black_name: str,
    white_name: str,
) -> None:
    """
    1局分の解析結果を DB に保存する。
    念のため INSERT OR IGNORE にして、重複IDなら保存を飛ばす。
    """
    conn = sqlite3.connect(DB_NAME)
    conn.execute(
        """
        INSERT OR IGNORE INTO games
        (id, opponent, start_datetime, sfen_list, scores, worst_move, moves_usi, black_name, white_name)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            game_id,
            opponent,
            start_datetime,
            json.dumps(sfen_list, ensure_ascii=False),
            json.dumps(scores, ensure_ascii=False),
            worst_move_idx,
            json.dumps(usi_moves, ensure_ascii=False),
            black_name,
            white_name,
        ),
    )
    conn.commit()
    conn.close()


# =========================
# 文字列処理
# =========================
def normalize_text(text: str | None) -> str:
    """
    改行コードや不可視文字を整える。
    """
    if not text:
        return ""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("\u00a0", " ")
    text = text.replace("\u200b", "")
    return text.strip()


def extract_game_id_from_href(href: str) -> str:
    """
    /swars/battles/xxx/?viewpoint=black のような href から battle ID を抜き出す。
    """
    full_url = urljoin(BASE_URL, href)
    path = urlparse(full_url).path.rstrip("/")
    if not path:
        return ""
    return path.split("/")[-1]


def extract_opponent_from_game_id(game_id: str, user_id: str) -> str:
    """
    game_id から相手名をざっくり推定する。
    例:
    namonaki1235813-opponent-20260422_080046
    """
    parts = game_id.split("-")
    if len(parts) >= 3:
        if parts[0] == user_id:
            return parts[1]
        if parts[1] == user_id:
            return parts[0]
    return "Opponent"


def extract_players_from_kif(kif_text: str):
    """
    KIF ヘッダから先手名と後手名を抜き出す。
    """
    black_match = re.search(r"^先手[：:](.+)$", kif_text, flags=re.MULTILINE)
    white_match = re.search(r"^後手[：:](.+)$", kif_text, flags=re.MULTILINE)

    black_name = black_match.group(1).strip() if black_match else "先手"
    white_name = white_match.group(1).strip() if white_match else "後手"
    return black_name, white_name


def extract_start_datetime_from_kif(kif_text: str) -> str:
    """
    KIF の開始日時を SQLite でソートしやすい文字列に変換する。
    例:
    2026/04/22 08:00:46 -> 2026-04-22 08:00:46
    """
    match = re.search(r"^開始日時[：:](.+)$", kif_text, flags=re.MULTILINE)
    if not match:
        return ""

    raw = match.group(1).strip()

    for fmt in ("%Y/%m/%d %H:%M:%S", "%Y-%m-%d %H:%M:%S"):
        try:
            dt = datetime.strptime(raw, fmt)
            return dt.strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            pass

    return raw


def is_kif_like(text: str) -> bool:
    """
    KIF らしいヘッダと指し手行があるかをざっくり判定する。
    """
    text = normalize_text(text)
    if not text:
        return False

    has_header = any(
        key in text
        for key in (
            "先手：", "後手：", "手合割：", "開始日時：",
            "先手:", "後手:", "手合割:", "開始日時:"
        )
    )
    has_moves = bool(re.search(r"^\s*\d+\s+\S+", text, flags=re.MULTILINE))
    return has_header and has_moves


def trim_to_kif(text: str) -> str:
    """
    ページ全体テキストの中から KIF 部分だけを切り出す。
    """
    text = normalize_text(text)
    if not text:
        return ""

    lines = [line.rstrip() for line in text.split("\n")]
    if not lines:
        return ""

    start_idx = None
    meta_prefixes = (
        "開始日時", "終了日時", "棋戦", "場所", "持ち時間",
        "手合割", "先手", "後手", "戦型", "備考"
    )

    for i, line in enumerate(lines):
        stripped = line.strip()
        if any(stripped.startswith(prefix) for prefix in meta_prefixes):
            start_idx = i
            break
        if re.match(r"^\d+\s+\S+", stripped):
            start_idx = i
            break

    if start_idx is None:
        return ""

    candidate_lines = lines[start_idx:]
    candidate = "\n".join(candidate_lines).strip()
    return candidate if is_kif_like(candidate) else ""


# =========================
# Playwright 補助
# =========================
def save_debug_files(page, prefix: str) -> None:
    """
    失敗時の調査用に HTML / 画像を保存する。
    """
    try:
        png_path = BASE_DIR / f"{prefix}.png"
        page.screenshot(path=str(png_path), full_page=True)
        log(f"DEBUG: スクリーンショット保存: {png_path}")
    except Exception as e:
        log(f"DEBUG: スクリーンショット保存失敗: {e}")

    try:
        html_path = BASE_DIR / f"{prefix}.html"
        html_path.write_text(page.content(), encoding="utf-8")
        log(f"DEBUG: HTML保存: {html_path}")
    except Exception as e:
        log(f"DEBUG: HTML保存失敗: {e}")


def collect_kif_candidates(page) -> list[str]:
    """
    adapter ページ内から KIF らしきテキスト候補を拾う。
    """
    candidates: list[str] = []
    seen: set[str] = set()

    selectors = [
        "textarea",
        "pre",
        "code",
        "[data-clipboard-text]",
        "[data-kif]",
        "[data-kifu]",
        "main",
        "article",
        "body",
    ]

    for selector in selectors:
        locator = page.locator(selector)
        try:
            count = min(locator.count(), 5)
        except Exception:
            continue

        for i in range(count):
            try:
                payload = locator.nth(i).evaluate(
                    """
                    (el) => {
                        const value = ('value' in el && typeof el.value === 'string') ? el.value : '';
                        return {
                            text: el.innerText || el.textContent || '',
                            value: value,
                            clipboard: el.getAttribute('data-clipboard-text') || '',
                            kif: el.getAttribute('data-kif') || '',
                            kifu: el.getAttribute('data-kifu') || ''
                        };
                    }
                    """
                )
            except Exception:
                continue

            for raw in payload.values():
                candidate = trim_to_kif(raw)
                if candidate and candidate not in seen:
                    seen.add(candidate)
                    candidates.append(candidate)

    candidates.sort(key=len, reverse=True)
    return candidates


def click_possible_reveal_buttons(page) -> None:
    """
    棋譜を開きそうなボタン類を押して再試行する。
    """
    texts = ["棋譜", "KIF", "Copy", "コピー", "表示", "開く"]

    for text in texts:
        for tag in ("button", "a"):
            locator = page.locator(tag, has_text=text)
            try:
                count = min(locator.count(), 2)
            except Exception:
                continue

            for i in range(count):
                try:
                    locator.nth(i).click(timeout=1500)
                    page.wait_for_timeout(500)
                    return
                except Exception:
                    continue


def collect_recent_battle_hrefs(
    page,
    user_id: str,
    target_pages: list[int],
    recent_game_limit: int,
) -> list[str]:
    """
    page=1 -> page=2 -> page=3 の順に検索ページを見て，
    画面に出ている順番のまま battle href を集める。

    ここで recent_game_limit 件に絞るので，
    「page1〜3にある対局のうち，最新順にN局だけ」が対象になる。

    注意:
    - DBに既にあるかどうかは，この段階では見ない。
    - つまり「直近N局を対象にして，そのうち未保存だけ保存する」動きになる。
    """
    recent_hrefs: list[str] = []
    seen_game_ids: set[str] = set()

    for page_no in target_pages:
        if len(recent_hrefs) >= recent_game_limit:
            break

        search_url = SEARCH_URL.format(query=quote(user_id, safe=""), page=page_no)
        log(f"DEBUG: 検索ページ取得 page={page_no}: {search_url}")
        page.goto(search_url, wait_until="domcontentloaded")

        try:
            page.wait_for_load_state("networkidle", timeout=10_000)
        except PlaywrightTimeoutError:
            pass

        detail_links = page.locator("a", has_text="詳細")
        link_count = detail_links.count()
        log(f"DEBUG: page={page_no} の詳細リンク数 = {link_count}")

        for i in range(link_count):
            if len(recent_hrefs) >= recent_game_limit:
                break

            href = detail_links.nth(i).get_attribute("href")
            if not href or href == "/":
                continue

            game_id = extract_game_id_from_href(href)
            if not game_id:
                continue

            # ページをまたいだ重複対策。
            if game_id in seen_game_ids:
                continue

            seen_game_ids.add(game_id)
            recent_hrefs.append(href)
            log(f"TARGET: {len(recent_hrefs)}/{recent_game_limit} page={page_no} -> {game_id}")

    return recent_hrefs


def scrape_single_game_from_href(page, battle_href: str, user_id: str):
    """
    1局分だけ adapter に飛んで KIF を取得する。
    """
    game_id = extract_game_id_from_href(battle_href)
    if not game_id:
        return None

    battle_url = urljoin(BASE_URL, battle_href)
    adapter_url = ADAPTER_URL.format(body=quote(battle_url, safe=""))

    log(f"DEBUG: KIF取得中 game_id={game_id}")
    page.goto(adapter_url, wait_until="domcontentloaded")

    try:
        page.wait_for_load_state("networkidle", timeout=10_000)
    except PlaywrightTimeoutError:
        pass

    kif_text = ""
    for _ in range(KIF_POLL_RETRY):
        candidates = collect_kif_candidates(page)
        if candidates:
            kif_text = candidates[0]
            break
        page.wait_for_timeout(KIF_POLL_INTERVAL_MS)

    if not kif_text:
        click_possible_reveal_buttons(page)
        for _ in range(KIF_POLL_RETRY):
            candidates = collect_kif_candidates(page)
            if candidates:
                kif_text = candidates[0]
                break
            page.wait_for_timeout(KIF_POLL_INTERVAL_MS)

    if not kif_text:
        log(f"Error: 棋譜データの取得に失敗しました: {game_id}")
        save_debug_files(page, f"debug_adapter_error_{game_id}")
        return None

    opponent = extract_opponent_from_game_id(game_id, user_id)
    return {
        "game_id": game_id,
        "opponent": opponent,
        "kif_text": kif_text,
    }


def scrape_recent_target_games(
    user_id: str,
    target_pages: list[int],
    recent_game_limit: int,
    existing_ids: set[str],
) -> list[dict]:
    """
    page1〜3から最新順に recent_game_limit 件だけを対象にし，
    その中でDBに存在しない局だけ KIF を取得して返す。
    """
    new_games: list[dict] = []

    with sync_playwright() as p:
        log("Playwright 起動")
        browser = p.chromium.launch(headless=HEADLESS)
        context = browser.new_context(
            user_agent=(
                "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/123.0.0.0 Safari/537.36"
            ),
            viewport={"width": 1280, "height": 1200},
        )
        context.set_default_timeout(PLAYWRIGHT_TIMEOUT_MS)
        page = context.new_page()

        try:
            # まず page1 -> page2 -> page3 の順にリンクだけを集める。
            # adapterページへ移動すると検索ページのlocatorが使いづらくなるため，先に一覧化する。
            recent_hrefs = collect_recent_battle_hrefs(
                page=page,
                user_id=user_id,
                target_pages=target_pages,
                recent_game_limit=recent_game_limit,
            )

            log(f"DEBUG: 今回の対象局数 = {len(recent_hrefs)}")

            for idx, href in enumerate(recent_hrefs, start=1):
                game_id = extract_game_id_from_href(href)
                if not game_id:
                    continue

                # 既存局なら，解析も保存もしない。
                if game_id in existing_ids:
                    log(f"SKIP: 既存対局なので保存しない {idx}/{len(recent_hrefs)} -> {game_id}")
                    continue

                log(f"DEBUG: 新規対局を取得中 {idx}/{len(recent_hrefs)} -> {game_id}")
                item = scrape_single_game_from_href(page, href, user_id)

                if item:
                    new_games.append(item)
                    # 同一実行内での重複防止。
                    existing_ids.add(game_id)

        finally:
            browser.close()
            log("Playwright 終了")

    return new_games


# =========================
# KIF 解析
# =========================
def parse_moves_from_kif(kif_text: str):
    """
    KIF文字列を cshogi で解析し、内部 move 配列を返す。
    """
    kif_obj = KIF.Parser.parse_str(kif_text)
    if hasattr(kif_obj, "moves"):
        return kif_obj.moves
    raise RuntimeError("KIFのパースに失敗しました。")


# =========================
# USI リスナー
# =========================
class Listener:
    """
    USI エンジンから流れてくる info 行を保持する簡易リスナー。
    """
    def __init__(self):
        self.info = ""
        self.bestmove = ""

    def __call__(self, line):
        self.info = self.bestmove
        self.bestmove = line


# =========================
# 1局解析
# =========================
def analyze_one_game(engine, listener, game_data: dict) -> None:
    """
    1局分を解析して DB に保存する。
    """
    game_id = game_data["game_id"]
    opponent = game_data["opponent"]
    kif_text = game_data["kif_text"]

    log(f"DEBUG: 解析開始 game_id={game_id}")

    moves = parse_moves_from_kif(kif_text)
    black_name, white_name = extract_players_from_kif(kif_text)
    start_datetime = extract_start_datetime_from_kif(kif_text)
    usi_moves = [cshogi.move_to_usi(move) for move in moves]

    board = cshogi.Board()
    sfen_list = [board.sfen()]
    scores = [0]

    engine.usinewgame(listener=listener)

    for i, move in enumerate(moves, start=1):
        board.push(move)
        sfen_list.append(board.sfen())

        try:
            listener.info = ""
            listener.bestmove = ""

            engine.position(sfen=board.sfen())
            engine.go(byoyomi=ENGINE_TIME_MS, listener=listener)

            raw_score = usi_info_to_score(listener.info)
            if raw_score is None:
                score = scores[-1]
            else:
                score = int(raw_score)

        except Exception as e:
            log(f"Warning: {game_id} / {i} 手目の評価値取得に失敗: {e}")
            score = scores[-1]

        # 先手視点で揃える
        if board.turn == cshogi.WHITE:
            score = -score

        scores.append(score)

    diffs = [abs(scores[i] - scores[i - 1]) for i in range(1, len(scores))]
    worst_move_idx = int(np.argmax(diffs)) + 1 if diffs else 0

    save_game_to_db(
        game_id=game_id,
        opponent=opponent,
        start_datetime=start_datetime,
        sfen_list=sfen_list,
        scores=scores,
        worst_move_idx=worst_move_idx,
        usi_moves=usi_moves,
        black_name=black_name,
        white_name=white_name,
    )

    log(f"NEW SAVED: {game_id} / start={start_datetime}")


# =========================
# 全ページ解析
# =========================
def analyze_and_save_recent_games(user_id: str) -> None:
    """
    page1〜3から最新順にRECENT_GAME_LIMIT局だけ見て、未保存の局だけ追加する。
    """
    log("analyze_and_save_recent_games 開始")

    existing_ids = get_existing_game_ids()
    log(f"DEBUG: 既存対局数 = {len(existing_ids)}")

    new_games = scrape_recent_target_games(user_id, TARGET_PAGES, RECENT_GAME_LIMIT, existing_ids)

    if not new_games:
        log("新しく保存すべき対局はありませんでした。")
        return

    if not ENGINE_PATH.exists():
        log(f"Error: エンジンが見つかりません: {ENGINE_PATH}")
        return

    if not EVAL_DIR.exists():
        log(f"Error: EvalDir が見つかりません: {EVAL_DIR}")
        return

    engine = None
    listener = Listener()

    try:
        log(f"DEBUG: エンジン起動開始: {ENGINE_PATH}")
        engine = Engine(str(ENGINE_PATH), connect=True)

        engine.setoption("EvalDir", str(EVAL_DIR), listener=listener)
        engine.setoption("BookFile", "no_book", listener=listener)
        engine.isready(listener=listener)

        log("DEBUG: engine.isready 完了")
        log(f"DEBUG: 新規保存対象の局数 = {len(new_games)}")

        for idx, game_data in enumerate(new_games, start=1):
            log(f"=== {idx}/{len(new_games)} 局目の解析開始 ===")
            analyze_one_game(engine, listener, game_data)

    finally:
        if engine is not None:
            try:
                engine.quit(listener=listener)
            except Exception:
                pass

    log("analyze_and_save_recent_games 完了")


if __name__ == "__main__":
    init_db()
    analyze_and_save_recent_games(DEFAULT_USER_ID)