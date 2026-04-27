import json
import re
import sqlite3
import time
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
BASE_DIR = Path(__file__).resolve().parent
DB_NAME = BASE_DIR / "shogi_wars.db"

ENGINE_PATH = BASE_DIR / "YaneuraOu" / "source" / "YaneuraOu-by-gcc"
EVAL_DIR = BASE_DIR / "YaneuraOu" / "source" / "eval"

DEFAULT_USER_ID = "namonaki1235813"

BASE_URL = "https://www.shogi-extend.com"
SEARCH_URL = f"{BASE_URL}/swars/search?query={{query}}"
ADAPTER_URL = f"{BASE_URL}/adapter?body={{body}}"

HEADLESS = True
PLAYWRIGHT_TIMEOUT_MS = 20_000
KIF_POLL_RETRY = 20
KIF_POLL_INTERVAL_MS = 500

ENGINE_TIME_MS = 500


# =========================
# ログ
# =========================
def log(message: str) -> None:
    now = time.strftime("%H:%M:%S")
    print(f"[{now}] {message}", flush=True)


# =========================
# DB
# =========================
def init_db() -> None:
    log(f"DB初期化: {DB_NAME}")
    conn = sqlite3.connect(DB_NAME)

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS games (
            id TEXT PRIMARY KEY,
            opponent TEXT,
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

    if "moves_usi" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN moves_usi TEXT")
    if "black_name" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN black_name TEXT")
    if "white_name" not in cols:
        conn.execute("ALTER TABLE games ADD COLUMN white_name TEXT")

    conn.commit()
    conn.close()
    log("DB初期化完了")


# =========================
# 文字列処理
# =========================
def normalize_text(text: str | None) -> str:
    if not text:
        return ""
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("\u00a0", " ")
    text = text.replace("\u200b", "")
    return text.strip()


def extract_game_id_from_href(href: str) -> str:
    full_url = urljoin(BASE_URL, href)
    path = urlparse(full_url).path.rstrip("/")
    if not path:
        return ""
    return path.split("/")[-1]


def extract_opponent_from_game_id(game_id: str, user_id: str) -> str:
    parts = game_id.split("-")
    if len(parts) >= 3:
        if parts[0] == user_id:
            return parts[1]
        if parts[1] == user_id:
            return parts[0]
    return "Opponent"


def extract_players_from_kif(kif_text: str):
    black_match = re.search(r"^先手[：:](.+)$", kif_text, flags=re.MULTILINE)
    white_match = re.search(r"^後手[：:](.+)$", kif_text, flags=re.MULTILINE)

    black_name = black_match.group(1).strip() if black_match else "先手"
    white_name = white_match.group(1).strip() if white_match else "後手"
    return black_name, white_name


def is_kif_like(text: str) -> bool:
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

        if count > 0:
            log(f"DEBUG: selector '{selector}' に {count} 件ヒット")

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
            except Exception as e:
                log(f"DEBUG: selector '{selector}' index {i} evaluate失敗: {e}")
                continue

            for key, raw in payload.items():
                candidate = trim_to_kif(raw)
                if candidate and candidate not in seen:
                    seen.add(candidate)
                    candidates.append(candidate)
                    log(
                        f"DEBUG: KIF候補追加 selector={selector} index={i} field={key} "
                        f"len={len(candidate)}"
                    )

    candidates.sort(key=len, reverse=True)
    return candidates


def click_possible_reveal_buttons(page) -> None:
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
                    log(f"DEBUG: {tag} '{text}' をクリック試行 index={i}")
                    locator.nth(i).click(timeout=1500)
                    page.wait_for_timeout(500)
                    log(f"DEBUG: {tag} '{text}' をクリック成功")
                    return
                except Exception as e:
                    log(f"DEBUG: {tag} '{text}' クリック失敗 index={i}: {e}")
                    continue


def scrape_with_playwright(user_id: str):
    start_all = time.perf_counter()

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
            log(f"Searching for {user_id} on Shogi-Extend...")
            search_url = SEARCH_URL.format(query=quote(user_id, safe=""))
            log(f"DEBUG: search URL: {search_url}")

            t0 = time.perf_counter()
            page.goto(search_url, wait_until="domcontentloaded")
            log(f"DEBUG: search page goto完了 ({time.perf_counter() - t0:.2f}秒)")

            try:
                t1 = time.perf_counter()
                page.wait_for_load_state("networkidle", timeout=10_000)
                log(f"DEBUG: search page networkidle完了 ({time.perf_counter() - t1:.2f}秒)")
            except PlaywrightTimeoutError:
                log("DEBUG: search page networkidle タイムアウト")

            detail_link = page.locator("a", has_text="詳細").first
            log("DEBUG: 詳細リンク待機開始")
            detail_link.wait_for(state="visible", timeout=PLAYWRIGHT_TIMEOUT_MS)
            log("DEBUG: 詳細リンク待機完了")

            battle_href = detail_link.get_attribute("href")
            log(f"DEBUG: 取得したURL: {battle_href}")

            if not battle_href or battle_href == "/":
                log("Error: 詳細リンクの取得に失敗しました。")
                save_debug_files(page, "debug_search_error")
                return None

            game_id = extract_game_id_from_href(battle_href)
            log(f"Found Game ID: {game_id}")

            if not game_id:
                log("Error: Game ID の抽出に失敗しました。")
                save_debug_files(page, "debug_game_id_error")
                return None

            battle_url = urljoin(BASE_URL, battle_href)
            adapter_url = ADAPTER_URL.format(body=quote(battle_url, safe=""))
            log(f"DEBUG: adapter URL: {adapter_url}")

            t2 = time.perf_counter()
            page.goto(adapter_url, wait_until="domcontentloaded")
            log(f"DEBUG: adapter page goto完了 ({time.perf_counter() - t2:.2f}秒)")

            try:
                t3 = time.perf_counter()
                page.wait_for_load_state("networkidle", timeout=10_000)
                log(f"DEBUG: adapter page networkidle完了 ({time.perf_counter() - t3:.2f}秒)")
            except PlaywrightTimeoutError:
                log("DEBUG: adapter page networkidle タイムアウト")

            log("DEBUG: adapter page 現在URL確認")
            log(f"DEBUG: page.url = {page.url}")

            kif_text = ""
            for attempt in range(KIF_POLL_RETRY):
                log(f"DEBUG: KIF候補探索開始 {attempt + 1}/{KIF_POLL_RETRY}")
                candidates = collect_kif_candidates(page)
                if candidates:
                    kif_text = candidates[0]
                    log(f"DEBUG: KIF候補取得成功 len={len(kif_text)}")
                    break

                log(f"Waiting for KIF text... (attempt {attempt + 1}/{KIF_POLL_RETRY})")
                page.wait_for_timeout(KIF_POLL_INTERVAL_MS)

            if not kif_text:
                log("DEBUG: KIF候補が見つからないため、表示系ボタンを試します...")
                click_possible_reveal_buttons(page)

                for attempt in range(KIF_POLL_RETRY):
                    log(f"DEBUG: クリック後 KIF候補探索開始 {attempt + 1}/{KIF_POLL_RETRY}")
                    candidates = collect_kif_candidates(page)
                    if candidates:
                        kif_text = candidates[0]
                        log(f"DEBUG: KIF候補取得成功(クリック後) len={len(kif_text)}")
                        break

                    log(f"Retry after click... (attempt {attempt + 1}/{KIF_POLL_RETRY})")
                    page.wait_for_timeout(KIF_POLL_INTERVAL_MS)

            if not kif_text:
                log("Error: 棋譜データの取得に失敗しました。")
                save_debug_files(page, "debug_adapter_error")
                return None

            debug_kif_path = BASE_DIR / "debug_last_kif.kif"
            debug_kif_path.write_text(kif_text, encoding="utf-8")
            log(f"DEBUG: KIF保存: {debug_kif_path}")

            preview = "\n".join(kif_text.splitlines()[:8])
            log("DEBUG: KIF先頭プレビュー開始")
            print(preview, flush=True)
            log("DEBUG: KIF先頭プレビュー終了")

            opponent = extract_opponent_from_game_id(game_id, user_id)
            log(f"DEBUG: 相手名推定: {opponent}")
            log(f"DEBUG: scrape_with_playwright 完了 ({time.perf_counter() - start_all:.2f}秒)")
            return game_id, opponent, kif_text

        finally:
            browser.close()
            log("Playwright 終了")


# =========================
# KIF 解析
# =========================
def parse_moves_from_kif(kif_text: str):
    try:
        log("DEBUG: KIF.Parser.parse_str を試行")
        kif_obj = KIF.Parser.parse_str(kif_text)
        if hasattr(kif_obj, "moves"):
            log(f"DEBUG: parse_str 成功 moves={len(kif_obj.moves)}")
            return kif_obj.moves
    except Exception as e:
        log(f"DEBUG: KIF.Parser.parse_str 失敗: {e}")

    raise RuntimeError("KIFのパースに失敗しました。")


# =========================
# USI リスナー
# =========================
class Listener:
    def __init__(self):
        self.info = ""
        self.bestmove = ""

    def __call__(self, line):
        self.info = self.bestmove
        self.bestmove = line


# =========================
# 解析
# =========================
def analyze_and_save(user_id: str) -> None:
    total_start = time.perf_counter()
    log("analyze_and_save 開始")

    result = scrape_with_playwright(user_id)
    if not result:
        log("scrape_with_playwright が None を返しました。終了します。")
        return

    game_id, opponent, kif_text = result
    log(f"DEBUG: 取得した game_id={game_id}, opponent={opponent}, kif_len={len(kif_text)}")

    try:
        moves = parse_moves_from_kif(kif_text)
        log(f"DEBUG: KIF parse succeeded, moves={len(moves)}")
    except Exception as e:
        debug_kif_path = BASE_DIR / "debug_parse_error.kif"
        debug_kif_path.write_text(kif_text, encoding="utf-8")
        log(f"Error: KIFのパースに失敗しました: {e}")
        log(f"DEBUG: {debug_kif_path} を確認してください。")
        return

    black_name, white_name = extract_players_from_kif(kif_text)
    usi_moves = [cshogi.move_to_usi(move) for move in moves]

    board = cshogi.Board()
    sfen_list = [board.sfen()]
    scores = [0]

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
        t0 = time.perf_counter()
        engine = Engine(str(ENGINE_PATH), connect=True)
        log(f"DEBUG: Engine生成完了 ({time.perf_counter() - t0:.2f}秒)")

        engine.setoption("EvalDir", str(EVAL_DIR), listener=listener)
        engine.setoption("BookFile", "no_book", listener=listener)

        t1 = time.perf_counter()
        engine.isready(listener=listener)
        log(f"DEBUG: engine.isready 完了 ({time.perf_counter() - t1:.2f}秒)")

        engine.usinewgame(listener=listener)

        log(f"AI解析中... (全 {len(moves)} 手)")
        analysis_start = time.perf_counter()

        for i, move in enumerate(moves, start=1):
            move_start = time.perf_counter()

            board.push(move)
            sfen_list.append(board.sfen())

            try:
                listener.info = ""
                listener.bestmove = ""

                engine.position(sfen=board.sfen())
                engine.go(
                    byoyomi=ENGINE_TIME_MS,
                    listener=listener,
                )

                raw_score = usi_info_to_score(listener.info)
                if raw_score is None:
                    score = scores[-1]
                else:
                    score = int(raw_score)

            except Exception as e:
                log(f"Warning: {i} 手目の評価値取得に失敗しました: {e}")
                score = scores[-1]

            # 先手視点に揃える
            if board.turn == cshogi.WHITE:
                score = -score

            scores.append(score)

            if i <= 5 or i == len(moves) or i % 10 == 0:
                elapsed = time.perf_counter() - move_start
                total_elapsed = time.perf_counter() - analysis_start
                log(
                    f"DEBUG: {i}/{len(moves)} 手目 評価値={score} "
                    f"(この手 {elapsed:.2f}秒 / 累計 {total_elapsed:.2f}秒)"
                )

    finally:
        if engine is not None:
            try:
                engine.quit(listener=listener)
                log("DEBUG: engine.quit 完了")
            except Exception as e:
                log(f"DEBUG: engine.quit 失敗: {e}")

    diffs = [abs(scores[i] - scores[i - 1]) for i in range(1, len(scores))]
    worst_move_idx = int(np.argmax(diffs)) + 1 if diffs else 0
    log(f"DEBUG: worst_move_idx={worst_move_idx}")

    conn = sqlite3.connect(DB_NAME)
    conn.execute(
        """
        INSERT OR REPLACE INTO games
        (id, opponent, sfen_list, scores, worst_move, moves_usi, black_name, white_name)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            game_id,
            opponent,
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

    log(f"Successfully saved: {game_id}")
    log(f"analyze_and_save 完了 ({time.perf_counter() - total_start:.2f}秒)")


if __name__ == "__main__":
    init_db()
    analyze_and_save(DEFAULT_USER_ID)