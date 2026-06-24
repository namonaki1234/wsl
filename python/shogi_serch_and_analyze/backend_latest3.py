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
LATEST_GAME_COUNT = 3

ENGINE_TIME_MS_PER_MOVE = 500
ENGINE_TIME_MS_FOR_CANDIDATES = 1200
CANDIDATE_TOP_N = 3
CANDIDATE_PV_PLIES = 5


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

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS latest3_games (
            slot INTEGER PRIMARY KEY,
            game_id TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
        """
    )

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS game_candidates (
            game_id TEXT PRIMARY KEY,
            worst_prev_sfen TEXT,
            actual_move_usi TEXT,
            actual_move_jp TEXT,
            candidates_json TEXT,
            updated_at TEXT NOT NULL
        )
        """
    )

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

    candidate = "\n".join(lines[start_idx:]).strip()
    return candidate if is_kif_like(candidate) else ""


# =========================
# 日本語変換補助
# =========================
PIECE_JP = {
    "P": "歩",
    "L": "香",
    "N": "桂",
    "S": "銀",
    "G": "金",
    "B": "角",
    "R": "飛",
    "K": "玉",
    "+P": "と",
    "+L": "杏",
    "+N": "圭",
    "+S": "全",
    "+B": "馬",
    "+R": "龍",
}
HAND_JP = {
    "P": "歩",
    "L": "香",
    "N": "桂",
    "S": "銀",
    "G": "金",
    "B": "角",
    "R": "飛",
}
FULLWIDTH_NUM = {
    "1": "１", "2": "２", "3": "３", "4": "４", "5": "５",
    "6": "６", "7": "７", "8": "８", "9": "９",
}
RANK_KANJI = {
    "a": "一", "b": "二", "c": "三", "d": "四", "e": "五",
    "f": "六", "g": "七", "h": "八", "i": "九",
}
PROMOTABLE_JP = {"歩", "香", "桂", "銀", "角", "飛"}


def parse_sfen_board(sfen: str):
    board_part, turn, hand_part, move_no = sfen.split()

    grid = []
    for rank in board_part.split("/"):
        row = []
        i = 0
        while i < len(rank):
            ch = rank[i]
            if ch.isdigit():
                row.extend([None] * int(ch))
                i += 1
                continue

            promoted = False
            if ch == "+":
                promoted = True
                i += 1
                ch = rank[i]

            owner = "black" if ch.isupper() else "white"
            key = ("+" + ch.upper()) if promoted else ch.upper()
            row.append({"owner": owner, "piece": PIECE_JP[key]})
            i += 1

        grid.append(row)

    hands = {"black": {}, "white": {}}
    if hand_part != "-":
        i = 0
        while i < len(hand_part):
            count_str = ""
            while i < len(hand_part) and hand_part[i].isdigit():
                count_str += hand_part[i]
                i += 1

            piece_char = hand_part[i]
            i += 1
            count = int(count_str) if count_str else 1

            side = "black" if piece_char.isupper() else "white"
            code = piece_char.upper()
            hands[side][code] = hands[side].get(code, 0) + count

    return grid, turn, hands, int(move_no)


def square_to_jp(square: str) -> str:
    return FULLWIDTH_NUM[square[0]] + RANK_KANJI[square[1]]


def usi_move_to_jp(board: cshogi.Board, usi: str) -> str:
    sfen = board.sfen()
    grid, _, _, _ = parse_sfen_board(sfen)

    m_drop = re.match(r"^([PLNSGBRK])\*([1-9][a-i])$", usi)
    if m_drop:
        piece_code = m_drop.group(1)
        dest = m_drop.group(2)
        piece_jp = HAND_JP.get(piece_code, piece_code)
        return f"{square_to_jp(dest)}{piece_jp}打"

    m_move = re.match(r"^([1-9][a-i])([1-9][a-i])(\+?)$", usi)
    if not m_move:
        return usi

    src = m_move.group(1)
    dest = m_move.group(2)
    promote = m_move.group(3) == "+"

    src_file = int(src[0])
    src_rank = ord(src[1]) - ord("a")
    src_row = src_rank
    src_col = 9 - src_file

    piece_name = "?"
    if 0 <= src_row < 9 and 0 <= src_col < 9:
        cell = grid[src_row][src_col]
        if cell is not None:
            piece_name = cell["piece"]

    dest_jp = square_to_jp(dest)

    if promote and piece_name in PROMOTABLE_JP:
        return f"{dest_jp}{piece_name}成"

    return f"{dest_jp}{piece_name}"


def pv_to_japanese(prev_sfen: str, pv_usi: list[str]) -> list[str]:
    board = cshogi.Board(prev_sfen)
    jp_moves = []

    for mv in pv_usi:
        try:
            jp = usi_move_to_jp(board, mv)
        except Exception:
            jp = mv
        jp_moves.append(jp)

        try:
            board.push_usi(mv)
        except Exception:
            break

    return jp_moves


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

            for _, raw in payload.items():
                candidate = trim_to_kif(raw)
                if candidate and candidate not in seen:
                    seen.add(candidate)
                    candidates.append(candidate)

    candidates.sort(key=len, reverse=True)
    return candidates


def scrape_latest_game_kifs(user_id: str, limit: int = 3):
    results = []

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
            search_url = SEARCH_URL.format(query=quote(user_id, safe=""))
            log(f"Searching latest {limit} games for {user_id}")
            log(f"DEBUG: search URL: {search_url}")

            page.goto(search_url, wait_until="domcontentloaded")

            try:
                page.wait_for_load_state("networkidle", timeout=10_000)
            except PlaywrightTimeoutError:
                log("DEBUG: search page networkidle タイムアウト")

            detail_links = page.locator("a", has_text="詳細")
            count = detail_links.count()

            hrefs = []
            seen = set()
            for i in range(count):
                href = detail_links.nth(i).get_attribute("href")
                if not href or href in seen:
                    continue
                seen.add(href)
                hrefs.append(href)
                if len(hrefs) >= limit:
                    break

            if not hrefs:
                log("Error: 詳細リンクが取得できませんでした。")
                save_debug_files(page, "debug_latest3_search_error")
                return []

            log(f"DEBUG: 取得できた詳細リンク数 = {len(hrefs)}")

            for idx, href in enumerate(hrefs, start=1):
                game_id = extract_game_id_from_href(href)
                if not game_id:
                    log(f"Warning: game_id 抽出失敗 href={href}")
                    continue

                battle_url = urljoin(BASE_URL, href)
                adapter_url = ADAPTER_URL.format(body=quote(battle_url, safe=""))
                log(f"[{idx}/{len(hrefs)}] game_id={game_id}")
                log(f"DEBUG: adapter URL: {adapter_url}")

                page.goto(adapter_url, wait_until="domcontentloaded")
                try:
                    page.wait_for_load_state("networkidle", timeout=10_000)
                except PlaywrightTimeoutError:
                    log("DEBUG: adapter page networkidle タイムアウト")

                kif_text = ""
                for attempt in range(KIF_POLL_RETRY):
                    candidates = collect_kif_candidates(page)
                    if candidates:
                        kif_text = candidates[0]
                        break
                    page.wait_for_timeout(KIF_POLL_INTERVAL_MS)

                if not kif_text:
                    log(f"Warning: KIF取得失敗 game_id={game_id}")
                    save_debug_files(page, f"debug_kif_error_{idx}")
                    continue

                opponent = extract_opponent_from_game_id(game_id, user_id)
                results.append({
                    "slot": idx,
                    "game_id": game_id,
                    "opponent": opponent,
                    "kif_text": kif_text,
                })

        finally:
            browser.close()
            log("Playwright 終了")

    return results


# =========================
# KIF解析
# =========================
def parse_moves_from_kif(kif_text: str):
    try:
        kif_obj = KIF.Parser.parse_str(kif_text)
        if hasattr(kif_obj, "moves"):
            return kif_obj.moves
    except Exception as e:
        raise RuntimeError(f"KIFのパースに失敗しました: {e}") from e

    raise RuntimeError("KIFのパースに失敗しました。")


# =========================
# エンジン用リスナー
# =========================
class ScoreListener:
    def __init__(self):
        self.last_info = ""

    def __call__(self, line: str):
        if line.startswith("info "):
            self.last_info = line


class PVListener:
    def __init__(self):
        self.lines = []

    def __call__(self, line: str):
        if line.startswith("info ") and " pv " in line:
            self.lines.append(line)


def parse_info_lines(lines, pv_plies=5, top_n=3):
    best = {}

    for line in lines:
        mp_match = re.search(r"\bmultipv\s+(\d+)", line)
        multipv = int(mp_match.group(1)) if mp_match else 1

        depth_match = re.search(r"\bdepth\s+(\d+)", line)
        depth = int(depth_match.group(1)) if depth_match else 0

        score_type = None
        score_value = None

        cp_match = re.search(r"\bscore\s+cp\s+(-?\d+)", line)
        mate_match = re.search(r"\bscore\s+mate\s+(-?\d+)", line)

        if cp_match:
            score_type = "cp"
            score_value = int(cp_match.group(1))
        elif mate_match:
            score_type = "mate"
            score_value = int(mate_match.group(1))

        if " pv " not in line:
            continue

        pv = line.split(" pv ", 1)[1].strip().split()

        if multipv not in best or depth >= best[multipv]["depth"]:
            best[multipv] = {
                "multipv": multipv,
                "depth": depth,
                "score_type": score_type,
                "score_value": score_value,
                "pv": pv,
            }

    results = []
    for mp in sorted(best.keys())[:top_n]:
        entry = best[mp]
        pv5 = entry["pv"][:pv_plies]
        results.append(
            {
                "multipv": entry["multipv"],
                "depth": entry["depth"],
                "score_type": entry["score_type"],
                "score_value": entry["score_value"],
                "pv5": pv5,
            }
        )
    return results


# =========================
# 保存
# =========================
def save_game_and_candidates(
    slot: int,
    game_id: str,
    opponent: str,
    sfen_list: list[str],
    scores: list[int],
    worst_move_idx: int,
    moves_usi: list[str],
    black_name: str,
    white_name: str,
    worst_prev_sfen: str,
    actual_move_usi: str,
    actual_move_jp: str,
    candidates: list[dict],
) -> None:
    updated_at = datetime.now().isoformat(timespec="seconds")
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
            json.dumps(moves_usi, ensure_ascii=False),
            black_name,
            white_name,
        ),
    )

    conn.execute(
        """
        INSERT OR REPLACE INTO game_candidates
        (game_id, worst_prev_sfen, actual_move_usi, actual_move_jp, candidates_json, updated_at)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        (
            game_id,
            worst_prev_sfen,
            actual_move_usi,
            actual_move_jp,
            json.dumps(candidates, ensure_ascii=False),
            updated_at,
        ),
    )

    conn.execute(
        """
        INSERT OR REPLACE INTO latest3_games (slot, game_id, updated_at)
        VALUES (?, ?, ?)
        """,
        (slot, game_id, updated_at),
    )

    conn.commit()
    conn.close()


def cleanup_latest3_unused_slots(valid_slots: int) -> None:
    conn = sqlite3.connect(DB_NAME)
    conn.execute("DELETE FROM latest3_games WHERE slot > ?", (valid_slots,))
    conn.commit()
    conn.close()


# =========================
# 解析本体
# =========================
def analyze_single_game(item: dict) -> None:
    slot = item["slot"]
    game_id = item["game_id"]
    opponent = item["opponent"]
    kif_text = item["kif_text"]

    log(f"解析開始: slot={slot}, game_id={game_id}")

    moves = parse_moves_from_kif(kif_text)
    black_name, white_name = extract_players_from_kif(kif_text)
    usi_moves = [cshogi.move_to_usi(mv) for mv in moves]

    board = cshogi.Board()
    sfen_list = [board.sfen()]
    scores = [0]

    if not ENGINE_PATH.exists():
        raise FileNotFoundError(f"エンジンが見つかりません: {ENGINE_PATH}")
    if not EVAL_DIR.exists():
        raise FileNotFoundError(f"EvalDir が見つかりません: {EVAL_DIR}")

    engine = Engine(str(ENGINE_PATH), connect=True)
    score_listener = ScoreListener()

    try:
        engine.setoption("EvalDir", str(EVAL_DIR), listener=score_listener)
        engine.setoption("BookFile", "no_book", listener=score_listener)
        engine.isready(listener=score_listener)
        engine.usinewgame(listener=score_listener)

        log(f"AI解析中... ({game_id}, 全 {len(moves)} 手)")
        for i, move in enumerate(moves, start=1):
            board.push(move)
            sfen_list.append(board.sfen())

            try:
                score_listener.last_info = ""
                engine.position(sfen=board.sfen())
                engine.go(byoyomi=ENGINE_TIME_MS_PER_MOVE, listener=score_listener)

                raw_score = usi_info_to_score(score_listener.last_info)
                score = int(raw_score) if raw_score is not None else scores[-1]
            except Exception as e:
                log(f"Warning: {game_id} {i}手目の評価値取得失敗: {e}")
                score = scores[-1]

            if board.turn == cshogi.WHITE:
                score = -score

            scores.append(score)

            if i <= 3 or i == len(moves) or i % 20 == 0:
                log(f"DEBUG: {game_id} {i}/{len(moves)} 手目 評価値={score}")

        diffs = [abs(scores[i] - scores[i - 1]) for i in range(1, len(scores))]
        worst_move_idx = int(np.argmax(diffs)) + 1 if diffs else 0

        if worst_move_idx <= 0 or worst_move_idx > len(moves):
            worst_move_idx = 1

        worst_prev_sfen = sfen_list[worst_move_idx - 1]
        actual_move_usi = usi_moves[worst_move_idx - 1]
        actual_move_jp = pv_to_japanese(worst_prev_sfen, [actual_move_usi])[0]

        pv_listener = PVListener()
        engine.setoption("MultiPV", CANDIDATE_TOP_N, listener=pv_listener)
        engine.isready(listener=pv_listener)
        engine.usinewgame(listener=pv_listener)
        engine.position(sfen=worst_prev_sfen)
        engine.go(byoyomi=ENGINE_TIME_MS_FOR_CANDIDATES, listener=pv_listener)

        raw_candidates = parse_info_lines(
            pv_listener.lines,
            pv_plies=CANDIDATE_PV_PLIES,
            top_n=CANDIDATE_TOP_N,
        )

        candidates = []
        for cand in raw_candidates:
            end_board = cshogi.Board(worst_prev_sfen)
            pv5 = cand["pv5"]
            pv5_jp = pv_to_japanese(worst_prev_sfen, pv5)

            for mv in pv5:
                try:
                    end_board.push_usi(mv)
                except Exception:
                    break

            candidates.append(
                {
                    "multipv": cand["multipv"],
                    "depth": cand["depth"],
                    "score_type": cand["score_type"],
                    "score_value": cand["score_value"],
                    "pv5": pv5,
                    "pv5_jp": pv5_jp,
                    "end_sfen": end_board.sfen(),
                }
            )

        save_game_and_candidates(
            slot=slot,
            game_id=game_id,
            opponent=opponent,
            sfen_list=sfen_list,
            scores=scores,
            worst_move_idx=worst_move_idx,
            moves_usi=usi_moves,
            black_name=black_name,
            white_name=white_name,
            worst_prev_sfen=worst_prev_sfen,
            actual_move_usi=actual_move_usi,
            actual_move_jp=actual_move_jp,
            candidates=candidates,
        )

        log(f"保存完了: {game_id}")

    finally:
        try:
            engine.quit(listener=score_listener)
        except Exception:
            pass


def analyze_latest3(user_id: str) -> None:
    init_db()

    items = scrape_latest_game_kifs(user_id, limit=LATEST_GAME_COUNT)
    if not items:
        log("最新対局を取得できませんでした。")
        return

    for item in items:
        try:
            analyze_single_game(item)
        except Exception as e:
            log(f"Error: {item['game_id']} の解析に失敗しました: {e}")

    cleanup_latest3_unused_slots(len(items))
    log("latest3 解析完了")


if __name__ == "__main__":
    analyze_latest3(DEFAULT_USER_ID)