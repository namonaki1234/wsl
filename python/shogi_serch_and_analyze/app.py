import html
import json
import re
import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib as mpl
import cshogi
import streamlit as st
from cshogi.usi import Engine


# =========================
# 設定
# =========================
BASE_DIR = Path(__file__).resolve().parent
DB_NAME = BASE_DIR / "shogi_wars.db"

ENGINE_PATH = BASE_DIR / "YaneuraOu" / "source" / "YaneuraOu-by-gcc"
EVAL_DIR = BASE_DIR / "YaneuraOu" / "source" / "eval"

st.set_page_config(layout="wide", page_title="将棋解析GUI")


# =========================
# 日本語フォント設定
# =========================
import matplotlib as mpl
import matplotlib.font_manager as fm

def configure_matplotlib_japanese():
    candidates = [
        "Noto Sans CJK JP",
        "Noto Serif CJK JP",
        "IPAexGothic",
        "IPAGothic",
        "Yu Gothic",
        "Meiryo",
    ]

    installed = {f.name for f in fm.fontManager.ttflist}
    chosen = None
    for name in candidates:
        if name in installed:
            chosen = name
            break

    if chosen is not None:
        mpl.rcParams["font.family"] = "sans-serif"
        mpl.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
        mpl.rcParams["axes.unicode_minus"] = False
        return chosen

    mpl.rcParams["axes.unicode_minus"] = False
    return None


JP_FONT_NAME = configure_matplotlib_japanese()


# =========================
# CSS
# =========================
st.markdown(
    """
    <style>
    .shogi-board-shell {
        display: inline-flex;
        flex-direction: column;
        align-items: center;
        gap: 10px;
        background: #f3e0b5;
        padding: 18px 20px;
        border-radius: 16px;
        border: 2px solid #c9a86a;
        box-shadow: 0 6px 18px rgba(0,0,0,0.08);
        margin: 0 auto;
    }

    .shogi-player-label {
        font-weight: 700;
        color: #5a3818;
        font-size: 22px;
    }

    .shogi-handline {
        color: #6b4a25;
        font-size: 15px;
    }

    .shogi-board-table {
        border-collapse: collapse;
        background: #e1bc75;
        margin: 0 auto;
    }

    .shogi-file-label, .shogi-rank-label {
        color: #5a3818;
        font-weight: 700;
        text-align: center;
        font-size: 16px;
        background: transparent;
        border: none;
        padding: 2px 4px;
    }

    .shogi-cell {
        border: 1px solid #5c4525;
        text-align: center;
        vertical-align: middle;
        background: #e7c98f;
        position: relative;
    }

    .shogi-cell.lastmove {
        background: #ffd76a !important;
        box-shadow: inset 0 0 0 3px rgba(210, 72, 35, 0.55);
    }

    .shogi-piece {
        font-weight: 700;
        color: #222;
        line-height: 1;
        display: inline-block;
        user-select: none;
    }

    .shogi-piece.white {
        transform: rotate(180deg);
    }

    .shogi-piece.black-flipped {
        transform: rotate(180deg);
    }

    .pv-card {
        border: 1px solid #ddd;
        border-radius: 12px;
        padding: 14px;
        background: #fafafa;
    }

    .pv-line {
        margin-top: 8px;
        padding: 8px 10px;
        border-radius: 8px;
        background: #f4f4f4;
        line-height: 1.7;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# =========================
# 定数
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
HAND_ORDER = ["R", "B", "G", "S", "N", "L", "P"]
HAND_JP = {
    "P": "歩",
    "L": "香",
    "N": "桂",
    "S": "銀",
    "G": "金",
    "B": "角",
    "R": "飛",
}
RANK_LABELS = ["一", "二", "三", "四", "五", "六", "七", "八", "九"]

FULLWIDTH_NUM = {
    "1": "１", "2": "２", "3": "３", "4": "４", "5": "５",
    "6": "６", "7": "７", "8": "８", "9": "９",
}
RANK_KANJI = {
    "a": "一", "b": "二", "c": "三", "d": "四", "e": "五",
    "f": "六", "g": "七", "h": "八", "i": "九",
}
PROMOTABLE_JP = {"歩", "香", "桂", "銀", "角", "飛"}


# =========================
# DB
# =========================
def get_connection():
    return sqlite3.connect(DB_NAME)


def get_table_columns():
    conn = get_connection()
    cols = [row[1] for row in conn.execute("PRAGMA table_info(games)").fetchall()]
    conn.close()
    return cols


def get_game_rows():
    if not DB_NAME.exists():
        return []

    conn = get_connection()
    rows = conn.execute(
        "SELECT id, opponent FROM games ORDER BY rowid DESC"
    ).fetchall()
    conn.close()
    return rows


def load_game(game_id: str):
    cols = get_table_columns()

    select_cols = ["id", "opponent", "sfen_list", "scores", "worst_move"]
    for col in ["moves_usi", "black_name", "white_name"]:
        if col in cols:
            select_cols.append(col)

    conn = get_connection()
    row = conn.execute(
        f"SELECT {', '.join(select_cols)} FROM games WHERE id = ?",
        (game_id,),
    ).fetchone()
    conn.close()

    return row, cols


# =========================
# 補助
# =========================
def format_game_label(row):
    game_id, opponent = row
    return f"{game_id}  vs {opponent}"


def clamp_step(step: int, max_step: int) -> int:
    return max(0, min(step, max_step))


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


def format_hands_text(hands: dict, side: str) -> str:
    items = []
    for code in HAND_ORDER:
        cnt = hands[side].get(code, 0)
        if cnt > 0:
            jp = HAND_JP[code]
            items.append(f"{jp}{cnt}" if cnt > 1 else jp)
    return "なし" if not items else " ".join(items)


def usi_to_dest_square(usi: str):
    if not usi or usi in ("開始局面", "(USI未保存)"):
        return None

    m_drop = re.match(r"^[PLNSGBRK]\*([1-9][a-i])$", usi)
    if m_drop:
        sq = m_drop.group(1)
    else:
        m_move = re.match(r"^[1-9][a-i]([1-9][a-i])\+?$", usi)
        if not m_move:
            return None
        sq = m_move.group(1)

    file_num = int(sq[0])
    rank_idx = ord(sq[1]) - ord("a")
    col_idx = 9 - file_num
    row_idx = rank_idx
    return row_idx, col_idx


def infer_names_from_game_id(game_id: str):
    parts = game_id.split("-")
    if len(parts) >= 3:
        return parts[0], parts[1]
    return "先手", "後手"


def piece_class(owner: str, flipped: bool) -> str:
    base = "shogi-piece"
    if not flipped:
        if owner == "white":
            return base + " white"
        return base
    else:
        if owner == "black":
            return base + " black-flipped"
        return base


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


def render_board_html(
    sfen: str,
    black_name: str,
    white_name: str,
    lastmove_usi: str | None = None,
    cell_px: int = 64,
    show_names: bool = True,
    flipped: bool = False,
):
    grid, turn, hands, _ = parse_sfen_board(sfen)
    highlight = usi_to_dest_square(lastmove_usi)

    row_order = list(range(9)) if not flipped else list(range(8, -1, -1))
    col_order = list(range(9)) if not flipped else list(range(8, -1, -1))

    file_labels = list(range(9, 0, -1)) if not flipped else list(range(1, 10))
    rank_labels = [RANK_LABELS[r] for r in row_order]

    top_name = white_name if not flipped else black_name
    bottom_name = black_name if not flipped else white_name
    top_side = "white" if not flipped else "black"
    bottom_side = "black" if not flipped else "white"

    html_parts = []
    html_parts.append('<div class="shogi-board-shell">')

    if show_names:
        html_parts.append(
            f'<div class="shogi-player-label">{html.escape("後手" if not flipped else "先手")}：{html.escape(top_name)}</div>'
        )
        html_parts.append(
            f'<div class="shogi-handline">持駒：{html.escape(format_hands_text(hands, top_side))}</div>'
        )

    html_parts.append('<table class="shogi-board-table">')

    html_parts.append("<tr><th></th>")
    for file_label in file_labels:
        html_parts.append(f'<th class="shogi-file-label">{file_label}</th>')
    html_parts.append("<th></th></tr>")

    for view_r, real_r in enumerate(row_order):
        html_parts.append(f'<tr><th class="shogi-rank-label">{rank_labels[view_r]}</th>')

        for real_c in col_order:
            cell = grid[real_r][real_c]
            classes = ["shogi-cell"]
            if highlight == (real_r, real_c):
                classes.append("lastmove")

            style = f'width:{cell_px}px; height:{cell_px}px;'
            html_parts.append(f'<td class="{" ".join(classes)}" style="{style}">')

            if cell is not None:
                owner = cell["owner"]
                piece = cell["piece"]
                cls = piece_class(owner, flipped)
                font_px = int(cell_px * (0.54 if len(piece) == 1 else 0.40))
                html_parts.append(
                    f'<span class="{cls}" style="font-size:{font_px}px;">{html.escape(piece)}</span>'
                )

            html_parts.append("</td>")

        html_parts.append(f'<th class="shogi-rank-label">{rank_labels[view_r]}</th></tr>')

    html_parts.append("<tr><th></th>")
    for file_label in file_labels:
        html_parts.append(f'<th class="shogi-file-label">{file_label}</th>')
    html_parts.append("<th></th></tr>")

    html_parts.append("</table>")

    if show_names:
        html_parts.append(
            f'<div class="shogi-handline">持駒：{html.escape(format_hands_text(hands, bottom_side))}</div>'
        )
        html_parts.append(
            f'<div class="shogi-player-label">{html.escape("先手" if not flipped else "後手")}：{html.escape(bottom_name)}</div>'
        )

    html_parts.append("</div>")
    return "".join(html_parts)


# =========================
# エンジン解析
# =========================
class PVListener:
    def __init__(self):
        self.lines = []

    def __call__(self, line):
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


@st.cache_data(show_spinner=False)
def analyze_best_lines_cached(prev_sfen: str, top_n: int = 3, pv_plies: int = 5, think_ms: int = 1200):
    listener = PVListener()
    engine = None

    try:
        engine = Engine(str(ENGINE_PATH), connect=True)
        engine.setoption("EvalDir", str(EVAL_DIR), listener=listener)
        engine.setoption("BookFile", "no_book", listener=listener)
        engine.setoption("MultiPV", top_n, listener=listener)
        engine.isready(listener=listener)
        engine.usinewgame(listener=listener)
        engine.position(sfen=prev_sfen)
        engine.go(byoyomi=think_ms, listener=listener)

        items = parse_info_lines(listener.lines, pv_plies=pv_plies, top_n=top_n)

        enriched = []
        for item in items:
            board = cshogi.Board(prev_sfen)
            jp_line = pv_to_japanese(prev_sfen, item["pv5"])

            for mv in item["pv5"]:
                try:
                    board.push_usi(mv)
                except Exception:
                    break

            enriched.append(
                {
                    **item,
                    "pv5_jp": jp_line,
                    "end_sfen": board.sfen(),
                }
            )
        return enriched

    finally:
        if engine is not None:
            try:
                engine.quit(listener=listener)
            except Exception:
                pass


def format_score_text(score_type, score_value):
    if score_type == "mate":
        if score_value is None:
            return "詰み評価"
        return f"Mate {score_value}"
    if score_type == "cp":
        if score_value is None:
            return "不明"
        sign = "+" if score_value > 0 else ""
        return f"{sign}{score_value}"
    return "不明"


# =========================
# コールバック
# =========================
def on_slider_change():
    st.session_state.current_step = clamp_step(
        int(st.session_state.step_slider),
        int(st.session_state.max_step),
    )


def go_prev():
    st.session_state.current_step = clamp_step(
        int(st.session_state.current_step) - 1,
        int(st.session_state.max_step),
    )
    st.session_state.step_slider = st.session_state.current_step


def go_next():
    st.session_state.current_step = clamp_step(
        int(st.session_state.current_step) + 1,
        int(st.session_state.max_step),
    )
    st.session_state.step_slider = st.session_state.current_step


# =========================
# サイドバー
# =========================
st.sidebar.title("将棋ウォーズ解析")

game_rows = get_game_rows()
if not game_rows:
    st.info("まだ対局データがありません。先に backend.py を実行してください。")
    st.stop()

game_id_list = [row[0] for row in game_rows]
game_label_map = {row[0]: format_game_label(row) for row in game_rows}

selected_id = st.sidebar.selectbox(
    "対局を選択",
    game_id_list,
    format_func=lambda game_id: game_label_map[game_id],
)

is_flipped = st.sidebar.checkbox("盤面を上下反転", value=False)

if st.session_state.get("selected_game_id") != selected_id:
    st.session_state.selected_game_id = selected_id
    st.session_state.current_step = 0
    st.session_state.step_slider = 0
    st.session_state.branch_move = ""


# =========================
# 読み込み
# =========================
row, cols = load_game(selected_id)
if row is None:
    st.error("対局データの読み込みに失敗しました。")
    st.stop()

idx = 0
game_id = row[idx]; idx += 1
opponent = row[idx]; idx += 1
sfen_json = row[idx]; idx += 1
scores_json = row[idx]; idx += 1
worst_move = row[idx]; idx += 1

moves_usi = []
black_name = ""
white_name = ""

if "moves_usi" in cols:
    moves_json = row[idx]; idx += 1
    moves_usi = json.loads(moves_json) if moves_json else []

if "black_name" in cols:
    black_name = row[idx] or ""
    idx += 1

if "white_name" in cols:
    white_name = row[idx] or ""
    idx += 1

if not black_name or not white_name:
    inferred_black, inferred_white = infer_names_from_game_id(game_id)
    black_name = black_name or inferred_black
    white_name = white_name or inferred_white

sfen_list = json.loads(sfen_json)
scores = json.loads(scores_json)

if not sfen_list or not scores:
    st.error("局面データまたは評価値データが空です。")
    st.stop()

max_step = len(sfen_list) - 1
st.session_state.max_step = max_step

if "current_step" not in st.session_state:
    st.session_state.current_step = 0
if "step_slider" not in st.session_state:
    st.session_state.step_slider = 0
if "branch_move" not in st.session_state:
    st.session_state.branch_move = ""

st.session_state.current_step = clamp_step(st.session_state.current_step, max_step)
st.session_state.step_slider = clamp_step(st.session_state.step_slider, max_step)

step = st.session_state.current_step

if step == 0:
    current_move_usi = "開始局面"
    current_move_jp = "開始局面"
else:
    if 0 <= step - 1 < len(moves_usi):
        current_move_usi = moves_usi[step - 1]
        prev_sfen_for_current = sfen_list[step - 1]
        try:
            current_move_jp = pv_to_japanese(prev_sfen_for_current, [current_move_usi])[0]
        except Exception:
            current_move_jp = current_move_usi
    else:
        current_move_usi = "(USI未保存)"
        current_move_jp = "(棋譜未保存)"

highlight_move_usi = None if step == 0 else current_move_usi


# =========================
# メイン表示
# =========================
st.title("将棋ウォーズ解析GUI")
st.caption(f"対局ID: {game_id} / 相手: {opponent}")

if JP_FONT_NAME is None:
    st.info("グラフの日本語表示が崩れる場合は、WSL 側で `sudo apt-get install -y fonts-noto-cjk` を実行すると改善しやすいです。")

st.subheader("評価値グラフ")
fig, ax = plt.subplots(figsize=(12, 3))
ax.plot(scores, lw=2)
ax.axhline(0, lw=0.8)

if 0 <= worst_move < len(scores):
    ax.scatter([worst_move], [scores[worst_move]], zorder=5, label="Worst Move")
    ax.legend()

ax.set_xlabel("手数")
ax.set_ylabel("評価値")
ax.grid(alpha=0.3)
st.pyplot(fig)

col1, col2 = st.columns([1.9, 1])

with col1:
    st.subheader("盤面")

    st.slider(
        "手数",
        min_value=0,
        max_value=max_step,
        key="step_slider",
        on_change=on_slider_change,
    )

    step = st.session_state.current_step
    current_sfen = sfen_list[step]

    st.markdown(f"**現在の指し手（USI）:** `{current_move_usi}`")
    st.markdown(f"**現在の指し手（日本語）:** `{current_move_jp}`")

    st.text_input(
        "別の手を検討（USI形式，例: 7g7f）",
        key="branch_move",
    )

    branch_move = st.session_state.branch_move.strip()
    board_sfen_to_show = current_sfen
    move_to_highlight = highlight_move_usi

    if branch_move:
        board = cshogi.Board(current_sfen)
        try:
            board.push_usi(branch_move)
            board_sfen_to_show = board.sfen()
            move_to_highlight = branch_move
            st.warning(f"分岐検討中: {branch_move}")
        except Exception:
            st.error("正しいUSI形式で入力してください。例: 7g7f")

    st.markdown(
        render_board_html(
            board_sfen_to_show,
            black_name=black_name,
            white_name=white_name,
            lastmove_usi=move_to_highlight,
            cell_px=64,
            show_names=True,
            flipped=is_flipped,
        ),
        unsafe_allow_html=True,
    )

    with st.expander("局面のSFENを表示"):
        st.code(current_sfen, language=None)

with col2:
    step = st.session_state.current_step

    st.subheader("情報")
    st.metric("現在の手数", f"{step} 手目")

    score_value = scores[step] if 0 <= step < len(scores) else 0
    st.metric("評価値", f"{score_value}")

    st.markdown(f"**現在の指し手（USI）:** `{current_move_usi}`")
    st.markdown(f"**現在の指し手（日本語）:** `{current_move_jp}`")
    st.markdown(f"**先手:** {black_name}")
    st.markdown(f"**後手:** {white_name}")
    st.markdown(f"**表示向き:** {'反転中（後手視点）' if is_flipped else '通常（先手視点）'}")

    if step == worst_move:
        st.error("🚨 ここが最大の悪手候補です")

    st.write(f"最大悪手候補: {worst_move} 手目")

    btn_col1, btn_col2 = st.columns(2)
    with btn_col1:
        st.button("◀ 前へ", use_container_width=True, on_click=go_prev)
    with btn_col2:
        st.button("次へ ▶", use_container_width=True, on_click=go_next)


# =========================
# 悪手の代替候補
# =========================
if step == worst_move and step > 0:
    st.subheader("悪手の代替候補（3候補・各5手読み）")

    prev_sfen = sfen_list[step - 1]
    actual_move = moves_usi[step - 1] if 0 <= step - 1 < len(moves_usi) else "(USI未保存)"
    try:
        actual_move_jp = pv_to_japanese(prev_sfen, [actual_move])[0]
    except Exception:
        actual_move_jp = actual_move

    st.markdown(f"**実際に指した手（USI）:** `{actual_move}`")
    st.markdown(f"**実際に指した手（日本語）:** `{actual_move_jp}`")

    with st.spinner("最善手候補を解析中..."):
        candidates = analyze_best_lines_cached(
            prev_sfen=prev_sfen,
            top_n=3,
            pv_plies=5,
            think_ms=1200,
        )

    if candidates:
        tabs = st.tabs([f"候補 {c['multipv']}" for c in candidates])

        for tab, cand in zip(tabs, candidates):
            with tab:
                first_move = cand["pv5"][0] if cand["pv5"] else "(なし)"
                first_move_jp = cand["pv5_jp"][0] if cand["pv5_jp"] else "(なし)"
                seq_text_usi = " → ".join(cand["pv5"]) if cand["pv5"] else "(読み筋なし)"
                seq_text_jp = " → ".join(cand["pv5_jp"]) if cand["pv5_jp"] else "(読み筋なし)"
                score_text = format_score_text(cand["score_type"], cand["score_value"])

                st.markdown(
                    f"""
                    <div class="pv-card">
                        <div><strong>候補手（USI）:</strong> <code>{html.escape(first_move)}</code></div>
                        <div><strong>候補手（日本語）:</strong> <code>{html.escape(first_move_jp)}</code></div>
                        <div><strong>評価値:</strong> {html.escape(score_text)}</div>
                        <div class="pv-line"><strong>5手の読み筋（USI）:</strong><br><code>{html.escape(seq_text_usi)}</code></div>
                        <div class="pv-line"><strong>5手の読み筋（日本語）:</strong><br><code>{html.escape(seq_text_jp)}</code></div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

                st.markdown("#### 5手後の局面")
                st.markdown(
                    render_board_html(
                        cand["end_sfen"],
                        black_name=black_name,
                        white_name=white_name,
                        lastmove_usi=cand["pv5"][-1] if cand["pv5"] else None,
                        cell_px=40,
                        show_names=False,
                        flipped=is_flipped,
                    ),
                    unsafe_allow_html=True,
                )
    else:
        st.info("代替候補を取得できませんでした。")


if "black_name" not in cols or "white_name" not in cols:
    st.info("先手・後手のユーザー名を正確に表示するには、backend.py を更新して再解析してください。")