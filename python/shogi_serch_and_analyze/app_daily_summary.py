import html
import json
import re
import sqlite3
from pathlib import Path

import streamlit as st


# =========================
# 設定
# =========================
BASE_DIR = Path(__file__).resolve().parent
DB_NAME = BASE_DIR / "shogi_wars.db"

st.set_page_config(layout="wide", page_title="将棋ウォーズ 3局まとめ")


# =========================
# CSS
# =========================
st.markdown(
    """
    <style>
    .game-card {
        border: 1px solid #ddd;
        border-radius: 16px;
        padding: 18px;
        margin-bottom: 22px;
        background: #fff;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
    }

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
        font-size: 20px;
    }

    .shogi-handline {
        color: #6b4a25;
        font-size: 14px;
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
        font-size: 15px;
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
        padding: 12px;
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


# =========================
# DB
# =========================
def get_connection():
    return sqlite3.connect(DB_NAME)


def load_latest3_ids():
    if not DB_NAME.exists():
        return []

    conn = get_connection()
    rows = conn.execute(
        """
        SELECT slot, game_id, updated_at
        FROM latest3_games
        ORDER BY slot
        """
    ).fetchall()
    conn.close()
    return rows


def load_game(game_id: str):
    conn = get_connection()
    row = conn.execute(
        """
        SELECT id, opponent, sfen_list, scores, worst_move, moves_usi, black_name, white_name
        FROM games
        WHERE id = ?
        """,
        (game_id,),
    ).fetchone()
    conn.close()
    return row


def load_candidate_summary(game_id: str):
    conn = get_connection()
    row = conn.execute(
        """
        SELECT worst_prev_sfen, actual_move_usi, actual_move_jp, candidates_json, updated_at
        FROM game_candidates
        WHERE game_id = ?
        """,
        (game_id,),
    ).fetchone()
    conn.close()
    return row


# =========================
# 補助
# =========================
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


def piece_class(owner: str, flipped: bool) -> str:
    base = "shogi-piece"
    if not flipped:
        return base + " white" if owner == "white" else base
    else:
        return base + " black-flipped" if owner == "black" else base


def render_board_html(
    sfen: str,
    black_name: str,
    white_name: str,
    lastmove_usi: str | None = None,
    cell_px: int = 56,
    show_names: bool = True,
    flipped: bool = False,
):
    grid, _, hands, _ = parse_sfen_board(sfen)
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
# 画面
# =========================
st.title("将棋ウォーズ 最新3局まとめ")
st.caption("各対局の最大悪手局面と，候補手3つ・各5手読みを一覧表示します。")

is_flipped = st.sidebar.checkbox("盤面を上下反転", value=False)

rows = load_latest3_ids()
if not rows:
    st.warning("まだ latest3 データがありません。先に `uv run python backend_latest3.py` を実行してください。")
    st.stop()

latest_updated = max(row[2] for row in rows)
st.info(f"最新更新: {latest_updated}")

for slot, game_id, updated_at in rows:
    game_row = load_game(game_id)
    cand_row = load_candidate_summary(game_id)

    if game_row is None or cand_row is None:
        st.error(f"game_id={game_id} のデータが不足しています。backend_latest3.py を再実行してください。")
        continue

    (
        _game_id,
        opponent,
        sfen_json,
        scores_json,
        worst_move,
        moves_json,
        black_name,
        white_name,
    ) = game_row

    (
        worst_prev_sfen,
        actual_move_usi,
        actual_move_jp,
        candidates_json,
        _cand_updated_at,
    ) = cand_row

    sfen_list = json.loads(sfen_json)
    scores = json.loads(scores_json)
    moves_usi = json.loads(moves_json) if moves_json else []
    candidates = json.loads(candidates_json) if candidates_json else []

    worst_sfen = sfen_list[worst_move] if 0 <= worst_move < len(sfen_list) else sfen_list[-1]
    score_before = scores[worst_move - 1] if 1 <= worst_move - 1 < len(scores) else scores[0]
    score_after = scores[worst_move] if 0 <= worst_move < len(scores) else scores[-1]
    score_diff = score_after - score_before

    st.markdown('<div class="game-card">', unsafe_allow_html=True)
    st.subheader(f"{slot}局目")
    st.markdown(f"**対局ID**: `{game_id}`")
    st.markdown(f"**相手**: {opponent}")
    st.markdown(f"**先手**: {black_name} / **後手**: {white_name}")
    st.markdown(f"**最大悪手**: {worst_move} 手目")
    st.markdown(f"**実際に指した手（USI）**: `{actual_move_usi}`")
    st.markdown(f"**実際に指した手（日本語）**: `{actual_move_jp}`")
    st.markdown(f"**評価値変化**: `{score_before} → {score_after}` （差分 `{score_diff}`）")

    left, right = st.columns([1.2, 1.8])

    with left:
        st.markdown("### 悪手局面")
        st.markdown(
            render_board_html(
                worst_sfen,
                black_name=black_name,
                white_name=white_name,
                lastmove_usi=actual_move_usi,
                cell_px=48,
                show_names=True,
                flipped=is_flipped,
            ),
            unsafe_allow_html=True,
        )

    with right:
        st.markdown("### 候補手3つ")
        if not candidates:
            st.info("候補手データがありません。")
        else:
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
                            cell_px=32,
                            show_names=False,
                            flipped=is_flipped,
                        ),
                        unsafe_allow_html=True,
                    )

    st.markdown('</div>', unsafe_allow_html=True)
    st.divider()