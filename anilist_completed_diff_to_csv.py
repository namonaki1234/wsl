# anilist_completed_diff_to_csv.py
import csv
import json
import sys
import requests

ENDPOINT = "https://graphql.anilist.co"
QUERY = """
query CompareCompleted($user1: String!, $user2: String!) {
  u1: MediaListCollection(userName: $user1, type: ANIME, status_in: [COMPLETED]) {
    lists {
      entries {
        media {
          id
          title { romaji english native }
          siteUrl
        }
      }
    }
  }
  u2: MediaListCollection(userName: $user2, type: ANIME, status_in: [COMPLETED]) {
    lists {
      entries {
        media {
          id
          title { romaji english native }
          siteUrl
        }
      }
    }
  }
}
"""

def collect_media(coll):
    """MediaListCollection → [media dict ...]"""
    if not coll or "lists" not in coll or coll["lists"] is None:
        return []
    medias = []
    for lst in coll["lists"]:
        for e in (lst.get("entries") or []):
            m = e.get("media")
            if m:
                medias.append(m)
    # 万一の重複に備えてidでユニーク化
    uniq = {}
    for m in medias:
        uniq[m["id"]] = m
    return list(uniq.values())

def main(user1="tayo6", user2="namonaki", out_path=None):
    if out_path is None:
        out_path = f"anilist_diff_{user1}_vs_{user2}.csv"

    resp = requests.post(
        ENDPOINT,
        json={"query": QUERY, "variables": {"user1": user1, "user2": user2}},
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    if "errors" in data:
        raise RuntimeError(json.dumps(data["errors"], ensure_ascii=False, indent=2))

    u1_media = collect_media(data["data"].get("u1"))
    u2_media = collect_media(data["data"].get("u2"))

    set1 = {m["id"] for m in u1_media}
    set2 = {m["id"] for m in u2_media}

    only_in_user1_ids = set1 - set2
    only_in_user2_ids = set2 - set1

    # id→media の辞書
    idx1 = {m["id"]: m for m in u1_media}
    idx2 = {m["id"]: m for m in u2_media}

    rows = []

    for mid in sorted(only_in_user1_ids):
        m = idx1[mid]
        rows.append({
            "side": "only_in_tayo6",
            "media_id": m["id"],
            "title_romaji": (m["title"]["romaji"] or "").strip() if m.get("title") else "",
            "title_english": (m["title"]["english"] or "").strip() if m.get("title") else "",
            "title_native": (m["title"]["native"] or "").strip() if m.get("title") else "",
            "site_url": m.get("siteUrl", "")
        })

    for mid in sorted(only_in_user2_ids):
        m = idx2[mid]
        rows.append({
            "side": "only_in_namonaki",
            "media_id": m["id"],
            "title_romaji": (m["title"]["romaji"] or "").strip() if m.get("title") else "",
            "title_english": (m["title"]["english"] or "").strip() if m.get("title") else "",
            "title_native": (m["title"]["native"] or "").strip() if m.get("title") else "",
            "site_url": m.get("siteUrl", "")
        })

    # CSVへ出力
    fieldnames = ["side", "media_id", "title_romaji", "title_english", "title_native", "site_url"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Done: {out_path}  ({len(rows)} rows)")

if __name__ == "__main__":
    # 使い方: python script.py [user1] [user2] [out_path]
    u1 = sys.argv[1] if len(sys.argv) > 1 else "tayo6"
    u2 = sys.argv[2] if len(sys.argv) > 2 else "namonaki"
    outp = sys.argv[3] if len(sys.argv) > 3 else None
    main(u1, u2, outp)
