// anilist_completed_diff_to_csv.mjs
import { writeFile } from "node:fs/promises";

const ENDPOINT = "https://graphql.anilist.co";
const QUERY = /* GraphQL */ `
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
`;

function csvEscape(val) {
  if (val === null || val === undefined) return "";
  const s = String(val);
  // ダブルクォートと改行を考慮してRFC4180風に
  if (/[",\n]/.test(s)) return `"${s.replace(/"/g, '""')}"`;
  return s;
}

function toCSV(rows, header) {
  const head = header.join(",");
  const body = rows
    .map((r) => header.map((h) => csvEscape(r[h])).join(","))
    .join("\n");
  return `${head}\n${body}\n`;
}

function collectMedia(coll) {
  if (!coll?.lists) return [];
  const medias = [];
  for (const lst of coll.lists) {
    for (const e of lst.entries ?? []) {
      if (e?.media) medias.push(e.media);
    }
  }
  // idでユニーク化
  const uniq = new Map();
  for (const m of medias) uniq.set(m.id, m);
  return [...uniq.values()];
}

async function fetchCompleted(user1, user2) {
  const token = process.env.ANILIST_TOKEN?.trim();
  const headers = {
    "Content-Type": "application/json",
    Accept: "application/json",
  };
  if (token) headers.Authorization = `Bearer ${token}`;

  const res = await fetch(ENDPOINT, {
    method: "POST",
    headers,
    body: JSON.stringify({ query: QUERY, variables: { user1, user2 } }),
    // Node の fetch には timeout option がないため AbortController を使う場合は別途実装
  });

  if (!res.ok) {
    const text = await res.text().catch(() => "");
    throw new Error(`AniList API error: ${res.status} ${res.statusText}\n${text}`);
  }

  const data = await res.json();
  if (data.errors) {
    throw new Error(`AniList GraphQL errors:\n${JSON.stringify(data.errors, null, 2)}`);
  }
  return data.data;
}

async function main() {
  const user1 = process.argv[2] || "tayo6";
  const user2 = process.argv[3] || "namonaki";
  const outPath =
    process.argv[4] || `anilist_diff_${user1}_vs_${user2}.csv`;

  const { u1, u2 } = await fetchCompleted(user1, user2);
  const u1Media = collectMedia(u1);
  const u2Media = collectMedia(u2);

  const set1 = new Set(u1Media.map((m) => m.id));
  const set2 = new Set(u2Media.map((m) => m.id));

  const idx1 = new Map(u1Media.map((m) => [m.id, m]));
  const idx2 = new Map(u2Media.map((m) => [m.id, m]));

  const onlyInUser1 = [...set1].filter((id) => !set2.has(id)).map((id) => idx1.get(id));
  const onlyInUser2 = [...set2].filter((id) => !set1.has(id)).map((id) => idx2.get(id));

  // 出力行作成
  const rows = [];
  for (const m of onlyInUser1.sort((a, b) => a.id - b.id)) {
    rows.push({
      side: `only_in_${user1}`,
      media_id: m.id,
      title_romaji: m.title?.romaji ?? "",
      title_english: m.title?.english ?? "",
      title_native: m.title?.native ?? "",
      site_url: m.siteUrl ?? "",
    });
  }
  for (const m of onlyInUser2.sort((a, b) => a.id - b.id)) {
    rows.push({
      side: `only_in_${user2}`,
      media_id: m.id,
      title_romaji: m.title?.romaji ?? "",
      title_english: m.title?.english ?? "",
      title_native: m.title?.native ?? "",
      site_url: m.siteUrl ?? "",
    });
  }

  const header = ["side", "media_id", "title_romaji", "title_english", "title_native", "site_url"];
  const csv = toCSV(rows, header);
  await writeFile(outPath, csv, "utf8");
  console.log(`Done: ${outPath} (${rows.length} rows)`);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
