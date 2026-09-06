"""Build the shared labelling page for the closing-accuracy test: one HTML file with every chip embedded as a JPEG
data URI and a small app that records 1 / 2 / 3 per case per labeller in the artifact's database.
Reads out/11_chips/manifest.json + out/11_chips/jpg/*.jpg; writes out/11_label_page.html.
Run:  $(cat .python_env) scripts/11d_label_page.py"""
import os, json, base64, hashlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.join(ROOT, "out"); CH = os.path.join(OUT, "11_chips")
TILE = os.environ.get("TILE", "19_39"); man = {k: v for k, v in json.load(open(os.path.join(CH, "manifest.json"))).items() if v.get("tile") == TILE}
cases = []
for cid, m in man.items():
    jp = os.path.join(CH, "jpg", f"{cid}.jpg")
    if not os.path.exists(jp): continue
    cases.append(dict(id=cid, tile=m["tile"], pool=m["pool"], pieces=m["n_pieces"], gap=m["gap_m"], sites=m["site_ids"], dates=[m["dates"][str(y)] for y in m["seasons"]],
                      side=int(m["side_m"]), img="data:image/jpeg;base64," + base64.b64encode(open(jp, "rb").read()).decode()))
cases.sort(key=lambda c: hashlib.md5(c["id"].encode()).hexdigest())  # one fixed shuffled order for every labeller
html = r"""<title>Lake Closing Test __TILE__</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{--ground:#f1f5f9;--surface:#ffffff;--line:#d5dee7;--ink:#16232e;--muted:#5b6b79;--accent:#1e6fa8;--accent-ink:#ffffff;--one:#1f7a5a;--two:#b8641a;--unsure:#6b7280;--btn:#f7f9fb;--focus:#9cc4e4;}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--ground:#0f171e;--surface:#17222b;--line:#2b3946;--ink:#e6edf3;--muted:#93a3b1;--accent:#4a9fdc;--accent-ink:#0b1218;--one:#4fbf93;--two:#e39a4c;--unsure:#9aa3ad;--btn:#1d2a34;--focus:#2f6a94;}}
:root[data-theme="dark"]{--ground:#0f171e;--surface:#17222b;--line:#2b3946;--ink:#e6edf3;--muted:#93a3b1;--accent:#4a9fdc;--accent-ink:#0b1218;--one:#4fbf93;--two:#e39a4c;--unsure:#9aa3ad;--btn:#1d2a34;--focus:#2f6a94;}
body{background:var(--ground);color:var(--ink);font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:15px;line-height:1.45;margin:0}
.mono{font-family:"IBM Plex Mono",ui-monospace,Menlo,monospace;font-variant-numeric:tabular-nums}
header{display:flex;align-items:center;gap:20px;padding:12px 22px;border-bottom:1px solid var(--line);background:var(--surface)}
header h1{font-size:17px;font-weight:600;margin:0;letter-spacing:.01em}
.who{margin-left:auto;display:flex;align-items:center;gap:10px;color:var(--muted);font-size:13px}
.who b{color:var(--ink);font-weight:500}
.prog{display:flex;align-items:center;gap:10px;font-size:13px;color:var(--muted)}
.bar{width:180px;height:6px;background:var(--line);border-radius:3px;overflow:hidden}.bar i{display:block;height:100%;background:var(--accent);width:0}
main{max-width:1560px;margin:0 auto;padding:16px 22px 40px;display:flex;flex-direction:column;gap:14px}
.chip{background:var(--surface);border:1px solid var(--line);border-radius:6px;padding:8px;overflow-x:auto}
.chip img{display:block;width:100%;min-width:900px;height:auto}
.ask{display:grid;grid-template-columns:1fr auto;gap:18px;align-items:center}
.q{font-size:18px;font-weight:500;text-wrap:balance}
.q small{display:block;font-size:13px;font-weight:400;color:var(--muted);margin-top:4px}
.btns{display:flex;gap:10px}
button{font:inherit;cursor:pointer;border:1px solid var(--line);background:var(--btn);color:var(--ink);border-radius:6px;padding:12px 16px;display:flex;align-items:center;gap:10px;min-width:170px}
button:focus-visible{outline:3px solid var(--focus);outline-offset:1px}
button kbd{font-family:"IBM Plex Mono",monospace;font-size:13px;border:1px solid var(--line);border-radius:4px;padding:1px 7px;background:var(--surface);color:var(--muted)}
button.on{border-width:2px}
.b1.on{border-color:var(--one)}.b2.on{border-color:var(--two)}.b3.on{border-color:var(--unsure)}
.b1 kbd{color:var(--one)}.b2 kbd{color:var(--two)}.b3 kbd{color:var(--unsure)}
.meta{display:flex;flex-wrap:wrap;gap:6px 22px;font-size:13px;color:var(--muted)}
.meta span b{color:var(--ink);font-weight:500}
nav{display:flex;align-items:center;gap:10px;font-size:13px;color:var(--muted)}
nav button{min-width:0;padding:6px 12px}
.dots{display:flex;flex-wrap:wrap;gap:3px;margin-left:auto}
.dots i{width:9px;height:9px;border-radius:2px;background:var(--line);display:block;cursor:pointer}
.dots i.l1{background:var(--one)}.dots i.l2{background:var(--two)}.dots i.l3{background:var(--unsure)}.dots i.cur{outline:2px solid var(--accent);outline-offset:1px}
.gate{max-width:560px;margin:60px auto;background:var(--surface);border:1px solid var(--line);border-radius:8px;padding:26px}
.gate h2{margin:0 0 6px;font-size:18px}.gate p{margin:0 0 16px;color:var(--muted)}
.gate .btns button{min-width:120px}
.note{font-size:13px;color:var(--muted)}
.notice{font-size:13px;color:var(--two)}
textarea{font:inherit;font-size:13px;width:100%;box-sizing:border-box;border:1px solid var(--line);border-radius:6px;background:var(--surface);color:var(--ink);padding:8px;resize:vertical;min-height:38px}
.help{background:var(--surface);border:1px solid var(--line);border-radius:6px;padding:12px 16px;font-size:14px;color:var(--muted);max-width:70ch}
.help b{color:var(--ink);font-weight:500}
</style>
<header><h1>Lake closing test, tile __TILE__</h1><div class="prog"><span id="pcount" class="mono">0 / 0</span><div class="bar"><i id="pbar"></i></div></div><div class="who" id="who"></div></header>
<div id="gate" class="gate"><h2>Who is labelling?</h2><p>Each person's answers are kept separately, and nobody sees the other's until the end.</p>
<div class="btns"><button data-n="josh">Josh</button><button data-n="claude">Claude</button></div></div>
<main id="app" hidden>
<div class="help">Every case is a place where the rule had to decide. <b>Red</b> is the site polygon (or two polygons, for a pair). The left panel is the terrain with each season's outline in colour; the four photos are the wettest clear day of the four wettest seasons. Judge the <b>basin</b>, not one day's water: an ice lid splitting a lake is still one lake. Keys <b>1</b>, <b>2</b>, <b>3</b> answer and advance; <b>←</b> <b>→</b> move; <b>n</b> jumps to the next unanswered.</div>
<div class="chip"><img id="img" alt="case chip"></div>
<div class="ask"><div class="q" id="q"></div>
<div class="btns"><button class="b1" data-l="1"><kbd>1</kbd>One lake</button><button class="b2" data-l="2"><kbd>2</kbd>Two or more lakes</button><button class="b3" data-l="3"><kbd>3</kbd>Can't tell</button></div></div>
<textarea id="note" placeholder="Optional note (saved with the answer)"></textarea>
<div class="meta" id="meta"></div>
<nav><button id="prev">← Prev</button><button id="next">Next →</button><button id="skip">Next unanswered (n)</button><span id="status" class="note"></span><div class="dots" id="dots"></div></nav>
</main>
<script>
const CASES = __CASES__;
let i = 0, who = null, db = null, labels = {}, notes = {};
const $ = s => document.querySelector(s);
const key = c => who + "__" + c.id;
function render(){
  const c = CASES[i]; $("#img").src = c.img;
  const isPair = c.pool === "apart";
  $("#q").innerHTML = (isPair ? "These two sites are <b>" + c.gap + " m</b> apart. Are they one lake basin, or separate lakes?" : "The rule made this <b>one site</b>. Is all of this water one lake basin, or more than one lake?") +
    "<small>" + {joined: "The 150 m closing glued " + c.pieces + " separate pieces of water into this site.", lid: "In at least one season this site held " + c.pieces + " separate water bodies.", apart: "The closing left these apart; they are the nearest neighbours in the tile."}[c.pool] + "</small>";
  $("#meta").innerHTML = "<span>case <b class=mono>" + c.id + "</b></span><span>window <b>" + (c.side/1000).toFixed(1) + " km</b></span><span>photos <b>" + c.dates.join(", ") + "</b></span><span>" + (i+1) + " of " + CASES.length + "</span>";
  document.querySelectorAll(".btns button[data-l]").forEach(b => b.classList.toggle("on", labels[key(c)] === +b.dataset.l));
  $("#note").value = notes[key(c)] || "";
  const done = CASES.filter(x => labels[key(x)]).length; $("#pcount").textContent = done + " / " + CASES.length; $("#pbar").style.width = (100*done/CASES.length) + "%";
  $("#dots").innerHTML = CASES.map((x, k) => "<i class='" + (labels[key(x)] ? "l" + labels[key(x)] : "") + (k === i ? " cur" : "") + "' data-k='" + k + "' title='" + x.id + "'></i>").join("");
  document.querySelectorAll("#dots i").forEach(d => d.onclick = () => { i = +d.dataset.k; render(); });
  window.scrollTo(0, 0);
}
async function save(l){
  const c = CASES[i]; labels[key(c)] = l; const note = $("#note").value.trim(); notes[key(c)] = note;
  const doc = {labeller: who, case_id: c.id, tile: c.tile, pool: c.pool, label: l, note, ts: new Date().toISOString()};
  try { localStorage.setItem("lbl_" + key(c), JSON.stringify(doc)); } catch(e){}
  if (db) { try { await db.doc("labels/" + key(c)).set(doc); $("#status").textContent = "saved"; } catch(e){ $("#status").textContent = "not saved to the shared store (" + (e.code || e.message) + "); kept in this browser"; $("#status").className = "notice"; } }
  const nxt = CASES.findIndex((x, k) => k > i && !labels[key(x)]); i = nxt >= 0 ? nxt : CASES.findIndex(x => !labels[key(x)]); if (i < 0) i = CASES.length - 1; render();
}
async function start(name){
  who = name; try { localStorage.setItem("labeller", who); } catch(e){}
  $("#who").innerHTML = "labelling as <b>" + who + "</b> <button id='sw' style='min-width:0;padding:3px 9px;font-size:12px'>switch</button>"; $("#sw").onclick = () => { try { localStorage.removeItem("labeller"); } catch(e){} location.reload(); };
  for (const c of CASES) { try { const s = localStorage.getItem("lbl_" + key(c)); if (s) { const d = JSON.parse(s); labels[key(c)] = d.label; notes[key(c)] = d.note || ""; } } catch(e){} }
  $("#gate").hidden = true; $("#app").hidden = false; i = Math.max(0, CASES.findIndex(c => !labels[key(c)])); render();
  db = await claude.use("db");
  if (!db) { $("#status").textContent = "shared store unavailable here; answers stay in this browser only"; $("#status").className = "notice"; return; }
  try { const snap = await db.collection("labels").where("labeller", "==", who).get();
    snap.docs.forEach(d => { const x = d.data(); if (x && x.case_id) { labels[who + "__" + x.case_id] = x.label; notes[who + "__" + x.case_id] = x.note || ""; } });
    $("#status").textContent = "shared store connected"; i = Math.max(0, CASES.findIndex(c => !labels[key(c)])); render();
  } catch(e) { $("#status").textContent = "could not read earlier answers (" + (e.code || e.message) + ")"; $("#status").className = "notice"; }
}
document.querySelectorAll("#gate button").forEach(b => b.onclick = () => start(b.dataset.n));
document.querySelectorAll(".btns button[data-l]").forEach(b => b.onclick = () => save(+b.dataset.l));
$("#prev").onclick = () => { i = (i + CASES.length - 1) % CASES.length; render(); }; $("#next").onclick = () => { i = (i + 1) % CASES.length; render(); };
$("#skip").onclick = () => { const n = CASES.findIndex((x, k) => k > i && !labels[key(x)]); if (n >= 0) { i = n; render(); } };
document.addEventListener("keydown", e => { if (!who || e.target.tagName === "TEXTAREA") return; if (e.key === "1" || e.key === "2" || e.key === "3") save(+e.key); else if (e.key === "ArrowLeft") $("#prev").click(); else if (e.key === "ArrowRight") $("#next").click(); else if (e.key === "n") $("#skip").click(); });
try { const s = localStorage.getItem("labeller"); if (s) start(s); } catch(e){}
</script>
"""
out = os.path.join(OUT, f"11_label_page_{TILE}.html"); open(out, "w").write(html.replace("__CASES__", json.dumps(cases)).replace("__TILE__", TILE))
print(f"{len(cases)} cases, {os.path.getsize(out)/1e6:.1f} MB -> {out}")
