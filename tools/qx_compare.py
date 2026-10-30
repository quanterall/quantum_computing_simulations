"""Does our code share code with qx (github.com/richarc/qx)? The method behind CREDITS.md and the site's note.

  1. dependencies: does any qx release depend on Nx, and does our service depend on qx? (mix.exs)
  2. identifiers: our distinctive module, function and field names, searched through all of qx's history
     (every commit, every file under lib/ and examples/)
  3. similarity: each of our Elixir files cut into overlapping 8-token shingles (comments and strings removed,
     identifiers and literals kept as written); the share of our shingles found anywhere in qx's lib/ at each
     release (containment). A positive control runs the same measure between consecutive qx releases, which
     share most of their code — the measure must score those high, or it could not detect shared code.

    ~/quanterall/.venv-tools/bin/python tools/qx_compare.py site/data
"""
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

OUT = Path(sys.argv[1])
OURS = Path(__file__).resolve().parents[1] / "quanterall_quantum" / "lib"
IDENTS = ["QxSim", "QxDemo", "QxStab", "QxProcesses", "visualisation_hint", "counts_before", "syndrome_counts",
          "corrected_counts", "logical_fidelity", "majority_decode", "build_cdf", "bisect", "sample_low",
          "measure_reset", "cmulmod", "frame_batch", "qec_cycle", "run_bitflip_repetition", "tuple_state_vector"]


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True).stdout


def tokens(src):
    src = re.sub(r'"""[\s\S]*?"""', " ", src)                  # heredocs (docs)
    src = re.sub(r'"(?:\\.|[^"\\])*"', ' "" ', src)            # strings
    src = re.sub(r"#[^\n]*", " ", src)                          # comments
    return re.findall(r"[A-Za-z_][A-Za-z0-9_?!]*|\d+(?:\.\d+)?|\|>|->|<-|<<<|>>>|&&&|\|\|\||[^\sA-Za-z0-9_]", src)


def shingles(toks, k=8):
    return {tuple(toks[i:i + k]) for i in range(len(toks) - k + 1)}


def qx_shingles(repo, tag):
    files = [f for f in git(repo, "ls-tree", "-r", "--name-only", tag, "lib").split() if f.endswith((".ex", ".exs"))]
    out = set()
    for f in files:
        out |= shingles(tokens(git(repo, "show", f"{tag}:{f}")))
    return out


with tempfile.TemporaryDirectory() as tmp:
    repo = Path(tmp) / "qx"
    subprocess.run(["git", "clone", "-q", "https://github.com/richarc/qx.git", str(repo)], check=True)
    tags = sorted(git(repo, "tag").split(), key=lambda t: [int(x) for x in t.lstrip("v").split(".")])
    commits = git(repo, "rev-list", "--all").split()

    deps = {t: "{:nx" in git(repo, "show", f"{t}:mix.exs") for t in tags}
    our_mix = (OURS.parent / "mix.exs").read_text()
    hits = {}
    for ident in IDENTS:
        found = subprocess.run(["git", "-C", str(repo), "grep", "-l", "-F", ident, *commits, "--", "lib", "examples"],
                               capture_output=True, text=True).stdout.strip()
        hits[ident] = len(found.splitlines()) if found else 0

    ours = {f.name: shingles(tokens(f.read_text())) for f in sorted(OURS.glob("*.ex"))}
    per_tag = {}
    for t in tags:
        qs = qx_shingles(repo, t)
        per_tag[t] = {name: round(len(s & qs) / max(len(s), 1), 4) for name, s in ours.items()}
    control = {}
    for a, b in zip(tags, tags[1:]):
        sa, sb = qx_shingles(repo, a), qx_shingles(repo, b)
        control[f"{a} in {b}"] = round(len(sa & sb) / max(len(sa), 1), 4)

    result = {
        "qx": "https://github.com/richarc/qx", "tags": tags, "commits_searched": len(commits),
        "qx_depends_on_nx_at_every_tag": all(deps.values()),
        "our_service_depends_on_qx": ("qx" in re.findall(r"\{:(\w+)", our_mix)),
        "identifier_hits_in_qx_history": hits,
        "our_shingles_found_in_qx_max": {name: max(per_tag[t][name] for t in tags) for name in ours},
        "positive_control_consecutive_qx_releases": {"min": min(control.values()), "median": sorted(control.values())[len(control) // 2]},
        "method": "8-token shingles after removing comments, docs and string contents; containment = share of our shingles present in qx's lib/ at a release",
        "per_tag": per_tag,
    }
    json.dump(result, open(OUT / "qx-comparison.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items() if k not in ("per_tag", "tags")}, indent=1))
