#!/usr/bin/env python3
"""Render the pipeline figure and barcode-design schematic as SVG, PNG and PDF.

Colours come from the LTC palette ``alger`` bundled in barbac
(``barbac_palette("alger")``): teal, sage, amber, red.

Usage:
    python3 manuscript/figures/render_pipeline_figures.py
    python3 manuscript/figures/render_pipeline_figures.py --design "ATTA-NNNNNN-ATCCC" \
        --design "ATGC-NNN-ATTA-NNN-AAATGC"

Rendering PNG/PDF uses the local Google Chrome in headless mode; SVG files are
written regardless.
"""
import argparse
import html
import pathlib
import subprocess
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

TEAL, SAGE, AMBER, RED = "#1A5B5B", "#ACC8BE", "#F4AB5C", "#D1422F"
INK, BODY, MUTED, RULE = "#14161a", "#2f2c28", "#6f6a63", "#dcd6cb"
FONT = "Arial,Helvetica,sans-serif"
MONO = "Menlo,Consolas,'DejaVu Sans Mono',monospace"
FILE_COLOURS = {"fastq": TEAL, "bam": AMBER, "csv": RED, "html": SAGE, "png": SAGE, "pdf": AMBER}


def esc(s):
    return html.escape(s, quote=True)


# ----------------------------------------------------------------- pipeline
def file_icon(x, y, ext, scale=1.0):
    c = FILE_COLOURS[ext]
    w, h = 30 * scale, 37 * scale
    fold = 9 * scale
    return (
        f'<path d="M{x} {y}h{w - fold}l{fold} {fold}v{h - fold}H{x}z" fill="{c}"/>'
        f'<path d="M{x + w - fold} {y}v{fold}h{fold}z" fill="#fff" fill-opacity="0.6"/>'
        f'<text x="{x + w / 2}" y="{y + h + 15 * scale}" text-anchor="middle" font-family="{MONO}" '
        f'font-size="{12 * scale}" fill="{BODY}">.{ext}</text>'
    )


def stage(x, y, w, h, title, fn=None, tool=None, lines=(), files=(), files_y=None):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" fill="#fff" stroke="{TEAL}" stroke-width="2.5"/>']
    out.append(f'<path d="M{x + 14} {y}H{x + w - 14}Q{x + w} {y} {x + w} {y + 14}V{y + 34}H{x}V{y + 14}Q{x} {y} {x + 14} {y}z" fill="{TEAL}"/>')
    out.append(f'<text x="{x + 18}" y="{y + 23}" font-family="{FONT}" font-size="17" font-weight="bold" fill="#fff">{esc(title)}</text>')
    cy = y + 64
    if fn:
        fw = 10.5 * len(fn) + 18
        out.append(f'<rect x="{x + 16}" y="{cy - 17}" width="{fw}" height="26" rx="6" fill="#FBE7E3"/>')
        out.append(f'<text x="{x + 25}" y="{cy + 1}" font-family="{MONO}" font-size="16" font-weight="bold" fill="{RED}">{esc(fn)}()</text>')
        cy += 34
    if tool:
        out.append(f'<text x="{x + 18}" y="{cy}" font-family="{FONT}" font-size="16" font-style="italic" font-weight="bold" fill="{INK}">{esc(tool)}</text>')
        cy += 28
    for line in lines:
        if isinstance(line, tuple):  # (emph, rest)
            emph, rest = line
            out.append(f'<circle cx="{x + 24}" cy="{cy - 5}" r="2.6" fill="{TEAL}"/>')
            out.append(f'<text x="{x + 34}" y="{cy}" font-family="{FONT}" font-size="15" fill="{BODY}">'
                       f'<tspan font-style="italic" font-weight="bold">{esc(emph)}</tspan> {esc(rest)}</text>')
        else:
            out.append(f'<text x="{x + 18}" y="{cy}" font-family="{FONT}" font-size="15" fill="{MUTED}">{esc(line)}</text>')
        cy += 23
    if files:
        fy = files_y if files_y is not None else y + h - 64
        fx = x + 20
        for ext in files:
            out.append(file_icon(fx, fy, ext))
            fx += 54
    return "\n".join(out)


def arrow_right(x, y):
    return f'<path d="M{x} {y - 16}L{x + 22} {y}L{x} {y + 16}z" fill="{TEAL}"/>'


def arrow_down(x, y):
    return f'<path d="M{x - 16} {y}L{x + 16} {y}L{x} {y + 22}z" fill="{TEAL}"/>'


def pipeline_svg():
    W, H = 1300, 860
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="{178 * H / W:.1f}mm" viewBox="0 0 {W} {H}" role="img" aria-labelledby="t d">',
         '<title id="t">barbac pipeline: Illumina reads or barcode counts to lineage trajectories</title>',
         '<desc id="d">Reads pass quality control, optional pair merging, mapping and barcode extraction. '
         'Extracted counts or an independent count table enter native abundance-aware clustering. '
         'Clusters feed statistics and lineage-frequency plots over time.</desc>',
         f'<rect width="{W}" height="{H}" fill="#fff"/>']

    # Inputs
    s.append(f'<rect x="40" y="40" width="215" height="92" rx="14" fill="#fff" stroke="{TEAL}" stroke-width="2.5"/>')
    s.append(file_icon(58, 62, "fastq"))
    s.append(f'<text x="104" y="76" font-family="{FONT}" font-size="16" font-weight="bold" fill="{INK}">Illumina reads</text>')
    s.append(f'<text x="104" y="98" font-family="{FONT}" font-size="14" fill="{MUTED}">R1-only or overlapping</text>')
    s.append(f'<text x="104" y="116" font-family="{FONT}" font-size="14" fill="{MUTED}">R1 + R2</text>')

    s.append(f'<rect x="805" y="40" width="200" height="92" rx="14" fill="#fff" stroke="{TEAL}" stroke-width="2.5"/>')
    s.append(file_icon(823, 62, "csv"))
    s.append(f'<text x="869" y="76" font-family="{FONT}" font-size="16" font-weight="bold" fill="{INK}">Barcode counts</text>')
    s.append(f'<text x="869" y="98" font-family="{FONT}" font-size="14" fill="{MUTED}">enter clustering</text>')
    s.append(f'<text x="869" y="116" font-family="{FONT}" font-size="14" fill="{MUTED}">without FASTQ</text>')

    # Stages
    s.append(stage(40, 180, 215, 300, "Quality control", fn="barbac_map", tool="FastQC",
                   lines=["Merge overlapping pairs", "PEAR; R1-only skips this"], files=("html", "fastq")))
    s.append(stage(300, 240, 190, 200, "Map to reference", tool="minimap2 · samtools",
                   lines=["primary alignments", "sorted, indexed"], files=("bam",)))
    s.append(stage(300, 505, 200, 215, "Mapping stats", tool="samtools · MultiQC",
                   lines=["mapped, unmapped reads", "command timings"], files=("csv",)))
    s.append(stage(545, 180, 205, 300, "Extract barcodes", fn="barbac_extract",
                   lines=[("Extract", "by coordinates"), "or flank pattern",
                          ("Keep", "observed indels"),
                          ("Count", "distinct barcodes"), "and record their length"], files=("csv",)))
    s.append(stage(545, 515, 205, 285, "Barcode stats", fn="barbac_extract_stats",
                   lines=[("Length", "distribution"), ("Count", "distribution"), ("Entropy", "distribution"),
                          "optional length filter"], files=("csv", "png")))
    s.append(stage(805, 180, 200, 365, "Cluster barcodes", fn="barbac_cluster",
                   lines=[("Pool", "timepoints per"), "population (optional)",
                          ("Rank", "by abundance;"), "assign to best centroid", "within distance d",
                          ("Guard", "with count ratio"), "and error rate",
                          ("Report", "centroids and"), "every member"], files=("csv",)))
    s.append(stage(805, 582, 205, 210, "Clustering stats", fn="cluster_stats",
                   lines=["clusters, total reads", "singletons, largest cluster", "top-K abundance share"], files=("csv",)))
    s.append(stage(1060, 205, 200, 270, "Time series", fn="barbac_ts_area",
                   lines=[("Plot", "lineage"), "frequencies over time",
                          ("Choose", "one of 32 LTC"), "palettes"], files=("png", "pdf")))

    # Arrows
    for x in (262, 503, 763, 1018):
        s.append(arrow_right(x, 330))
    for x, y in ((147, 136), (905, 136), (395, 444), (647, 484), (905, 549)):
        s.append(arrow_down(x, y))

    s.append(f'<text x="40" y="{H - 22}" font-family="{FONT}" font-size="14" fill="{MUTED}">'
             'Distance d uses Hamming or Levenshtein. Colours: LTC palette alger.</text>')
    s.append("</svg>")
    return "\n".join(s)


# ------------------------------------------------------------------ designs
def parse_design(design):
    """Split 'ATTA-NNNNNN-ATCCC' into [(seq, is_variable), ...]."""
    blocks = []
    for part in design.strip().upper().split("-"):
        if not part:
            continue
        blocks.append((part, set(part) == {"N"}))
    return blocks


def design_svg(designs, tile=36, gap=4, block_gap=16):
    rows = [parse_design(d) for d in designs]
    longest = max(sum(len(b) for b, _ in r) * (tile + gap) + (len(r) - 1) * block_gap for r in rows)
    label_w = 60
    W = label_w + longest + 60
    row_h = 118
    H = 30 + row_h * len(rows) + 70
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="178mm" height="{178 * H / W:.1f}mm" viewBox="0 0 {W} {H}" role="img" aria-labelledby="t d">',
         '<title id="t">Designed barcodes with fixed anchors and variable positions</title>',
         '<desc id="d">' + esc("; ".join(designs)) + '. Fixed anchor bases are shown in teal, variable N positions in amber.</desc>',
         f'<rect width="{W}" height="{H}" fill="#fff"/>']
    for ri, (design, blocks) in enumerate(zip(designs, rows)):
        y = 40 + ri * row_h
        letter = chr(ord("A") + ri)
        s.append(f'<text x="20" y="{y + tile / 2 + 7}" font-family="{FONT}" font-size="22" font-weight="bold" fill="{INK}">{letter}</text>')
        x = label_w
        n_var = sum(len(b) for b, v in blocks if v)
        n_fix = sum(len(b) for b, v in blocks if not v)
        for seq, variable in blocks:
            bx = x
            fill = AMBER if variable else TEAL
            txt = INK if variable else "#fff"
            for base in seq:
                s.append(f'<rect x="{x}" y="{y}" width="{tile}" height="{tile}" rx="6" fill="{fill}"/>')
                s.append(f'<text x="{x + tile / 2}" y="{y + tile / 2 + 7}" text-anchor="middle" font-family="{MONO}" '
                         f'font-size="20" font-weight="bold" fill="{txt}">{base}</text>')
                x += tile + gap
            bw = x - gap - bx
            # bracket + label under the block
            by = y + tile + 10
            s.append(f'<path d="M{bx} {by}v6h{bw}v-6" fill="none" stroke="{MUTED}" stroke-width="1.5"/>')
            lab = f"{len(seq)} variable" if variable else f"anchor {seq}"
            s.append(f'<text x="{bx + bw / 2}" y="{by + 24}" text-anchor="middle" font-family="{FONT}" font-size="13" fill="{MUTED}">{esc(lab)}</text>')
            x += block_gap
        total = n_var + n_fix
        s.append(f'<text x="{label_w}" y="{y + tile + 58}" font-family="{FONT}" font-size="14" fill="{BODY}">'
                 f'<tspan font-family="{MONO}" font-weight="bold">{esc(design.upper())}</tspan>'
                 f'   {total} nt: {n_fix} fixed, {n_var} variable</text>')
    # legend
    ly = H - 34
    s.append(f'<rect x="{label_w}" y="{ly - 14}" width="18" height="18" rx="4" fill="{TEAL}"/>')
    s.append(f'<text x="{label_w + 26}" y="{ly}" font-family="{FONT}" font-size="14" fill="{BODY}">fixed anchor base (same in every barcode)</text>')
    s.append(f'<rect x="{label_w + 330}" y="{ly - 14}" width="18" height="18" rx="4" fill="{AMBER}"/>')
    s.append(f'<text x="{label_w + 356}" y="{ly}" font-family="{FONT}" font-size="14" fill="{BODY}">N: variable position (A, C, G or T)</text>')
    s.append("</svg>")
    return "\n".join(s)


# ---------------------------------------------------------------- rendering
def svg_size_mm(svg_text):
    import re
    w = float(re.search(r'width="([\d.]+)mm"', svg_text).group(1))
    h = float(re.search(r'height="([\d.]+)mm"', svg_text).group(1))
    return w, h


def render(svg_path, dpi=400):
    svg_text = svg_path.read_text()
    w_mm, h_mm = svg_size_mm(svg_text)
    with tempfile.TemporaryDirectory(prefix="barbac-fig-") as tmp:
        tmp = pathlib.Path(tmp)
        page = tmp / "page.html"
        page.write_text(
            "<!DOCTYPE html><html><head><meta charset='utf-8'><style>"
            f"@page{{size:{w_mm}mm {h_mm}mm;margin:0}}"
            "html,body{margin:0;padding:0;background:#fff}"
            f"svg{{display:block;width:{w_mm}mm;height:{h_mm}mm}}"
            "</style></head><body>" + svg_text + "</body></html>")
        base = [CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars", "--no-pdf-header-footer"]
        subprocess.run(base + [f"--print-to-pdf={svg_path.with_suffix('.pdf')}", f"file://{page}"],
                       check=True, capture_output=True)
        px_w = round(w_mm / 25.4 * dpi)
        px_h = round(h_mm / 25.4 * dpi)
        css_w = round(w_mm / 25.4 * 96)
        css_h = round(h_mm / 25.4 * 96)
        scale = px_w / css_w
        subprocess.run(base + [f"--force-device-scale-factor={scale:.4f}", f"--window-size={css_w},{css_h}",
                               f"--screenshot={svg_path.with_suffix('.png')}", f"file://{page}"],
                       check=True, capture_output=True)
        return px_w, px_h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", action="append",
                    help="Barcode design such as ATTA-NNNNNN-ATCCC; repeat for several rows")
    ap.add_argument("--outdir", default=str(HERE))
    ap.add_argument("--no-render", action="store_true", help="Write SVG only")
    args = ap.parse_args()
    designs = args.design or ["ATTA-NNNNNN-ATCCC", "ATGC-NNN-ATTA-NNN-AAATGC"]
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    targets = {outdir / "pipeline.svg": pipeline_svg(), outdir / "barcode_designs.svg": design_svg(designs)}
    for path, text in targets.items():
        path.write_text(text)
        print("wrote", path)
        if not args.no_render:
            px = render(path)
            print("rendered", path.with_suffix(".pdf"), path.with_suffix(".png"), f"{px[0]}x{px[1]} px")


if __name__ == "__main__":
    main()
