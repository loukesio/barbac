"""Export the benchmark designs as a publication schematic (not observed reads)."""
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/barbac-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
DEST = HERE.parents[1] / 'manuscript/media/benchmark_datasets'


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'svg.fonttype': 'none', 'svg.hashsalt': 'barbac-datasets', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(12, 7.2))
    fig.subplots_adjust(left=.025, right=.985, top=.97, bottom=.03)
    ax.set(xlim=(0, 120), ylim=(0, 100)); ax.axis('off')
    ink, blue, fixed = '#242424', '#d9e7f1', '#e5e5e5'
    ax.text(1, 96, 'Five benchmark datasets', fontsize=18, weight='bold', color=ink)
    ax.text(1, 91, 'Designed barcode architecture, error regimes and input scale', color=ink)
    ax.text(1, 83, 'Dataset', weight='bold')
    ax.text(25, 83, 'Barcode design (N = variable base)', weight='bold')
    ax.text(76, 83, 'Errors and scale', weight='bold')
    designs = [
        ('R-S', 'Random / substitutions', [('N'*20, True)],
         'Substitutions: 0.5% per base\n10,000 truth barcodes; 1 million reads / seed'),
        ('R-I', 'Random / mixed errors', [('N'*20, True)],
         'Substitutions + insertions + deletions\n0.5% each; same scale as R-S'),
        ('A-S', 'Anchored / substitutions', [('N'*8, True), ('ATGC', False), ('N'*8, True), ('ATCGTTAA', False)],
         'Substitutions: 0.5% per base\n10,000 truth barcodes; 1 million reads / seed'),
        ('A-I', 'Anchored / mixed errors', [('N'*8, True), ('ATGC', False), ('N'*8, True), ('ATCGTTAA', False)],
         'Substitutions + insertions + deletions\n0.5% each; same scale as A-S'),
        ('Milos', 'Johnson reference simulation', [('N'*20, True)],
         'Substitutions + sparse length-changing errors\n100,000 truth barcodes; 24,996,128 reads'),
    ]
    for i, (code, label, blocks, notes) in enumerate(designs):
        y = 74 - i*13
        ax.axhline(y+5, xmin=.01, xmax=.99, color='#dddddd', lw=.7)
        ax.text(1, y, code, weight='bold', fontsize=12)
        ax.text(1, y-3.7, label, fontsize=8.5)
        x = 25
        for seq, variable in blocks:
            width = len(seq)*1.65
            ax.add_patch(Rectangle((x, y-1), width, 5, facecolor=blue if variable else fixed,
                                   edgecolor=ink, linewidth=.65))
            ax.text(x+width/2, y+1.5, seq, ha='center', va='center', family='monospace', fontsize=8.5)
            x += width
        length = sum(len(seq) for seq, _ in blocks)
        variable = sum(len(seq) for seq, v in blocks if v)
        ax.text(25, y-4, f'{length} bases; {variable} variable positions', fontsize=9)
        ax.text(76, y-.1, notes, va='center', fontsize=9, linespacing=1.5)
    ax.text(1, 11, 'R-S, R-I, A-S, A-I: independent seeds 42, 43, 44; lognormal abundance (sigma = 1.5).', fontsize=9)
    ax.text(1, 7, 'Milos: one fixed input; 1,544,850 observed unique sequences. 1,001 reads differ in length from their parent.', fontsize=9)
    ax.text(1, 3, 'Schematics show the truth designs. Indels change observed lengths; matching length does not exclude an indel history.', fontsize=9)
    for ext in ('png', 'svg', 'pdf'):
        metadata = {'Date': None} if ext == 'svg' else (
            {'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
        fig.savefig(DEST.with_suffix('.'+ext), dpi=240, facecolor='white', metadata=metadata)
    svg = DEST.with_suffix('.svg')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__ == '__main__':
    main()
