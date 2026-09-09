"""Insert the validated time-series application without changing benchmark results."""
import json
from pathlib import Path
import shutil

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    validation=json.loads((HERE/'results/validation.json').read_text())
    assert validation['status']=='numerical_checks_passed'
    assert validation['visual_inspection']=='passed'
    paper=ROOT/'manuscript/barbac_manuscript.md';text=paper.read_text()
    heading='**3.6 Barcode time-series application**'
    end=text.index('**4 Discussion**')
    start=text.find(heading)
    if start==-1:start=end
    paragraph=(HERE/'results/paper_time_series_section.md').read_text().split('\n',1)[1].strip()
    section=heading+'\n\n'+paragraph+'\n\n![](media/time_series_trajectories.png)\n\n'+(
        '**Figure 6. Barcode frequency trajectories in the Chen et al. (2023) hBFA1 YPD assay.** '
        'L01–L03 identify the three barcode pairs with the largest pooled published counts, '
        'shown in both biological assay replicates. Frequencies are conditional on the same '
        '2,314 published pair identifiers for the publication and each barbac configuration. '
        'Markers denote generations 8, 16, 24 and 40; lines connect measurements without fitting '
        'a growth model. Published-pair molecule coverage and unmatched counts are reported '
        'separately in the accompanying time-series analysis. Publication agreement is not '
        'a ground-truth accuracy measurement.\n\n')
    text=text[:start]+section+text[end:]
    reference=('[Chen,V. *et al.* (2023) Evolution of haploid and diploid populations reveals common, '
        'strong, and variable pleiotropic effects in non-home environments. *eLife*, **12**, e92899.]'
        '(https://doi.org/10.7554/eLife.92899)\n\n')
    if 'https://doi.org/10.7554/eLife.92899' not in text:
        text=text.replace('**References**\n\n','**References**\n\n'+reference,1)
    paper.write_text(text)
    for ext in ['png','pdf','svg']:
        shutil.copyfile(HERE/f'results/barcode_trajectories.{ext}',ROOT/f'manuscript/media/time_series_trajectories.{ext}')
    print('Inserted Section 3.6, Figure 6 and the Chen 2023 reference.')


if __name__=='__main__':main()
