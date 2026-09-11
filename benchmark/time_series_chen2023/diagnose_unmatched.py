"""Describe exact-ID differences without changing any clustering assignments."""
import json
from pathlib import Path

import pandas as pd
from rapidfuzz.distance import Levenshtein

from extract_barcodes import sha256

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
CATEGORIES = ['published_pair', 'unpublished_pair_both_components_known',
              'BC2_absent_BC1_known', 'BC2_known_BC1_absent', 'both_components_absent']


def main():
    provenance = json.loads((OUT / 'provenance.json').read_text())
    rows = []
    for method in provenance['methods']:
        series = pd.read_csv(OUT / f'time_series_{method}.csv.gz')
        published = series.loc[series.is_published_pair, 'Barcode'].unique()
        known2 = {b.split('_')[0] for b in published}
        known1 = {b.split('_')[1] for b in published}
        categories = {}
        for barcode in series.Barcode.unique():
            bc2, bc1 = barcode.split('_')
            categories[barcode] = ('published_pair' if barcode in published else
                'unpublished_pair_both_components_known' if bc2 in known2 and bc1 in known1 else
                'BC2_absent_BC1_known' if bc1 in known1 else
                'BC2_known_BC1_absent' if bc2 in known2 else 'both_components_absent')
        series['category'] = series.Barcode.map(categories)
        for sample, frame in series.groupby('sample'):
            totals = frame.groupby('category').counts.sum().reindex(CATEGORIES, fill_value=0)
            assert totals.sum() == frame.assigned_molecules.iloc[0]
            for category, count in totals.items():
                rows.append(dict(method=method, sample=sample, category=category,
                    molecules=int(count), input_molecules=int(frame.input_molecules.iloc[0]),
                    percent_input=100 * count / frame.input_molecules.iloc[0]))
        if method == 'barbac_lv':
            top = series[~series.is_published_pair].groupby('Barcode').counts.sum().nlargest(20)
            distances = []
            for barcode, count in top.items():
                bc2, bc1 = barcode.split('_')
                distances.append(dict(Barcode=barcode, molecules=int(count), category=categories[barcode],
                    nearest_published_BC2_LV=min(Levenshtein.distance(bc2, b) for b in known2),
                    nearest_published_BC1_LV=min(Levenshtein.distance(bc1, b) for b in known1)))
            pd.DataFrame(distances).to_csv(OUT / 'largest_unmatched_LV_pairs.csv', index=False)
            unmatched = int(series.loc[~series.is_published_pair, 'counts'].sum())
            total = int(series.counts.sum())
            diagnosis = dict(total_LV_molecules=total, unmatched_LV_molecules=unmatched,
                pooled_unmatched_percent=100 * unmatched / total,
                largest_two_unmatched_molecules=int(top.iloc[:2].sum()),
                largest_two_percent_of_unmatched=100 * int(top.iloc[:2].sum()) / unmatched,
                largest_two_nearest_published_BC2_LV=[d['nearest_published_BC2_LV'] for d in distances[:2]],
                interpretation='The largest differences include abundant pairs whose BC2 sequences are outside '
                    'the distance-three radius of every published BC2. This does not establish their biological '
                    'validity or the reason the publication omits them. No diagnostic remapping was applied.')
    pd.DataFrame(rows).to_csv(OUT / 'published_identity_categories.csv', index=False)
    diagnosis['inputs_sha256'] = {f'time_series_{m}.csv.gz': sha256(OUT / f'time_series_{m}.csv.gz')
                                  for m in provenance['methods']}
    diagnosis['script_sha256'] = sha256(Path(__file__))
    (OUT / 'unmatched_diagnosis.json').write_text(json.dumps(diagnosis, indent=2) + '\n')
    print(json.dumps({k: v for k, v in diagnosis.items() if k != 'inputs_sha256'}, indent=2))


if __name__ == '__main__':
    main()
