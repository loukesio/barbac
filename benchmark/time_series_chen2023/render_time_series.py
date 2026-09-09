"""Render scientific figures and an evidence-backed time-series results report."""
import json
import os
from pathlib import Path

HERE=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'generated/time_series/cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT=HERE/'results'
LABELS={'barbac_lv':'barbac LV + Poisson','barbac_hamming':'barbac Hamming',
        'shepherd':'Shepherd (e=0.005)','starcode_sphere':'Starcode sphere',
        'starcode_mp':'Starcode message passing','bartender':'Bartender'}
BLUE='#255D8C';GOLD='#B07D22';GRAY='#626262'


def save(fig,name):
    for ext in ['png','pdf','svg']:
        path=OUT/f'{name}.{ext}'
        fig.savefig(path,dpi=200,facecolor='white',bbox_inches='tight')
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def style(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',color='#e6e6e6',linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    provenance=json.loads((OUT/'provenance.json').read_text())
    assert provenance['status']=='complete'
    data=pd.read_csv(OUT/'selected_trajectories.csv')
    summary=pd.read_csv(OUT/'method_summary.csv').set_index('method')
    extraction=pd.read_csv(OUT/'extraction_summary.csv')
    diagnosis=json.loads((OUT/'unmatched_diagnosis.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,
        'axes.titlesize':10,'text.color':'#242424','axes.labelcolor':'#242424',
        'xtick.color':'#424242','ytick.color':'#424242','svg.hashsalt':'chen2023-time-series'})
    top=['L01','L02','L03'];fig,axes=plt.subplots(3,2,figsize=(10.4,8),sharex=True)
    for i,lineage in enumerate(top):
        for j,replicate in enumerate(['R1','R2']):
            ax=axes[i,j];style(ax)
            sub=data[data.display_id.eq(lineage)&data.replicate.eq(replicate)]
            pub=sub[sub.method.eq('barbac_lv')].sort_values('generation')
            ax.plot(pub.generation,100*pub.published_frequency,color=GRAY,linestyle='--',
                    marker='s',markerfacecolor='white',linewidth=1.4,markersize=4,label='Publication')
            for method,color,marker,linestyle in [('barbac_lv',BLUE,'o','-'),('barbac_hamming',GOLD,'^',':')]:
                x=sub[sub.method.eq(method)].sort_values('generation')
                assert len(x)==4
                ax.plot(x.generation,100*x.frequency_published_set,color=color,linestyle=linestyle,
                        marker=marker,markersize=4,linewidth=1.4,label=LABELS[method])
            ax.set_title(f'{lineage} · assay replicate {replicate}',loc='left')
            ax.set_ylim(bottom=0);ax.set_xticks([8,16,24,40])
            if j==0:ax.set_ylabel('Frequency (%)')
            if i==2:ax.set_xlabel('Generation')
        # Both replicates of a lineage use the same vertical scale.
        maximum=max(axes[i,j].get_ylim()[1] for j in (0,1))
        for j in (0,1):axes[i,j].set_ylim(0,maximum)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper left',bbox_to_anchor=(.085,.927),ncol=3,frameon=False)
    fig.suptitle('Barcode frequencies across four sampled generations',x=.085,y=.995,ha='left',fontsize=15)
    fig.text(.085,.95,'Chen et al. 2023 · hBFA1 / YPD · normalized within the same 2,314 published barcode pairs',fontsize=9)
    fig.text(.085,.015,'Three pairs with the highest pooled published counts. Lines connect observations; no growth model is fitted.',fontsize=8,color=GRAY)
    fig.subplots_adjust(top=.855,bottom=.09,left=.085,right=.98,hspace=.37,wspace=.23)
    save(fig,'barcode_trajectories')

    lv=pd.read_csv(OUT/'time_series_barbac_lv.csv.gz')
    abundance=lv.groupby('Barcode').counts.sum().sort_values(ascending=False)
    largest=abundance.head(4).index.tolist()
    known=dict(zip(data.Barcode,data.display_id));lookup=[]
    for i,bc in enumerate(largest):
        lookup.append(dict(Barcode=bc,display_id=known.get(bc,f'U{i+1:02}'),selection='Top four pooled LV pairs'))
    for id in top:
        bc=data.loc[data.display_id.eq(id),'Barcode'].iloc[0]
        if bc not in largest:lookup.append(dict(Barcode=bc,display_id=id,selection='Top three pooled published pairs'))
    pd.DataFrame(lookup).to_csv(OUT/'figure_lineages.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(10.4,4.6),sharey=True)
    colors=[BLUE,GOLD,'#7B8245','#BD668A','#d4d4d4']
    ids={x['Barcode']:x['display_id'] for x in lookup}
    for ax,replicate in zip(axes,['R1','R2']):
        style(ax);sub=lv[lv.replicate.eq(replicate)]
        pivot=sub.pivot(index='Barcode',columns='generation',values='frequency_assigned').fillna(0)
        values=pivot.reindex(largest,fill_value=0).to_numpy()*100
        rest=100-values.sum(axis=0);assert (rest>=-1e-8).all()
        bottom=np.zeros(4)
        for i,heights in enumerate([*values,rest]):
            label=ids[largest[i]] if i<4 else 'Other barcode pairs'
            ax.bar(range(4),heights,bottom=bottom,color=colors[i],edgecolor='white',linewidth=.65,width=.7,label=label)
            bottom+=heights
        assert np.allclose(bottom,100)
        ax.set_xticks(range(4),[8,16,24,40]);ax.set_ylim(0,100)
        ax.set_title(f'Assay replicate {replicate}',loc='left');ax.set_xlabel('Generation (sampled)')
    axes[0].set_ylabel('Assigned molecules (%)')
    fig.suptitle('Composition of the LV barbac time series',x=.085,y=.99,ha='left',fontsize=15)
    fig.text(.085,.91,'All LV-assigned molecules · four largest pooled barcode pairs plus Other · distance 3',fontsize=9)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.53,.01),ncol=5,frameon=False)
    fig.subplots_adjust(top=.8,bottom=.23,left=.085,right=.98,wspace=.2)
    save(fig,'barcode_composition')

    times=summary.sort_values('clustering_seconds_median')
    fig,ax=plt.subplots(figsize=(9,4.7));style(ax);ax.grid(False)
    positions=np.arange(len(times));med=times.clustering_seconds_median.to_numpy()
    ax.barh(positions,med,color=BLUE,height=.62)
    ax.errorbar(med,positions,xerr=np.vstack([med-times.clustering_seconds_min,
        times.clustering_seconds_max-med]),fmt='none',ecolor='#292929',capsize=3,linewidth=1)
    ax.set_yticks(positions,[LABELS[x] for x in times.index]);ax.invert_yaxis()
    for y,v,maximum in zip(positions,med,times.clustering_seconds_max):
        ax.text(maximum+times.clustering_seconds_max.max()*.015,y,f'{v:.2f} s',va='center',fontsize=9)
    ax.set_xlim(0,times.clustering_seconds_max.max()*1.2)
    ax.set_xlabel('Combined component clustering workflow (seconds)')
    fig.suptitle('Clustering time for the pooled barcode components',x=.28,y=.985,ha='left',fontsize=13)
    fig.text(.28,.925,'Three serial repeats · bars: median · whiskers: min–max · preprocessing excluded',fontsize=8.5)
    fig.subplots_adjust(left=.28,right=.98,top=.855,bottom=.15)
    save(fig,'clustering_time')

    rows=[]
    for method in LABELS:
        r=summary.loc[method]
        rows.append(f'| {LABELS[method]} | {r.spearman_median:.6f} | {r.published_mass_percent_median:.2f}% | {r.clustering_seconds_median:.2f} | {int(r.unassigned_molecules):,} |')
    total=int(extraction.total_pairs.sum());retained=int(extraction.barbac_input_molecules.sum())
    r=summary.loc['barbac_lv'];h=summary.loc['barbac_hamming']
    report=f'''# Chen 2023 barcode time series

The complete hBFA1 / YPD subset is processed: **{total:,} read pairs across eight
samples**, at generations 8, 16, 24 and 40 in two biological assay replicates.
The shared mapping/extraction workflow retains **{retained:,} UMI-deduplicated
molecules ({100*retained/total:.2f}% of input pairs)** for the method comparison.
This percentage is extraction retention, not clustering accuracy.

For barbac LV with the Poisson option, median count-rank agreement with the
publication is **Spearman ρ = {r.spearman_median:.6f}** across the eight samples.
The median fraction of input molecules assigned to exact published pair IDs is
**{r.published_mass_percent_median:.2f}%**. Both barcode components together take
**{r.clustering_seconds_median:.2f} seconds** to cluster at the median of three
fresh serial repeats. Shared FASTQ preprocessing and trajectory reconstruction
are outside that clustering time.

| Method | Median Spearman ρ | Median molecule mass on published pair IDs | Clustering, median seconds | Unassigned molecules, all samples |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

The correlation uses all 2,314 published pair IDs at each sample, inserting zero
for absent method calls. Exact identifier matching is sensitive to different
centroid choices. Publication counts have additional filtering and correction
steps; they are a comparison reference, not known biological truth. These
measurements do not establish that one method has the highest true accuracy.
Starcode has slightly higher reference agreement; Hamming barbac has the lowest
median clustering time. LV barbac combines very high agreement with a similar
median time to Starcode.

The sample-median overlap above differs from the pooled molecule-weighted
calculation: **{diagnosis['pooled_unmatched_percent']:.2f}%** of all LV molecules
fall outside the published pair set. Two abundant unpublished pairs account for
**{diagnosis['largest_two_percent_of_unmatched']:.2f}% of this unmatched mass**.
Each has BC2 at least seven edits from every published BC2, outside the tested
distance-three correction radius. Thus the discrepancy includes abundant
identities absent from the reference. This does not establish whether those
pairs are biological lineages, artifacts, or removed by the author's filters.
The paper describes lane-intersection and chimera filtering; reproducing the
full author filtering and lineage-selection process is separate from this
clustering comparison. No diagnostic remapping was applied.

Shepherd could not automatically estimate an error rate for the low-diversity
BC1 component. Both components therefore use its documented `-e 0.005` option,
the rate already configured for barbac. This is a supplied assumption, not a
measured sequencing error rate or a value selected against publication counts.
The failed automatic attempt is retained in the generated audit files and is
excluded from the three successful timing repeats.

![Barcode trajectories](barcode_trajectories.png)

The trajectory comparison normalizes each method and the publication over the
same 2,314 published pair IDs. The table above exposes each method's molecule
mass outside that set rather than hiding it through normalization. L01–L03 are
chosen by pooled published abundance, before inspecting method agreement.

![LV barcode composition](barcode_composition.png)

The composition plot instead uses **all LV-assigned molecules**, including pairs
absent from the publication. Its four largest pairs are selected by pooled LV
counts. [figure_lineages.csv](figure_lineages.csv) identifies every labeled pair.

![Clustering workflow time](clustering_time.png)

## Reproducible outputs and interpretation

- `time_series_METHOD.csv.gz`: complete eight-sample tables for every observed
  method pair and every published pair, including explicit zero observations,
  counts, replicate, generation, method and both normalization denominators.
- `sample_method_agreement.csv`: all 48 sample/method comparisons, including
  unmatched and unassigned molecule counts, count-rank agreement and total
  variation between frequencies conditioned on the published pair set.
- `extraction_summary.csv`: reconciled quality, BAM extraction and UMI outcomes.
- `clustering_runs.csv`, `combined_clustering_times.csv`, `repeat_stability.csv`:
  all timing repeats, reported maximum child RSS, and membership stability.
- `published_identity_categories.csv`: molecule mass by exact component/pair
  membership in the published set, for every sample and method.
- `largest_unmatched_LV_pairs.csv`, `unmatched_diagnosis.json`: largest absent
  pairs and their nearest published component edit distances, without remapping.
- `provenance.json`: parameters, source hashes, tool hashes and input hashes.
- Every figure is also available as PDF and SVG for manuscript use.

The workflow is FastQC → PEAR → minimap2 → indexed BAM → `barbac_xtr()` flank
extraction → original-pair Q30 filtering and UMI deduplication. It uses the
existing `barbac_env`. Observed 24–28-base components remain variable length;
BC1 is reverse-complemented into the author's orientation. PEAR consensus can
change the original-mate barcode sequence; those changes are counted in the
extraction audit. Both components are clustered separately and their pair
identities preserved. All methods receive the same inputs at distance 3.

Clustering is retrospective: pooling all four generations uses information
from later samples. Frequency trajectories describe this selected fitness
assay; no fitness coefficients, causal effects or prospective performance are
claimed. Publication-specific GC, lane, chimera and abundance filters are not
re-created here, so differences may arise before or after clustering.

Source: [Chen, Johnson, Hérissant et al. (2023)](https://doi.org/10.7554/eLife.92899),
BioProject PRJNA912754, with author analysis files pinned to revision
`a375a116bf69160f634a3d6a1b0ecc53ad62d142`.
'''
    (OUT/'README.md').write_text(report)
    (OUT/'paper_time_series_section.md').write_text(f'''### Barcode time-series application

We analyzed eight barcode-sequencing samples from the hBFA1 YPD fitness assay
of Chen, Johnson, Hérissant et al. (2023), covering generations 8, 16, 24 and 40
in two biological assay replicates. Of {total:,} paired reads, the common
mapping, variable-length barcode extraction and UMI-deduplication workflow
retained {retained:,} molecules. We pooled each barcode component across the
eight samples for retrospective clustering at distance three and reconstructed
sample-specific paired-barcode counts using a consistent membership map.
Barbac LV with support ordering and the optional Poisson indel model achieved
median Spearman correlation {r.spearman_median:.6f} with published counts over
2,314 published pair identifiers; the corresponding Hamming value was
{h.spearman_median:.6f}. Median combined-component clustering workflow times
were {r.clustering_seconds_median:.2f} s and {h.clustering_seconds_median:.2f} s,
respectively, over three serial repeats. These are agreement measurements
against a differently processed experimental reference, not ground-truth
accuracy estimates. Starcode showed slightly higher reference agreement,
while Hamming barbac had the lowest median clustering time (Table 3).
The median proportion of LV molecules assigned to exact published pair IDs
was {r.published_mass_percent_median:.2f}%; two abundant absent pairs account for
{diagnosis['largest_two_percent_of_unmatched']:.2f}% of the pooled unmatched mass,
and both have BC2 at least seven edits from every published BC2. Their absence
cannot be attributed to clustering error from this comparison alone.

**Table 3. Publication agreement and clustering workflow time for the Chen 2023 time series.**
Spearman correlations use all 2,314 published pair IDs, filling absent calls with
zero. Molecule coverage is the percentage of method-input molecules assigned to
those exact IDs; both agreement columns are medians across eight samples.
Times are medians of three serial repeats, summing the BC2 and BC1 workflows
including startup and required exports, excluding shared FASTQ preprocessing.
Shepherd uses the supplied error rate 0.005 because automatic estimation failed
on BC1; the value matches the pre-existing barbac configuration and was not
selected using publication agreement. Unassigned totals span all eight samples.

| Method | Median Spearman ρ | Median molecule mass on published pair IDs | Clustering, median seconds | Unassigned molecules, all samples |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}
''')
    print('Rendered three figures in PNG/PDF/SVG and the results report.')


if __name__=='__main__':main()
