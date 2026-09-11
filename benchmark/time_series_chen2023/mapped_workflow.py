"""Run or verify the mapping and BAM-extraction stages for one selected sample."""
import hashlib
import json
from pathlib import Path
import subprocess

from extract_barcodes import sha256
from extract_mapped_pairs import join_pairs
from run_sample import checked, digest

HERE=Path(__file__).resolve().parent


def run(stage, work, index, row):
    work=Path(work).resolve()
    raw=work/'raw'/row['run']
    files=[raw/f'{row["run"]}_{i}.fastq.gz' for i in (1,2)]
    if not all(checked(p,row,i) for i,p in enumerate(files,1)):
        raise ValueError('Full input FASTQs must pass manifest size/MD5 checks')
    mapping=work/'mapping'/row['sample']
    if stage=='map' and not (mapping/'mapping.json').exists():
        subprocess.run(['Rscript',str(HERE/'map_sample.R'),str(work),str(index)],check=True)
    receipt=json.loads((mapping/'mapping.json').read_text())
    if receipt['status']!='complete': raise ValueError('Mapping is not complete')
    reference=HERE/'reference/chen2023_masked_amplicon.fasta'
    if 'reference_md5' in receipt:
        if receipt['reference_md5']!=digest(reference): raise ValueError('Reference changed since mapping')
        if receipt['input_md5']!=[digest(p) for p in files]: raise ValueError('Mapped input hashes differ')
    elif receipt.get('reference_sha256')!=sha256(reference):
        raise ValueError('Cached reference hash differs')
    bam=mapping/'merged/bam'/f'{row["run"]}_ANC.assembled_sorted.bam'
    if not bam.is_file() or not Path(str(bam)+'.bai').is_file(): raise ValueError('BAM or index missing')
    if stage=='map': return receipt
    out=work/'extracted'/row['sample']
    if (out/'extraction.json').exists():
        saved=json.loads((out/'extraction.json').read_text())
        if saved['bam_extraction']['bam_md5']!=digest(bam):
            raise ValueError('BAM changed since completed paired extraction')
        if list(saved['inputs'].values())!=[sha256(p) for p in files]:
            raise ValueError('FASTQs changed since completed paired extraction')
        for name,expected in saved['outputs'].items():
            if sha256(out/name)!=expected: raise ValueError('Completed extraction output changed')
        return saved
    components=work/'bam_components'/row['sample']
    if not (components/'bam_components.json').exists():
        subprocess.run(['Rscript',str(HERE/'extract_mapped_components.R'),str(bam),str(components)],check=True)
    extraction=json.loads((components/'bam_components.json').read_text())
    if extraction['bam_md5']!=digest(bam): raise ValueError('BAM changed since component extraction')
    for name,expected in extraction['output_md5'].items():
        if digest(components/Path(name).name)!=expected: raise ValueError('BAM component CSV changed')
    return join_pairs(*files,components,out,row)
