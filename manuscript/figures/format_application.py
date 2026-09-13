"""Repair a clipped axis label in the saved vector plot; preserve every data path."""
import hashlib,json,pathlib
from pypdf import PdfReader,PdfWriter
from pypdf.generic import DecodedStreamObject,NameObject
HERE=pathlib.Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'benchmark/workflow_runtime/generated/run-01/lineage_trajectories.pdf'
r=PdfReader(SOURCE);page=r.pages[0];data=page.get_contents().get_data()
old=b'0.00 20.00 -20.00 0.00 21.33 -12.82 Tm [(Fr) 10 (action of e) 30 (xtr) 10 (acted barcode reads)] TJ'
new=b'0.00 16.00 -16.00 0.00 20.00 79.50 Tm (Barcode frequency) Tj'
assert data.count(old)==1
updated=data.replace(old,new)
assert updated.replace(new,old)==data
stream=DecodedStreamObject();stream.set_data(updated)
page[NameObject('/Contents')]=stream.flate_encode()
page.mediabox.upper_right=(510,259)
w=PdfWriter();w.add_page(page)
with (HERE/'workflow_application.pdf').open('wb') as f:w.write(f)
(HERE/'application_formatting.json').write_text(json.dumps(dict(
  source_pdf_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
  formatted_pdf_sha256=hashlib.sha256((HERE/'workflow_application.pdf').read_bytes()).hexdigest(),
  data_paths_unchanged=True,changes=['Shorter centred y-axis label','Six-point right margin for last tick label'],
  timed_output_preserved=True,formatting_excluded_from_workflow_timing=True),indent=2)+'\n')
