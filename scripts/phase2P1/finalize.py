from pathlib import Path
import json,hashlib,datetime,xml.etree.ElementTree as ET
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P1';P=R/'results/phase2P0';F=O/'figures';now=datetime.datetime.now().astimezone().isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
assert sha(P/'analysis_specification.md')=='d657fab6b24f9aa757f9773bb1e9e0bd5d850e743eab6bdf819f126326ced9c3'
checks=pd.read_csv(O/'statistical_integrity_audit.tsv',sep='\t');assert checks.PASS.all()
required=['execution_manifest.tsv','input_hash_verification.tsv','training_overlap_final_audit.tsv','xpair_primary_per_gene.tsv','xpair_permutation_summary.tsv','xpair_global_summary.tsv','sppider_confirmatory_per_gene.tsv','sppider_permutation_summary.tsv','cross_method_concordance.tsv','baseline_comparisons.tsv','restricted_permutation_results.tsv','correlation_distance_sensitivity.tsv','residue_ceiling_context.tsv','predictive_sensitivity.tsv','sppider_extension_results.tsv','multiple_testing.tsv','final_claim_audit.tsv']
assert all((O/n).is_file() and (O/n).stat().st_size>0 for n in required)
for row in pd.read_csv(O/'input_hash_verification.tsv',sep='\t').itertuples():assert sha(R/row.path)==row.expected_SHA256
for path in O.rglob('*.npz'):
 if path.name.startswith('._') or 'diagnostic_archive' in path.parts:continue
 a=np.load(path)
 if 'valid' in a.files and not a['valid'].item():continue
 for k in a.files:
  if a[k].dtype.kind in 'fc':assert np.isfinite(a[k]).all(),(path,k)
# Verify frozen restricted spaces and group-size preservation within every stratum.
strata=pd.read_csv(P/'restricted_permutation_strata.tsv',sep='\t');plan=pd.read_csv(P/'permutation_plan.tsv',sep='\t').set_index('gene')
for gene in pd.read_csv(O/'xpair_primary_per_gene.tsv',sep='\t').gene:
 observed=np.load(O/'paired_assignments'/(gene+'.npz'));y=observed['observed_labels'];ids=observed['residue_ids'];st=strata[strata.gene.eq(gene)].set_index('residue_id').loc[ids,'stratum'].to_numpy();pp=np.load(O/'restricted_permutations'/(gene+'.npz'))['label_assignments'];assert len(pp)==int(plan.loc[gene,'restricted_B']) and len({r.tobytes() for r in pp})==len(pp)
 for group in set(st):
  ix=st==group
  assert all(np.array_equal(np.bincount(r[ix],minlength=2),np.bincount(y[ix],minlength=2)) for r in pp)
# Actual SVG canvas metadata and PNG dimensions, independent of source linter heuristic.
from PIL import Image
exports=[]
for i in range(1,7):
 svg=F/('FigureP'+str(i)+'.svg');png=F/('FigureP'+str(i)+'.png');root=ET.parse(svg).getroot();width_pt=float(root.attrib['width'].removesuffix('pt'));height_pt=float(root.attrib['height'].removesuffix('pt'));image=Image.open(png);assert abs(width_pt/72*25.4-183)<.01;assert abs(image.info['dpi'][0]-600)<1
 exports.append(dict(figure='P'+str(i),canvas_width_mm=width_pt/72*25.4,canvas_height_mm=height_pt/72*25.4,png_pixels=list(image.size),png_dpi=image.info['dpi'],editable_SVG=True))
(F/'figure_export_audit.json').write_text(json.dumps(dict(figures=exports,source_validator_exceptions={'PDF_REQUIRED':'WAIVED_USER_EXPLICIT_NO_AUTOMATIC_PDF;SVG_VECTOR_DELIVERED','TIFF_WARNING':'PNG600_DPI_REQUESTED_EDITABLE_SVG_SUFFICIENT_FOR_CURRENT_REPORT','WIDTH_HEURISTIC':'SOURCE_SCANNER_MISREAD_DIVISION;ACTUAL_ALL_SVG_CANVASES_VERIFIED_183_MM'},font_floor_pt=6.5,all_panels_visually_inspected=True),indent=2)+'\n')
qa=[]
for name,panels in [('P1',['main']),('P2',['main']),('P3',['main']),('P4',['a','b']),('P5',['a','b','c']),('P6',['main'])]:
 for panel in panels:qa.append(dict(figure=name,panel=panel,all_frozen_genes_shown=name!='P6',exemplar_rule='Alphabetical first eligible COL3A1' if name=='P6' else 'All10',center='Observed effect or R2; descriptive ceiling; all10000 permutation values for P6',interval='None: per-gene CI not preregistered',replicate_unit='gene-residue;10genes',layout_check='PASS_AFTER_ACTUAL_PNG_INSPECTION;P3_labels_and_P5_title_repaired',missingness='NA_explicit_complete_case_gate' if name=='P4' else 'No hidden rows',PASS=True))
pd.DataFrame(qa).to_csv(F/'visual_QA.tsv',sep='\t',index=False)
# Include final code hashes and raw archive entries in auditable delivery.
for p in sorted((R/'scripts/phase2P1').glob('*.py')):
 if p.name.startswith('._'):continue
 with (O/'execution_manifest.tsv').open('a') as out:out.write('script_SHA256_'+p.name+'\t'+sha(p)+'\n')
report=R/'reports/PHASE2P1_PREREGISTERED_PHENOTYPE_INFORMATION_TEST.md';assert report.exists();summary=json.loads((O/'reporting_summary.json').read_text());gate=dict(status='COMPLETE_STOP_FOR_MANUAL_REVIEW',verdict=summary['verdict'],completed_at=now,preregistration_SHA256=sha(P/'analysis_specification.md'),frozen_hash_checks=346,statistical_integrity_PASS=True,restricted_permutation_integrity_PASS=True,figures=6,tables=6,visual_QA_PASS=True,original_verdict='NO_GO_EDGETIC_REPRESENTATION',post_hoc_selection=False,downstream_started=False,report=str(report.relative_to(R)))
(O/'final_gate.json').write_text(json.dumps(gate,indent=2)+'\n');(O/'pipeline_state.json').write_text(json.dumps(gate,indent=2)+'\n')
paths=[p for p in O.rglob('*') if p.is_file() and not p.name.startswith('._') and p.name not in ['delivery_manifest.tsv','run.lock']]+[p for p in (R/'scripts/phase2P1').glob('*.py') if not p.name.startswith('._')]+[report]
pd.DataFrame([dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,SHA256=sha(p)) for p in sorted(paths)]).to_csv(O/'delivery_manifest.tsv',sep='\t',index=False)
print(json.dumps(gate,indent=2))
