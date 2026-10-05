from pathlib import Path
import sys,json,hashlib,math,datetime
from collections import Counter
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'results/phase2P0';O.mkdir(exist_ok=True);M=O/'matrices';M.mkdir(exist_ok=True)
def read(p):return pd.read_csv(R/p,sep='\t',dtype=str,keep_default_na=False)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(rows,n,columns=None):pd.DataFrame(rows,columns=columns).to_csv(O/n,sep='\t',index=False,na_rep='NA')
def parts(s):return set(x for x in str(s).split(';') if x and x not in ['NA','nan'])
shared='COL3A1 COL4A1 COL4A3 FGFR2 GUCY2D KCNQ2 LHCGR MEN1 MFN2 MMUT PTEN RHO RPGR RS1 TUBA1A'.split();ext='COL7A1 DYNC1H1 MYO7A NF1 RYR1 SCN1A USH2A VWF'.split();genes=shared+ext
# Recheck all original frozen inputs and Phase2S reviewed outputs before construction.
freeze=read('metadata/PHASE2SR_INPUT_FREEZE_MANIFEST.tsv');audit=[]
for r in freeze.itertuples(index=False):
 observed=sha(R/r.path);audit.append(dict(path=r.path,SHA256=observed,status='PASS' if observed==r.SHA256 else 'FAIL'))
assert len(audit)==124 and all(x['status']=='PASS' for x in audit),'INPUT_FREEZE_VIOLATION'
write(audit,'source_freeze_validation.tsv')
gate=json.loads((R/'results/phase2S_full/final_gate.json').read_text());assert gate['figure_visual_QA']=='PASS' and gate['verdict']=='PARTIAL_SEQUENCE_RESCUE'
g=read('results/phase2S_full/rescued_gene_classification.tsv');assert set(g.loc[g.BOTH_RESCUED.eq('True'),'gene'])==set(shared);assert set(g.loc[g.SPPIDER_only_rescued.eq('True'),'gene'])==set(ext)
write([dict(gene=x,universe='BOTH_SEQUENCE_BALANCED_15') for x in shared],'frozen_shared15_genes.tsv');write([dict(gene=x,universe='SPPIDER_ONLY_EXTENSION_8') for x in ext],'sppider_extension8_genes.tsv')
v=read('results/phase2S/frozen_variants.tsv');v=v.loc[v.gene.isin(genes)].copy();groups=read('results/phase2S/frozen_disease_groups.tsv');allowed={x:set(z.loc[z.strict_disease_ge10.eq('YES'),'MONDO']) for x,z in groups.groupby('gene')};v['labels']=[parts(r.strict_MONDO_groups)&allowed.get(r.gene,set()) for r in v.itertuples()]
status={m:read('results/phase2S_full/pair_status_'+m+'.tsv') for m in ['sppider','xpair']};tracks={};common=[];common_by={};seqs={}
for gene in genes:
 good={m:status[m].loc[status[m].gene.eq(gene)&status[m].status.eq('SCORED')].set_index('partner') for m in status};pars=sorted(set(good['sppider'].index)&set(good['xpair'].index));common_by[gene]=pars
 for m in status:
  tracks[m,gene]={}
  selected=pars if gene in shared else sorted(good['sppider'].index) if m=='sppider' else []
  for par in selected:
   r=good[m].loc[par];path=R/r.score_array_path;assert sha(path)==r.score_array_SHA256
   z=np.load(path);a=z['query'] if m=='xpair' else z['receptor'];assert len(a)==int(r.query_length) and np.isfinite(a).all() and ((a>=0)&(a<=1)).all();tracks[m,gene][par]=a
  if len(good[m]):
   acc=good[m].iloc[0].query_accession
   if acc not in seqs:seqs[acc]=''.join(s.strip() for s in (R/'external_models/canonical_sequences'/(acc+'.fasta')).read_text().splitlines() if not s.startswith('>'))
 for par in pars:
  a=good['xpair'].loc[par];b=good['sppider'].loc[par];assert a.query_sha256==b.query_sha256 and a.partner_sha256==b.partner_sha256
  common.append(dict(gene=gene,partner=par,query_accession=a.query_accession,partner_accession=a.partner_accession,xpair_track=a.score_array_path,xpair_SHA256=a.score_array_SHA256,sppider_track=b.score_array_path,sppider_SHA256=b.score_array_SHA256,selection_rule='BOTH_VALID_FULL_WT_TRACKS_NO_DISEASE_LABELS',scope='PRIMARY' if gene in shared else 'SECONDARY_EXTENSION_COMMON_AUDIT'))
write(common,'common_partner_set_by_gene.tsv')
registry=[];excluded=[];variant_rows={};counter=0
for row in v.itertuples():
 if not row.canonical_residue.isdigit():excluded.append(dict(gene=row.gene,VariationID=row.VariationID,reason='CANONICAL_RESIDUE_UNAVAILABLE'));continue
 variant_rows.setdefault((row.gene,row.canonical_UniProt_accession,int(row.canonical_residue)),[]).append(row)
for (gene,acc,pos),rs in sorted(variant_rows.items()):
 q=seqs.get(acc,'');safe=all(r.mapping_status in ['EXACT_CANONICAL','SAFE_ISOFORM_REMAP'] and r.primary_cohort_eligible=='YES' for r in rs) and 1<=pos<=len(q) and all(q[pos-1]==r.canonical_reference_aa for r in rs)
 labels=set().union(*(r.labels for r in rs));multi=any(len(r.labels)>1 for r in rs);ambig=any(not r.labels for r in rs)
 labelstatus='MULTILABEL_VARIANT_RESIDUE' if multi else 'AMBIGUOUS_DISEASE_RESIDUE' if ambig else 'CROSS_DISEASE_RESIDUE' if len(labels)>1 else 'SINGLE_DISEASE_RESIDUE' if len(labels)==1 else 'AMBIGUOUS_DISEASE_RESIDUE'
 sp=tracks['sppider',gene];xp=tracks['xpair',gene];count=len(common_by[gene]);rid=f'{gene}|{acc}|{pos}'
 registry.append(dict(residue_id=rid,gene=gene,universe='BOTH_SEQUENCE_BALANCED_15' if gene in shared else 'SPPIDER_ONLY_EXTENSION_8',canonical_accession=acc,canonical_residue=pos,WT_amino_acid=q[pos-1] if 1<=pos<=len(q) else 'NA',VariationIDs=';'.join(sorted({r.VariationID for r in rs})),substitutions=';'.join(sorted({r.protein_change for r in rs})),MONDO_labels=';'.join(sorted(labels)),n_substitutions=len({r.protein_change for r in rs}),n_Variation_IDs=len({r.VariationID for r in rs}),n_disease_labels=len(labels),ClinVar_RCVs=';'.join(sorted(set().union(*(parts(r.RCVs) for r in rs)))),phenotype_status=labelstatus,cross_disease_linked=len(labels)>1,multilabel_variant=multi,ambiguous_label=ambig,canonical_mapping_safe=safe,XPAIR_available=bool(xp) and safe,SPPIDER_available=bool(sp) and safe,common_partner_n=count,model_membership_eligible=safe and labelstatus=='SINGLE_DISEASE_RESIDUE' and (count>=3 if gene in shared else len(sp)>=3),protein_length=len(q)))
reg=pd.DataFrame(registry);write(reg,'residue_phenotype_registry.tsv');cross=[];ceil=[];groupcounts=[];ready=[];primary={}
for gene in genes:
 rr=reg.loc[reg.gene.eq(gene)];annot=rr.loc[rr.n_disease_labels.gt(0)];dis=annot.loc[annot.cross_disease_linked];combinations=Counter(dis.MONDO_labels);cross.append(dict(gene=gene,phenotype_annotated_residues=len(annot),cross_disease_linked_residues=len(dis),DiscordantResidueRate=len(dis)/len(annot) if len(annot) else np.nan,DiscordantResiduePct=100*len(dis)/len(annot) if len(annot) else np.nan,disease_combinations=json.dumps(dict(combinations),sort_keys=True)))
 mass_total=0.;correct=0.;macro=[]
 for row in rr.itertuples():
  rs=variant_rows[row.gene,row.canonical_accession,row.canonical_residue];mass=Counter()
  for r in rs:
   if r.labels:
    for label in r.labels:mass[label]+=1/len(r.labels)
  n=sum(mass.values())
  if n:mass_total+=n;correct+=max(mass.values());macro.append(max(mass.values())/n)
 ceil.append(dict(gene=gene,RESIDUE_INFORMATION_CEILING=correct/mass_total if mass_total else np.nan,variant_weighted_label_mass=mass_total,maximum_correct_mass=correct,macro_residue_ceiling=np.mean(macro) if macro else np.nan,definition='ONE_VOTE_PER_VARIANT_FRACTIONAL_ACROSS_MULTILABELS_MAJORITY_PER_RESIDUE',score_optimized=False))
 sub=rr.loc[rr.model_membership_eligible];counts=sub.MONDO_labels.value_counts();retained=sorted(counts.index[counts.ge(10)]);all_single=rr.loc[rr.phenotype_status.eq('SINGLE_DISEASE_RESIDUE')];allcounts=all_single.MONDO_labels.value_counts()
 for label in sorted(allowed.get(gene,set())):groupcounts.append(dict(gene=gene,MONDO=label,unique_single_disease_residues=int(allcounts.get(label,0)),mapped_joint_representation_residues=int(counts.get(label,0)),retained_for_testing=int(counts.get(label,0))>=10))
 eligible=len(retained)>=2 and (len(common_by[gene])>=3 if gene in shared else len(tracks['sppider',gene])>=3);prim=sub.loc[sub.MONDO_labels.isin(retained)] if eligible else sub.iloc[:0];primary[gene]=prim
 reasons=[]
 if len(counts)<2:reasons.append('FEWER_THAN_TWO_SINGLE_DISEASE_GROUPS_WITH_REPRESENTATION')
 if len(retained)<2:reasons.append('FEWER_THAN_TWO_DISEASE_GROUPS_WITH_GE10_UNIQUE_RESIDUES')
 if gene in shared and len(common_by[gene])<3:reasons.append('FEWER_THAN_THREE_COMMON_PARTNERS')
 ready.append(dict(gene=gene,universe='PRIMARY_SHARED15' if gene in shared else 'SECONDARY_SPPIDER_EXTENSION8',single_disease_groups=len(counts),retained_ge10_groups=len(retained),retained_MONDOs=';'.join(retained),common_partner_n=len(common_by[gene]),PHENOTYPE_READY=eligible,primary_residue_n=len(prim),exclusion_reason=';'.join(reasons) if reasons else 'PASS'))
 for row in rr.itertuples():
  reason='PRIMARY_RETAINED' if row.residue_id in set(prim.residue_id) else row.phenotype_status if row.phenotype_status!='SINGLE_DISEASE_RESIDUE' else 'UNSAFE_CANONICAL_MAPPING' if not row.canonical_mapping_safe else 'INSUFFICIENT_COMMON_REPRESENTATION' if not row.model_membership_eligible else 'DISEASE_GROUP_OR_GENE_SAMPLE_SIZE_GATE'
  excluded.append(dict(gene=gene,residue_id=row.residue_id,VariationID=row.VariationIDs,reason=reason))
write(cross,'cross_disease_residue_audit.tsv');write(ceil,'residue_information_ceiling.tsv');write(groupcounts,'phenotype_group_residue_counts.tsv');write(ready,'phenotype_ready_genes.tsv');write(excluded,'phenotype_exclusion_log.tsv')
# Feature files contain no association computation. All rows use frozen residue order.
manifest={k:[] for k in ['xpair_primary','sppider_confirmatory','xpair_unlabeled_control','sppider_unlabeled_control','baseline_position','baseline_domain','baseline_generic_interface','baseline_generic3d']};domain=read('results/phase1B0/baseline_domain.tsv');geo=read('results/phase1B0/baseline_generic3d.tsv');struct=read('results/phase1A/partner_specific_matrix_structural.tsv');struct=struct.loc[struct.gene.isin(genes)];baselines={}
for gene in shared:
 rr=primary[gene].sort_values('canonical_residue');partners=common_by[gene]
 if rr.empty:
  for key in manifest:manifest[key].append(dict(gene=gene,status='GENE_NOT_PHENOTYPE_READY',rows=0,columns=0,path='NA',SHA256='NA'))
  continue
 ids=rr.residue_id.to_numpy(dtype=str);positions=rr.canonical_residue.to_numpy(dtype=int);membership=M/(gene+'_residue_membership.tsv');rr[['residue_id','gene','canonical_accession','canonical_residue','MONDO_labels']].to_csv(membership,sep='\t',index=False)
 for m,key,control in [('xpair','xpair_primary','xpair_unlabeled_control'),('sppider','sppider_confirmatory','sppider_unlabeled_control')]:
  a=np.column_stack([tracks[m,gene][p][positions-1] for p in partners]);sorted_a=np.sort(a,axis=1)[:,::-1]
  for which,data,cols in [(key,a,partners),(control,sorted_a,[f'rank_{i+1}' for i in range(len(partners))])]:
   path=M/(gene+'_'+which+'.npz');np.savez(path,scores=data,residue_ids=ids,columns=np.array(cols,dtype=str));manifest[which].append(dict(gene=gene,status='FROZEN_RAW_PROBABILITIES',rows=len(rr),columns=len(cols),path=str(path.relative_to(R)),SHA256=sha(path),membership_path=str(membership.relative_to(R)),membership_SHA256=sha(membership),nonconstant_dimensions=int((data.std(axis=0,ddof=0)>0).sum()),disease_labels_used_to_transform=False))
  sums=pd.DataFrame(dict(residue_id=ids,mean=a.mean(axis=1),maximum=a.max(axis=1),SD=a.std(axis=1,ddof=0),fraction_ge0_5=(a>=.5).mean(axis=1)));sums.to_csv(M/(gene+'_'+m+'_collapsed_controls.tsv'),sep='\t',index=False)
 for key in ['baseline_position','baseline_domain','baseline_generic_interface','baseline_generic3d']:
  rows=[]
  for row in rr.itertuples():
   d=dict(residue_id=row.residue_id,gene=gene,canonical_residue=row.canonical_residue)
   if key=='baseline_position':d.update(residue_number=row.canonical_residue,normalized_position=row.canonical_residue/row.protein_length,status='AVAILABLE')
   elif key in ['baseline_domain','baseline_generic3d']:
    src=domain if key=='baseline_domain' else geo;ss=src.loc[src.gene.eq(gene)&src.canonical_residue.eq(str(row.canonical_residue))];cols=['domain_identity','relative_position_within_domain','distance_to_nearest_domain_boundary_aa','annotation_source','annotation_status'] if key=='baseline_domain' else ['AFDB_plddt','monomer_SASA_A2','intrachain_CA_neighbor_n_8A','generic3d_status']
    if len(ss):
     vals=ss[cols].drop_duplicates();assert len(vals)==1,'INCONSISTENT_FROZEN_RESIDUE_BASELINE';d.update(vals.iloc[0].to_dict());d['status']='FROZEN_AVAILABLE'
    else:d.update({c:'NA' for c in cols});d['status']='FROZEN_ANNOTATION_UNAVAILABLE'
   else:
    ss=struct.loc[struct.gene.eq(gene)&struct.canonical_residue.eq(str(row.canonical_residue))];known=ss.loc[ss.interface_status.isin(['INTERFACE','NON_INTERFACE_IN_OBSERVED_COMPLEX','NON_INTERFACE_IN_PIONEER_INTERFACE_MODEL'])];distance=pd.to_numeric(ss.interface_min_distance_A,errors='coerce');d.update(any_interface=bool(known.interface_status.eq('INTERFACE').any()) if len(known) else 'NA',known_interfaces_touching_residue=known.loc[known.interface_status.eq('INTERFACE'),'partner'].nunique() if len(known) else 'NA',minimum_interface_distance=distance.min() if distance.notna().any() else 'NA',status='FROZEN_AVAILABLE' if len(known) or distance.notna().any() else 'STRUCTURAL_INFORMATION_MISSING')
   rows.append(d)
  path=M/(gene+'_'+key+'.tsv');pd.DataFrame(rows).to_csv(path,sep='\t',index=False,na_rep='NA');manifest[key].append(dict(gene=gene,status='FROZEN_WITH_EXPLICIT_MISSINGNESS',rows=len(rows),columns=len(rows[0]),path=str(path.relative_to(R)),SHA256=sha(path),source='FROZEN_UNIPROTKB_FEATURE_PROXY_NOT_PFAM_INTERPRO' if key=='baseline_domain' else 'FROZEN_PHASE1_OUTPUT'))
for key,rows in manifest.items():write(rows,key+'_matrix_manifest.tsv' if key in ['xpair_primary','sppider_confirmatory'] else key+'_manifest.tsv')
# Audit query and pair scopes separately: public identity lists unavailable, never infer no overlap.
overlap=[]
for gene in genes:
 for method in ['XPAIR_PARTNER_INTERFACE','SPPIDER_SEQ_PARTNER_INTERFACE']:
  for scope in ['QUERY_PROTEIN','QUERY_PARTNER_PAIRS']:
   overlap.append(dict(gene=gene,method=method,scope=scope,status='TRAINING_OVERLAP_UNKNOWN',query_accession=reg.loc[reg.gene.eq(gene),'canonical_accession'].iloc[0],common_pairs_audited=len(common_by[gene]),reason='BUNDLED_PACKAGE_HAS_NO_COMPLETE_TRAINING_IDENTITY_LIST;OFFICIAL_RESOURCE_FETCH_503;CURRENT_SPPIDER_DATASETS_NOT_PROVEN_EXACT_V2_TRAINING_SPLIT',official_resource='https://gitlab.lcqb.upmc.fr/srescalli/X-PAIR' if method.startswith('XPAIR') else 'https://huggingface.co/datasets/aporollo-lab/SPPIDER-seq-datasets'))
write(overlap,'training_overlap_audit.tsv')
permutations=[]
for gene in shared:
 rr=primary[gene];sizes=rr.MONDO_labels.value_counts().tolist();logn=math.lgamma(len(rr)+1)-sum(math.lgamma(n+1) for n in sizes) if len(rr) else 0;space=math.factorial(len(rr))//math.prod(math.factorial(n) for n in sizes) if len(rr) else 0;seed=int(hashlib.sha256(('PHASE2P1|20261004|'+gene).encode()).hexdigest()[:16],16)%(2**32)
 permutations.append(dict(gene=gene,status='FROZEN_NO_PERMUTATIONS_EXECUTED' if len(rr) else 'GENE_INELIGIBLE',seed=seed,requested_B=10000,label_assignment_space=str(space),permutation_B=min(10000,space-1) if space else 0,group_sizes=json.dumps(sizes),group_size_preservation='EXACT',paired_across_labeled_sorted_and_methods=True,restricted_min_valid=1000))
write(permutations,'permutation_plan.tsv')
readydf=pd.DataFrame(ready);pr=readydf.loc[readydf.universe.eq('PRIMARY_SHARED15')&readydf.PHENOTYPE_READY];n=len(pr);verdict='GO_PHASE2P1_FOCUSED' if n>=10 else 'HOLD_LOW_GENE_COUNT' if n>=5 else 'NO_GO_PHENOTYPE_TEST';sharedreg=reg.loc[reg.gene.isin(shared)];crossdf=pd.DataFrame(cross);summary=dict(verdict=verdict,shared_genes=15,extension_genes=8,total_unique_residues=len(sharedreg),single_disease_residues=int(sharedreg.phenotype_status.eq('SINGLE_DISEASE_RESIDUE').sum()),cross_disease_residues=int(sharedreg.phenotype_status.eq('CROSS_DISEASE_RESIDUE').sum()),multilabel_ambiguous_residues=int(sharedreg.phenotype_status.isin(['MULTILABEL_VARIANT_RESIDUE','AMBIGUOUS_DISEASE_RESIDUE']).sum()),median_discordant_residue_rate=float(crossdf.loc[crossdf.gene.isin(shared),'DiscordantResidueRate'].median()),genes_ge2_groups=int(readydf.loc[readydf.universe.eq('PRIMARY_SHARED15'),'single_disease_groups'].ge(2).sum()),genes_ge2_ge10_groups=int(readydf.loc[readydf.universe.eq('PRIMARY_SHARED15'),'retained_ge10_groups'].ge(2).sum()),genes_ge3_common_partners=sum(len(common_by[g])>=3 for g in shared),phenotype_ready_n=n,phenotype_ready_genes=pr.gene.tolist(),total_primary_residues=int(pr.primary_residue_n.sum()),median_residues_per_disease=float(pd.DataFrame(groupcounts).loc[lambda d:d.gene.isin(pr.gene)&d.retained_for_testing,'mapped_joint_representation_residues'].median()),median_common_partners=float(pr.common_partner_n.median()),known_training_overlap_genes=0,unknown_training_overlap_genes=15,association_statistics_computed=False,phase2P1_started=False,original_verdict='NO_GO_EDGETIC_REPRESENTATION')
(O/'cohort_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
