#!/usr/bin/env python3
"""Write the Phase 2S technical gate from audited available artifacts."""
from pathlib import Path
import hashlib
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'results/phase2S'
def rd(name):return pd.read_csv(P/name,sep='\t',dtype=str,keep_default_na=False)

def main():
    mf=rd('input_freeze_manifest.tsv');methods=rd('method_manifest.tsv');summary=rd('method_comparison_summary.tsv')
    genes=rd('frozen_gene_universes.tsv');seq=rd('sequence_manifest.tsv');partners=rd('frozen_partner_sets.tsv')
    original=ROOT/'reports/PHASE1B0R3_OFFLINE_MODEL_INGESTION_AND_SCORING.md'
    expected=mf.loc[mf.path.eq(str(original.relative_to(ROOT))),'sha256'].iloc[0]
    assert hashlib.sha256(original.read_bytes()).hexdigest()==expected
    assert '`NO_GO_EDGETIC_REPRESENTATION`' in original.read_text()
    assert len(genes)==132 and genes.U1_STRUCTURAL_42.eq('True').sum()==42
    acquired=ROOT/'logs/phase2S/encoding_model_acquisition.log'
    assert acquired.is_file() and 'SSLError' in acquired.read_text()
    text='# PHASE 2S — Representation Coverage Benchmark\n\n'
    text+='## Current gate\n\n`TECHNICAL_HOLD_SEQUENCE_MODELS`\n\n'
    text+='The frozen inputs, partner set, prespecified analysis and available structure-dependent baselines have been prepared. Valid sequence inference is blocked by official Hugging Face TLS failures for both required encoding models. No rescue outcome can be inferred from this technical hold. The original edgetic verdict remains `NO_GO_EDGETIC_REPRESENTATION`.\n\n'
    text+='## Frozen benchmark scope\n\n'
    text+=f'- U1: 42 exact Phase 1A structural-ready genes. U2: 132 of 137 strict genes after the five known-incompatible exclusions.\n- Frozen credible gene–partner pairs: {len(partners):,}.\n- Required distinct canonical accessions: {len(seq):,}; existing frozen sequences available: {int(seq.status.eq("AVAILABLE_FROZEN").sum()):,}; the remainder are explicitly pending sequence retrieval. No ambiguous accession or variant mapping was forced.\n- Every consumed frozen source and available cached sequence is listed with actual SHA256 in `metadata/PHASE2S_INPUT_FREEZE_MANIFEST.tsv` and its results copy.\n\n'
    text+='## Methods and technical validation (Table S1)\n\n'
    text+='SPPIDER-seq is pinned to the [manuscript release](https://github.com/aporollo-lab/SPPIDER-seq/tree/v2.0-manuscript), not the retrained v2.2.0 model. The [published article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142516/) identifies v2.0-manuscript as its archived release. X-PAIR is pinned to the [July 2026 PyPI release 0.1.0](https://pypi.org/project/xpair/0.1.0/), with its bundled pretrained checkpoint.\n\n'
    text+=methods[['method','version','source_revision','size_bytes','sha256','encoding_model','technical_validation']].to_markdown(index=False)+'\n\n'
    text+='The two SPPIDER neural-head weights and X-PAIR checkpoint were acquired and hashed. SPPIDER requires `facebook/esm2_t33_650M_UR50D`; X-PAIR checkpoint hyperparameters specify `ankh_large`, requiring `ElnaggarLab/ankh-large`. Both official model API requests failed with `SSLEOFError` in Python requests and SSL handshake errors in curl. No matching local Hugging Face cache was found. Raw acquisition evidence is in `logs/phase2S/encoding_model_acquisition.log`.\n\n'
    text+='SPPIDER_EXAMPLE_PASS and XPAIR_EXAMPLE_PASS are **not established**. Technical success rates are **NA (not run)**, and all sequence representation, rescue, consensus and structural-validation outcomes are **NA**. Empty raw-prediction tables indicate that inference was not run; they do not encode zero biological coverage.\n\n'
    text+='SPPIDER_RESCUE_STATUS = TECHNICAL_HOLD; XPAIR_RESCUE_STATUS = TECHNICAL_HOLD. Neither model has an evaluable FULL, PARTIAL or NO_RESCUE result.\n\n'
    text+='The reproducible environment targets Python 3.12.13, xpair 0.1.0, torch 2.7, pytorch-lightning 2.5.1, torchmetrics 1.7.1, esm 3.2.3 and transformers 4.48.1. An initial transformers 4.48.3 constraint conflicted with esm requirements; this constraint was corrected. The compatible retry then exited after a PyPI zstd request timed out; evidence is recorded in `logs/phase2S/environment_install_compatible.log`. This environment is incomplete. Installed package versions must be verified before any example or study inference.\n\n'
    text+='SPPIDER official code processes the complete sequence through overlapping ESM windows (1024 tokens including BOS/EOS; stride 512), merging query logits and pooling across partner chunks. Receptor-centric query probabilities are prespecified as primary; peptide-centric output is secondary. X-PAIR checkpoint length limits are 50–2000 residues for each input; unsupported inputs must be explicit failures, with no truncation. PPI-existence output must not filter interface probabilities.\n\n'
    text+='## Available coverage baseline (Table S2)\n\n'
    text+=summary.to_markdown(index=False)+'\n\n'
    text+='The original technical scoring success was 10,917/11,150 rows (97.91%), representing 10,499 unique inference units. Original ≥70% disruption representation covered 9 of the 42 genes, and the rebuilt primary cohort contained FGFR2, MEN1 and TUBA1A. Those historical results remain unchanged.\n\n'
    text+='The table above uses the frozen credible variant–partner universe as a common benchmark denominator. The original ≥70% gate is copied directly from the immutable legacy outputs. Audit found that legacy score-row membership and the full-cohort single-label registry produce different primary denominators. The initial derived cross-method C4/C5 results were withdrawn; they must not be used for rescue conclusions. New structural C4/C5 are NA pending exact membership reconciliation. Details are in `logs/phase2S/denominator_reconciliation.log`. A row lacking a structure remains unavailable for structural methods; this denominator differs from the already structure-selected 11,150-row production scoring denominator. Binary interface membership and finite interface distance have separate availability. Missing distances are not derived from binary labels. Matrices are saved as sparse long coordinates keyed by gene, variant and frozen partner.\n\n'
    text+='## Coverage QC and planned comparisons\n\n'
    text+='`analysis_specification.md` was saved before study-model predictions. The 70% gene threshold and strict disease-group balance thresholds are unchanged. The common sequence-arm partner universe includes all frozen credible PPI partners without requiring solved complexes, with evidence strata retained. Ordinary availability does not depend on a probability cutoff. Secondary CONSENSUS_HIGH requires both probabilities ≥0.5 and remains separate from consensus availability. No mutant-versus-WT deltas are permitted.\n\n'
    text+='Only paired descriptive comparisons are planned; no formal paired hypothesis test was selected. Structural validity would use experimental/PDB labels, protein/complex-grouped resampling, and frozen interface thresholds. Partner specificity and a small deterministic scrambled-partner control remain unevaluated until both examples pass. Technical failure-predictor covariates are preserved, but no most-important predictor has been estimated.\n\n'
    text+='## Deliverables and unavailable analyses\n\n'
    text+='| Table | File | Status |\n|---|---|---|\n| S1 Method/version manifest | `method_manifest.tsv` | Actual acquired checkpoint hashes; encoding models unavailable |\n| S2 Coverage by method | `method_comparison_summary.tsv` | Row/variant structural availability evaluated; new C4/C5 deferred; sequence arms NA |\n| S3 Coverage by gene | `coverage_by_gene.tsv` | Same |\n| S4 Disease-group balance | `coverage_by_disease_group.tsv` | Coverage QC only |\n| S5 Rescued genes | `sequence_rescued_genes.tsv` | Not evaluable |\n| S6 Technical failures | `technical_failures_by_method.tsv` | Model-level acquisition errors |\n| S7 PDB validation | `structural_validation.tsv` | Not run |\n\n'
    text+='All requested machine-readable output paths are present. Unexecuted sequence prediction outputs contain headers only, with method-level status explained in the manifest and failure table. Figures S1–S6 are deferred because the required cross-method sequence results are unavailable; no zero-valued rescue bars or fabricated validation curves were produced. This is a recoverable technical gate, not a completed coverage comparison.\n\n'
    text+='## Recovery steps\n\n'
    text+='Acquire the two exact encoding models from their official repositories, pin their revisions and verify actual file hashes. Complete and audit the pinned environment. Reconcile full-cohort primary membership against the exact frozen legacy rules before releasing any new C4/C5 endpoint. Obtain missing canonical sequences only for the frozen accession requests, saving new files under this workspace and recording their provenance. Run the official example gates first. Then infer each unique WT protein pair once, cache full-length outputs and extract frozen variant positions. Only after valid sequence inference may rescue counts, phenotype-balanced coverage, structural validation, specificity controls and the planned figures be calculated. Do not begin Phase 2P.\n\n'
    text+='Original NO_GO_EDGETIC_REPRESENTATION was not altered.\n\nDisease-label phenotype association testing was NOT started.\n'
    report=ROOT/'reports/PHASE2S_REPRESENTATION_COVERAGE_BENCHMARK.md'
    report.write_text(text,encoding='utf-8')
    print('report='+str(report),flush=True)

if __name__=='__main__':main()
