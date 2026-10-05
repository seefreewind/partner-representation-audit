# Phase 0R — Full ClinVar RCV recovery and 965-gene re-audit

**PHASE 0R VERDICT: `GO_CONDITIONAL_MECHANISM`**

This report uses the frozen September 2026 ClinVar snapshot and a newly reconstructed 965-gene prescreen list. All counts below come from the Phase 0R workflow.

## Frozen release and seed checkpoint

- Frozen ClinVar monthly snapshot: 2026-09-03 (`variant_summary_2026-09.txt.gz`; archived checksum recorded in the checkpoint).
- Frozen RCV file: `ClinVarRCVRelease_2026-09.xml.gz`; official archive path `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_2026-09.xml.gz`, server last-modified 2026-09-03; the XML root identifies release date `2026-08-29`. Both dates are recorded because they describe the archive path and its internal release date. Official sidecar MD5 and gzip integrity checks are required before parsing.
- Global seed: 61,766 unique Variation IDs; 60,456 unique gene–protein substitutions; 114,506 candidate RCV accessions; 4,176 genes.
- Recovery: 114,506/114,506 RCVs (100.00%); 0 failed; recovery gate `CONTINUE`.
- The archive was resumed from the existing partial file, validated with `gzip -t` and official MD5 (`True`), then streamed with `lxml.etree.iterparse`; only seed-listed RCVs were retained. Full-archive SHA256/MD5 and recovery source for every RCV are in the manifests/table.

## Rebuilt 965-gene attrition

The preserved preliminary disease rows reproduced exactly 965 genes under the original summary-screen rule (≥2 MONDO groups each with ≥10 summary-screen substitutions). These are prescreen candidates, not exact-RCV findings.

| Stage | Genes |
|---|---:|
| previous strict pre-screen genes | 965 |
| after exact RCV P/LP missense | 941 |
| after unambiguous MONDO normalization | 763 |
| at least 2 diseases after parent-child collapse | 549 |
| at least 2 diseases with >=5 unique variants each | 253 |
| strict >=10 variants/disease and >=20 total | 137 |

- After exact RCV P/LP germline missense review: 941 genes.
- After at least one unambiguous MONDO disease mapping: 763 genes; with ≥2 disease groups after same-MONDO and parent-child collapsing: 549.
- Expanded (≥2 diseases, ≥5 unique substitutions per disease): 253 genes.
- Strict (≥2 diseases, ≥10 substitutions per disease, ≥20 total): 137 genes.
- Unambiguous parent terms were removed from independent disease counts when a more specific descendant was also present; variants attached only to a removed broad parent were not reassigned to a child. Relationship details are in `disease_parent_child_collapses.tsv`.
- Exact synonyms and duplicate variant–condition submissions are collapsed. Hotspot flag: at least one strict-sized disease group has a single residue with ≥20% of substitutions or its top three residues with ≥50%. Residue counts remain separately reported.

| Minimum substitutions per disease | Genes | Disease groups | Unique substitutions |
|---:|---:|---:|---:|
| 2 | 410 | 1340 | 17996 |
| 3 | 342 | 1032 | 17119 |
| 5 | 253 | 701 | 15792 |
| 10 | 155 | 390 | 13054 |
| 15 | 105 | 239 | 10934 |
| 20 | 76 | 174 | 9120 |

## Mechanism states and PPI coverage

- Strict exact-RCV genes: 137. Mechanism gene states: `CONFIRMED_COMPATIBLE` 0; `UNRESOLVED` 132; `CONFIRMED_INCOMPATIBLE` 5; `NOT_KNOWN_INCOMPATIBLE` 132. Counts are genes; disease-pair states are separately listed. Missing G2P/GenCC evidence is `UNRESOLVED`.
- Strict genes with ≥3 credible PPI partners: 133 of 137. Credibility includes experimental BioGRID/IntAct support or precomputed very-high-confidence PIONEER / Predictomes SPOC >0.86 coverage; Tier 1/2 experimental counts remain separate. HuRI was unavailable after a timeout.
- Partner-analysis-ready coverage proxy: 96.4% of strict-group variants are in genes with ≥3 credible partners; this is a gene-level PPI proxy, not variant-level interface readiness.
- Final high-quality analyzable genes: 0. A gene qualifies here only if it passes strict disease thresholds, is `CONFIRMED_COMPATIBLE`, and has ≥3 credible partners. Tier 1/2 experimental partner counts are reported separately.

## Top exact-RCV candidates

Ranked by mechanism state, disease distinction, variant counts, then credible/PIONEER partner coverage; no predicted perturbation outcome was used.

| Rank | Gene | Distinct diseases | Strict-group variants | Mechanism | Credible partners | PIONEER partners | Analyzable |
|---:|---|---:|---:|---|---:|---:|---|
| 1 | LMNA | 23 | 162 | UNRESOLVED | 803 | 444 | NO |
| 2 | TRPV4 | 17 | 40 | UNRESOLVED | 17 | 8 | NO |
| 3 | HBB | 17 | 54 | UNRESOLVED | 119 | 7 | NO |
| 4 | MYH7 | 15 | 165 | UNRESOLVED | 72 | 19 | NO |
| 5 | FBN1 | 15 | 750 | UNRESOLVED | 32 | 16 | NO |
| 6 | PTPN11 | 14 | 84 | UNRESOLVED | 466 | 89 | NO |
| 7 | SCN2A | 14 | 70 | UNRESOLVED | 22 | 7 | NO |
| 8 | ATM | 13 | 83 | UNRESOLVED | 226 | 83 | NO |
| 9 | SCN4A | 13 | 83 | UNRESOLVED | 170 | 4 | NO |
| 10 | FGFR2 | 12 | 42 | UNRESOLVED | 191 | 37 | NO |
| 11 | FGFR3 | 11 | 30 | UNRESOLVED | 363 | 324 | NO |
| 12 | ACTA1 | 11 | 133 | UNRESOLVED | 108 | 31 | NO |
| 13 | CACNA1A | 11 | 53 | UNRESOLVED | 166 | 88 | NO |
| 14 | FKRP | 11 | 56 | UNRESOLVED | 32 | 1 | NO |
| 15 | MYO7A | 10 | 83 | UNRESOLVED | 11 | 4 | NO |
| 16 | BRCA1 | 10 | 188 | UNRESOLVED | 472 | 257 | NO |
| 17 | PTEN | 10 | 172 | UNRESOLVED | 580 | 83 | NO |
| 18 | MPZ | 10 | 75 | UNRESOLVED | 4 | 3 | NO |
| 19 | TUBA1A | 9 | 83 | UNRESOLVED | 582 | 51 | NO |
| 20 | BRCA2 | 9 | 112 | UNRESOLVED | 119 | 27 | NO |

## Top five exclusion reasons among the 965

Reasons are mutually exclusive and assigned at the first failed attrition stage.

| Rank | Reason | Genes |
|---:|---|---:|
| 1 | FEWER_THAN_TWO_DISEASES_WITH_5_VARIANTS | 296 |
| 2 | FEWER_THAN_TWO_DISEASES_AFTER_PARENT_CHILD_COLLAPSE | 214 |
| 3 | NO_UNAMBIGUOUS_MONDO_DISEASE | 178 |
| 4 | FEWER_THAN_TWO_DISEASES_WITH_10_VARIANTS | 98 |
| 5 | NO_EXACT_RCV_PLP_MISSENSE | 24 |

## ACTB / ACTG1 exact-RCV pilot audit

| Gene | Exact RCVs | Unique P/LP missense variants | Distinct diseases | Variants per disease | Strict | Reason | Credible partners | PIONEER partners |
|---|---:|---:|---:|---|---|---|---:|---:|
| ACTB | 119 | 83 | 4 | MONDO:0009470=4 (4 residues);MONDO:0011823=1 (1 residues);MONDO:0100433=3 (3 residues);MONDO:0700092=2 (2 residues) | NO | FEWER_THAN_TWO_DISEASES_WITH_5_VARIANTS | 509 | 165 |
| ACTG1 | 87 | 51 | 3 | MONDO:0018838=2 (2 residues);MONDO:0019497=1 (1 residues);MONDO:0019719=1 (1 residues) | NO | FEWER_THAN_TWO_DISEASES_WITH_5_VARIANTS | 202 | 47 |

## Data-source failures and limitations

- Frozen RCV archive download resumed from the retained partial; transfer completion, byte size, SHA256, `gzip -t` result, and archive date are recorded in `metadata/phase0R_manifest.tsv` and the download log.
- Targeted E-utilities fallback was not needed: the frozen archive recovered every seed RCV. The fallback is restricted to accessions missing after archive parsing.
- HuRI timed out; IntAct and BioGRID experimental PPI support and precomputed PIONEER/Predictomes partner coverage were used where available. No new complex structures were downloaded from Predictomes.
- Parent-child and unresolved MedGen/MONDO mappings are conservatively excluded from independent disease counts and retained in the audit outputs.

## Phase boundary

Phase 0R only. No structure modeling, interface scoring, MutPred-PPI, FoldX, AlphaMissense comparison, PERMANOVA, classification, tissue matching, or machine learning was run. **Phase 1 was NOT started.**
