# Phase 2S-R offline encoder and denominator recovery

## Gate

`GO_RESUME_PHASE2S`

Both exact encoders loaded offline and both official example workflows passed. Canonical sequence recovery, exact legacy denominator reproduction and the disease-blind technical pilot passed. Recovery ends here. The full Phase 2S benchmark has not started.

## Preserved state

All 124 frozen files were rehashed and remain unchanged, including original results and neural heads. U1=42, U2=132 and 15,942 frozen gene–partner pairs remain fixed.

## Exact legacy reproduction

An independent recomputation used original scores only and traced `scripts/phase1b0r3_finalize.py::make_membership`. It reproduced every gene primary denominator, ≥3 partner count, rounded coverage, balance class and final gate. A separate audit reproduced all 109 disease-group denominators, ≥3 counts, rounded percentages and multilabel exclusions. Results: **9 genes ≥70%; 3 balanced: FGFR2, MEN1, TUBA1A. LEGACY_DENOMINATOR_EXACT_PASS=YES.**

The legacy rule first constructs disease membership from Phase 1A `primary_cohort_eligible=YES`. Its gene-level multilabel exclusion is then restricted to variants appearing in original score rows. Consequently 674 multilabel variants absent from score rows remain in legacy primary denominators (1,555 disease-membership rows). This exactly explains why an ordinary full-cohort single-label registry disagreed. We preserve this behavior explicitly rather than retrospectively fix the original result. Group balance uses original strict counts and excludes multilabel variants for the retained-group ≥10 gate, as in the original script.

COMMON_PRIMARY_MEMBERSHIP contains **9,978 included unique gene–variant units**. U1 retains exact immutable legacy membership. U2 outside U1 uses frozen mapping eligibility and one strict label, fixed before inference. Both models will evaluate the same frozen identities; no future method-specific membership or structure filter is permitted. This U2 extension is explicit and does not claim to have historical legacy scores. Pair-space denominators remain separate from variant ≥3 representation.

## Canonical sequences

7,300 requested; 7,300 recovered. {"EXACT": 7300}. All 267 prior sequences retain identical hashes.

| universe   |   distinct_accessions |   recovered_accessions |   accession_pct |   frozen_PPI_pairs |   complete_pairs |   pair_completeness_pct |
|:-----------|----------------------:|-----------------------:|----------------:|-------------------:|-----------------:|------------------------:|
| U1         |                  5844 |                   5844 |             100 |              11354 |            11180 |                 98.4675 |
| U2         |                  7300 |                   7300 |             100 |              15942 |            15471 |                 97.0455 |

UniProtKB release 2026_03 (2026-09-02) and authoritative batch responses are preserved. Partners without frozen accessions remain unavailable; no isoform or accession identity was substituted.

## Exact models and runtime

| repo                         | revision                                 | weight_file       |      bytes | sha256                                                           |
|:-----------------------------|:-----------------------------------------|:------------------|-----------:|:-----------------------------------------------------------------|
| facebook/esm2_t33_650M_UR50D | 08e4846e537177426273712802403f7ba8261b6c | model.safetensors | 2609506392 | a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0 |
| ElnaggarLab/ankh-large       | 74b371dbfa3ee0a05d32ae74df0c2e0b82d6b9a6 | pytorch_model.bin | 7515077421 | 517b6e8b279dedcb477af240b35c46bd6eb3307723eb281e60d4b2c8a87b889b |

The user authorized local download after switching networks on 2026-10-02. Official revision, byte sizes, SHA256 and tokenizer/config hashes are recorded in the manifests.

Independent Python 3.12.13 environments reside at `/Users/zy/.codex/tools/phase2s_sppider` and `/Users/zy/.codex/tools/phase2s_xpair`. Both imports and package consistency passed. Package pins are recorded in environment_sppider.txt and environment_xpair.txt. exFAT package copying failed; a workspace sparse-image mount was denied, so reusable user-level SSD environments were installed under the existing environment authorization. Task outputs remain in the workspace.

SPPIDER v2.0-manuscript uses unchanged ESM2, heads, MAX_TOKENS=1024 and STRIDE=512. No installed SPPIDER source patch was required. X-PAIR 0.1.0 uses the bundled checkpoint, Ankh T5EncoderModel/AutoTokenizer with dimension 1536 and the frozen 50–2000 length rule. Only its encoder path and local_files_only were patched; the exact diff is preserved. No unsupported protein was truncated.

## Example and pilot validation

SPPIDER used its official TP53 three-partner FASTAs and validated all six receptor/peptide tracks. X-PAIR ran the documented official prediction workflow on those same three pairs; no executable example FASTA was bundled in its documentation. Both output lengths matched WT sequences, all probabilities were finite and within [0,1].

The 20 genes were fixed before inference: ten deterministic U1 genes and ten additional U2 genes. The first three lexicographic frozen partners per gene were selected for technical checks only. ABCA4 has no frozen partners and contributes no pilot pairs; it remains selected and unavailable. Total pilot pairs=57. No full-universe coverage or rescue conclusion is inferred from this subset.

| method   |   selected_pairs |   scored_pairs | statuses                                       |
|:---------|-----------------:|---------------:|:-----------------------------------------------|
| sppider  |               57 |             57 | {"SCORED": 57}                                 |
| xpair    |               57 |             42 | {"SCORED": 42, "MODEL_LENGTH_UNSUPPORTED": 15} |

SPPIDER retains complete float32 receptor and peptide residue arrays in keyed npz caches. Its original window merge is used; short and long sequences, both endpoints and internal coverage were checked without padding or truncation. X-PAIR raw six-decimal tracks are archived before the unchanged official CLI rounds its final output to three decimals and removes its transient file. This wrapper preserves output only; inference and source remain unchanged. A second pass reuses the same embeddings and must reproduce every exported probability exactly. No full-precision tensor equivalence beyond the official six-decimal export is claimed.

Variant positions were extracted from cached WT tracks only when the frozen canonical mapping and reference residue matched. No mutant embeddings were generated. The ≥3 calculation counts distinct frozen partners with available scores and remains separate from pair-space availability. No coordinate-offset or sequence-mismatch failure remained. Evidence: sppider_window_validation.tsv, pilot_variant_scores_*.tsv, pilot_variant_partner_count_validation.tsv and inference_cache_manifest.tsv.

## Stop

Original NO_GO_EDGETIC_REPRESENTATION was not altered.

No sequence-model rescue verdict was issued.

Disease-label phenotype association testing was NOT started.

Wait for authorization before the full Phase 2S benchmark.
