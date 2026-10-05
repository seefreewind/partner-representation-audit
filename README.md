# Partner representation audit for allelic heterogeneity

Frozen code and derived evidence for a staged audit of technical scoreability, phenotype-balanced coverage and incremental phenotype information in partner-aware protein representations.

Versioned archive: [Zenodo v1.0.0](https://doi.org/10.5281/zenodo.23161598). All-version DOI: [10.5281/zenodo.23161597](https://doi.org/10.5281/zenodo.23161597). The release tag preserves the original frozen archive; this README adds its assigned DOI.

## Frozen conclusions

- Original disruption arm: `NO_GO_EDGETIC_REPRESENTATION`
- Sequence representation arm: `PARTIAL_SEQUENCE_RESCUE`
- Preregistered phenotype test: `NO_GLOBAL_PARTNER_IDENTITY_SUPPORT`

The phenotype specification was cryptographically frozen before disease-label association statistics: SHA256 `d657fab6b24f9aa757f9773bb1e9e0bd5d850e743eab6bdf819f126326ced9c3`.

## Contents and reuse

`results/phase2P0/analysis_specification.md` describes the frozen test; `results/phase2P0/matrices/` contains the observation units and representations. `results/phase2P1/` contains exact estimates, saved permutation/bootstrap outputs and sensitivity outputs. `results/phase2S_full/` contains reconciled coverage, pair statuses and internal structural-validation outputs. Historical results retain separate denominators under their original phase directories. `submission_BIB/PartnerRepresentation_BIB_Supplementary_Tables.xlsx` retains full-precision tables. `submission_BIB/figure_source_map.tsv` maps every figure panel to source data; `submission_BIB/figures/source/` provides presentation scripts.

All numerical files are archival copies of completed outputs. ARCHIVE_MANIFEST.tsv records source and archive SHA256 hashes. Tables above 25 MB are losslessly gzipped; decompress to their original path before using historical scripts. macOS resource forks, credentials, personal logs, manuscripts, correspondence and raw licensed database copies are excluded.

## Reproducibility routes

1. Inspect saved exact results and the specification without running any inference or statistics.
2. To recreate presentation figures, decompress `.tsv.gz` tables, install the plotting dependencies and run the selected scripts in `submission_BIB/figures/source/`. They read the completed outputs, not model inference.
3. Historical analysis scripts are retained for methodological transparency. They expect the documented project layout and some original local model/environment paths. Full end-to-end reruns require acquisition of third-party resources from their official sources and adjustment of local paths. This archive is not a containerized one-command model rerun.

Environment snapshots are `results/phase2SR/environment_sppider.txt` and `environment_xpair.txt`; model revisions and hashes are under `external_models/`. Model weights, large encoder embeddings and full WT residue tracks are not redistributed. Their official version/revision identifiers are recorded in the manifests and source reports. ClinVar and molecular-resource records are referenced by identifiers; restricted source databases (including OMIM), third-party training datasets, model packages and copyrighted articles must be obtained independently under their own terms.

The authors' original code and derived tables/figures are released under CC BY 4.0. This licence does not override rights or licences for any referenced third-party software, models or source databases. Model authors retain ownership of their software; upstream packages and weights are not included.

## Authors

Sisi Xu; Yu Zhang; Tiantian Zeng; Jiechen Liu; Ling Qiu. Corresponding author: Ling Qiu. No journal name is included in the repository name. Archived files retain their historical paths for provenance.

## Reassemble the large pairwise table

GitHub stores the compressed pairwise table in three lossless parts below its per-file limit. Reassemble before decompressing:

```sh
cat results/phase2S_full/partner_specificity_pairwise.tsv.gz.part* > results/phase2S_full/partner_specificity_pairwise.tsv.gz
gzip -dk results/phase2S_full/partner_specificity_pairwise.tsv.gz
```

The SHA256 of the reassembled gzip stream is b56d407aa9a796a54fee822e514a7f750d9902765c5a1af4c5574364aa3ad0f7. Each part and the uncompressed original retain hashes in ARCHIVE_MANIFEST.tsv.
