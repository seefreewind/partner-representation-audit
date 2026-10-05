from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2];P=R/'results/phase2SR'
new=pd.read_csv(P/'legacy_disease_denominator_reproduction.tsv',sep='\t').set_index(['gene','MONDO'])
old=pd.read_csv(R/'results/phase1B0R3/differential_score_coverage.tsv',sep='\t').set_index(['gene','MONDO'])
assert set(new.index)==set(old.index)
for idx,r in new.iterrows():
    o=old.loc[idx]
    assert r.denominator==o.strict_unique_variants_n,idx
    assert r.ge3==o.variants_with_ge3_scored_Tier1_2_partners_n,idx
    assert r.coverage_pct==o.variants_with_ge3_scored_Tier1_2_partners_pct,idx
    assert r.multilabel_n==o.multilabel_variants_excluded_from_primary_n,idx
print('PASS: every frozen disease-group denominator, >=3 count, rounded coverage and multilabel exclusion reproduced exactly',len(new))
# The legacy code has a score-row-conditioned multilabel exclusion. Retain
# its frozen membership for U1; do not recreate it separately by future model.
m=pd.read_csv(P/'legacy_primary_membership_exact.tsv',sep='\t')
anomaly=m.loc[m.single_label_status.eq('MULTILABEL')&m.included_in_legacy_denominator]
print('Legacy multilabel entries retained because absent from original score-row membership:',len(anomaly),'rows;',len(anomaly[['gene','VariationID']].drop_duplicates()),'variants')
