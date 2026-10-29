from pathlib import Path
import numpy as np
import pandas as pd

R=Path('results'); check=R/'manuscript_final/reproducibility/dispersion_check'
out=R/'manuscript_final/corrected_final/source_data';out.mkdir(parents=True,exist_ok=True)
focus={'RTK':'REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES','mTORC1':'HALLMARK_MTORC1_SIGNALING'}
paths={
 'discovery_primary':check/'pathway_robustness/primary_all_pathways.csv',
 'discovery_lowcell':check/'pathway_robustness/exclude_61cell_all_pathways.csv',
 'discovery_composition':check/'pathway_robustness/composition_adjusted_all_pathways.csv',
 'validation_primary':R/'GSE327289_robustness/20_primary_all_pathways.csv',
 'validation_lowcell':R/'GSE327289_robustness/20_exclude_under100_all_pathways.csv',
 'validation_composition':R/'GSE327289_robustness/20_composition_adjusted_all_pathways.csv'
}
rows=[]
for model,path in paths.items():
 df=pd.read_csv(path)
 assert df.pathway.is_unique, f'Duplicate pathway: {path}'
 for label,gene_set in focus.items():
  hit=df.loc[df.pathway.eq(gene_set)]
  assert len(hit)==1,(path,gene_set)
  r=hit.iloc[0]
  rows.append(dict(label=label,model=model,NES=float(r.NES),FDR=float(r.padj),available=True))
res=pd.DataFrame(rows);res['significant']=res.FDR.lt(.05)
assert res.groupby('label').significant.sum().to_dict()=={'RTK':6,'mTORC1':4}
res.to_csv(out/'Figure4_robustness_input.csv',index=False)

cross=pd.read_csv(check/'corrected_pathway_comparison.csv')
assert len(cross)==1168 and cross.pathway.is_unique
assert {'pathway','NES_discovery','NES_validation','padj_discovery','padj_validation'}.issubset(cross.columns)
cross['replicated_corrected']=(cross.padj_discovery.lt(.05)&cross.padj_validation.lt(.05)&np.sign(cross.NES_discovery).eq(np.sign(cross.NES_validation)))
assert cross.replicated_corrected.sum()==1
assert cross.loc[cross.replicated_corrected,'pathway'].tolist()==[focus['RTK']]
cross.to_csv(out/'Figure3C_all_pathway_replication.csv',index=False)
summary=cross.loc[cross.pathway.isin(focus.values()),['pathway','NES_discovery','NES_validation','padj_discovery','padj_validation','replicated_corrected']].copy()
summary['label']=summary.pathway.map({v:k for k,v in focus.items()})
summary=summary.set_index('label').loc[['RTK','mTORC1']].reset_index()
summary.to_csv(out/'Figure3A_pathway_enrichment.csv',index=False)
print('PASS: 1/1168 replicated, RTK 6/6, mTORC1 4/6')
print('Saved:',out)
