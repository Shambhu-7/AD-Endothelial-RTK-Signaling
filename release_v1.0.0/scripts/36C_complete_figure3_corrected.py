from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('results/manuscript_final/corrected_final')
SRC=OUT/'source_data'; FIG=OUT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
summary=pd.read_csv(SRC/'Figure3A_pathway_enrichment.csv')
all_pathways=pd.read_csv(SRC/'Figure3C_all_pathway_replication.csv')
assert len(all_pathways)==1168
assert int(all_pathways.replicated_corrected.sum())==1
assert summary.set_index('label').loc['RTK','replicated_corrected']
assert not summary.set_index('label').loc['mTORC1','replicated_corrected']

fig=plt.figure(figsize=(13.5,9))
gs=fig.add_gridspec(2,2,height_ratios=[1.05,1],hspace=.42,wspace=.34)
a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[1,1])
colors=['#277DA1','#F8961E'];x=np.arange(2);width=.32
for offset,ds,nes_col,fdr_col,color in [(-width/2,'Discovery','NES_discovery','padj_discovery',colors[0]),(width/2,'Validation','NES_validation','padj_validation',colors[1])]:
 bars=a.bar(x+offset,summary[nes_col],width,color=color,label=ds)
 for bar,fdr in zip(bars,summary[fdr_col]):
  a.annotate(f'FDR = {fdr:.3g}',(bar.get_x()+bar.get_width()/2,bar.get_height()),xytext=(0,5),textcoords='offset points',ha='center',fontsize=9)
a.set_xticks(x,['RTK signaling','mTORC1 signaling']);a.set_ylim(0,3);a.set_ylabel('Normalized enrichment score');a.legend(frameon=False)
a.set_title('a  Cross-cohort pathway enrichment',loc='left',fontweight='bold')
# Panel b intentionally replaces stale leading-edge overlap; corrected leading-edge
# significance has NOT been recomputed and must not be reported as verified.
y=np.arange(2)
for offset,ds,col,color in [(-.15,'Discovery','padj_discovery',colors[0]),(.15,'Validation','padj_validation',colors[1])]:
 vals=-np.log10(summary[col].clip(lower=np.finfo(float).tiny))
 b.scatter(y+offset,vals,s=90,color=color,label=ds)
b.axhline(-np.log10(.05),color='gray',ls='--',lw=1,label='FDR 0.05')
b.set_xticks(y,['RTK','mTORC1']);b.set_ylabel('−log₁₀(FDR)');b.legend(frameon=False,fontsize=8)
b.set_title('b  Cohort-specific statistical support',loc='left',fontweight='bold')
c.text(.5,.64,'1 / 1,168',transform=c.transAxes,ha='center',va='center',fontsize=32,fontweight='bold')
c.text(.5,.40,'Same-direction pathways significant\nin both datasets (FDR < 0.05)',transform=c.transAxes,ha='center',va='center',fontsize=10)
c.text(.5,.20,'RTK signaling only',transform=c.transAxes,ha='center',fontsize=10)
c.set_xticks([]);c.set_yticks([]);c.set_title('c  Cross-dataset pathway replication',loc='left',fontweight='bold')
for ax in (a,b,c):
 ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
fig.savefig(FIG/'Figure3_Complete.pdf',bbox_inches='tight')
fig.savefig(FIG/'Figure3_Complete.png',dpi=600,bbox_inches='tight')
plt.close(fig)
print('Figure 3 saved (corrected; panel b now displays FDR rather than unverified leading-edge overlap).')
