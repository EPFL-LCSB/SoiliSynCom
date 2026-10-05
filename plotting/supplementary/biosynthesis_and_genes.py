"""Plot NIS and BBB-associated gene-set panels underlying Supplementary V5.

Inputs are read only; no solver calls are made.
"""
from pathlib import Path
import json
import hashlib
import subprocess
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from scipy.io import loadmat
from scipy.sparse import issparse
from common import ROOT, EXTERNAL_DATA_ROOT, OUT, FIG, DATA, ORDER, COLORS, GROUPS, GROUP_COLORS, save, manifest, source_table, sha256

NIS = EXTERNAL_DATA_ROOT / 'Code/Analysis/NIS/output'
figures = []
audit = []


def short(s):
    return re.sub(r'^Thermo|(?:PTYG|Succinate|Xylose)Vitconstrained$', '', str(s))


def frame(path):
    d = pd.read_excel(path, index_col=0)
    d.index = [short(x) for x in d.index]
    d.columns = [short(x) for x in d.columns]
    return d


def taxa_ticks(ax, axis='x', labels=None, size=6):
    labels = labels if labels is not None else ORDER
    if axis == 'x':
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=size)
        ticks = ax.get_xticklabels()
    else:
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=size)
        ticks = ax.get_yticklabels()
    for t, s in zip(ticks, labels):
        if s in COLORS:
            t.set_color(COLORS[s])
    ax.tick_params(length=0, pad=2)


def clean_axis(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)


def group_legend(fig, y=0.018, size=5.7):
    handles = [Patch(facecolor=c, label=g) for g,c in GROUP_COLORS.items()]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5,y),
               ncol=3, frameon=False, fontsize=size, handlelength=1.1, columnspacing=1.3)


def title(ax, letter, text):
    ax.set_title(text, fontsize=8.5, loc='left', pad=6)
    ax.text(-.11,1.04,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)


def nis_ptyg():
    classes = [('Cofactors','Cofactors'),('Nucleotides','Nucleotides'),('Lipids','Lipids'),('Amino_Acids','Amino acids')]
    fig=plt.figure(figsize=(6.5,6.4))
    axes=[fig.add_axes(p) for p in [[.085,.575,.34,.345],[.50,.575,.34,.345],[.2925,.145,.34,.345]]]
    rows=[]
    for i,(fn,label) in enumerate(classes):
        p=ROOT/f'data/graphing/NIS/average_similarity_{fn}PTYG.xlsx'
        d=frame(p).loc[ORDER,ORDER]
        assert d.shape==(21,21) and np.allclose(d,d.T)
        for a in ORDER:
            for b in ORDER:
                rows.append(dict(bbb_class=label,strain_a=a,strain_b=b,recorded_similarity=d.loc[a,b],
                                 interpretation='no_extracted_amino_acid_routes' if fn=='Amino_Acids' else 'class_average_with_missing_routes_zero'))
        if fn=='Amino_Acids':
            assert np.count_nonzero(d.values)==0
            continue
        ax=axes[i]
        z=d.values.copy();np.fill_diagonal(z,np.nan)
        cmap=plt.get_cmap('viridis').copy();cmap.set_bad('white')
        im=ax.imshow(z,vmin=0,vmax=1,cmap=cmap,interpolation='nearest')
        taxa_ticks(ax,'x',ORDER,5.5);taxa_ticks(ax,'y',ORDER,5.5)
        clean_axis(ax);title(ax,'ABC'[i],label)
    ca=fig.add_axes([.90,.34,.018,.3]);cb=fig.colorbar(im,cax=ca);cb.set_label('Class-average reaction-set similarity',fontsize=7);cb.ax.tick_params(labelsize=6)
    group_legend(fig,.015)
    stem=save(fig,'S12_ptyg_biosynthetic_similarity')
    source_table(pd.DataFrame(rows),'S12_ptyg_similarity.tsv')
    figures.append(dict(number=12,title='Biosynthetic-route similarity under PTYG',pages=[stem],
      sources=[f'RawData/NIS/average_similarity_{x[0]}PTYG.xlsx' for x in classes]+['Code/Analysis/NIS/compareBiosynthesis.m'],
      caption='Figure S12. Biosynthetic-route similarity under PTYG. (A–C) Class-average reaction-set similarity for cofactors, nucleotides and lipids across 21 strains. Each value is the mean Jaccard similarity of the reaction unions for the biomass building blocks (BBBs) assigned to that class. The archived comparator assigns zero when either strain lacks an extracted route for a BBB; values therefore combine route availability and reaction-set similarity. White diagonal cells omit self-comparisons. No amino-acid routes were recorded in the PTYG extraction; that all-zero matrix is retained in the source data but is not plotted as biological similarity. Strains follow a common taxonomic order, with label colours identifying taxonomic groups. The results describe modelled routes under the supplied-medium constraints, rather than experimental pathway activities.',
      checks=dict(strains=21,classes=4,amino_acid_entries_all_zero=True,all_matrices_symmetric=True)))


def essentiality():
    files={'PTYG':'data/graphing/GeneEssentiality/EssentialityDataCategorizedAllPTYG.xlsx','M9':'data/graphing/GeneEssentiality/EssentialityDataCategorizedAllSuccinate.xlsx'}
    categories=['Amino acids','Nucleotides','Cofactors','Fatty acids & Lipids','Cell wall components','Polyamines']
    genes={};bycat={};membership=[];summary=[]
    for med,path in files.items():
        sheets=pd.read_excel(ROOT/path,sheet_name=None);genes[med]={};bycat[med]={}
        for s,d in sheets.items():
            cm={str(c).strip():set(str(v).strip() for v in d[c].dropna() if str(v).strip()) for c in d.columns}
            bycat[med][s]=cm;genes[med][s]=set().union(*cm.values())
            summary.append(dict(medium=med,strain=s,unique_genes=len(genes[med][s]),category_memberships=sum(map(len,cm.values()))))
            for cat in categories:
                membership.append(dict(medium=med,strain=s,category=cat,unique_genes_in_category=len(cm[cat]) if cat in cm else np.nan,
                                       category_recorded=cat in cm))
    totals=pd.DataFrame(summary)
    assert [len(genes['PTYG']),len(genes['M9'])]==[21,19]
    assert (totals.query("medium=='PTYG'").unique_genes.min(),totals.query("medium=='PTYG'").unique_genes.median(),totals.query("medium=='PTYG'").unique_genes.max())==(48,58,65)
    assert (totals.query("medium=='M9'").unique_genes.min(),totals.query("medium=='M9'").unique_genes.median(),totals.query("medium=='M9'").unique_genes.max())==(67,85,96)
    shared=[s for s in ORDER if s in genes['M9']]
    overlap=pd.DataFrame([dict(strain=s,shared=len(genes['PTYG'][s]&genes['M9'][s]),PTYG_only=len(genes['PTYG'][s]-genes['M9'][s]),M9_only=len(genes['M9'][s]-genes['PTYG'][s])) for s in shared])
    for r in overlap.itertuples():
        assert r.shared+r.PTYG_only==len(genes['PTYG'][r.strain])
        assert r.shared+r.M9_only==len(genes['M9'][r.strain])
    source_table(totals,'S13_unique_bbb_associated_genes.tsv');source_table(pd.DataFrame(membership),'S13_bbb_category_memberships.tsv');source_table(overlap,'S13_medium_overlap.tsv')
    longgenes=[dict(medium=med,strain=s,gene=g) for med,ss in genes.items() for s,gs in ss.items() for g in sorted(gs)]
    source_table(pd.DataFrame(longgenes),'S13_unique_gene_lists.tsv')
    fig,axs=plt.subplots(2,1,figsize=(6.5,5.1),gridspec_kw={'height_ratios':[1,1.15]})
    fig.subplots_adjust(left=.1,right=.98,bottom=.17,top=.92,hspace=.55)
    ax=axs[0]
    for i,s in enumerate(ORDER):
        p=len(genes['PTYG'][s]);m=len(genes['M9'][s]) if s in genes['M9'] else None
        if m is not None:ax.plot([i,i],[p,m],color='#b0b0b0',lw=.8)
        ax.scatter(i,p,color='#3b6ea8',s=19,zorder=3)
        if m is not None:ax.scatter(i,m,color='#b85b45',s=20,marker='s',zorder=3)
        else:ax.text(i,92,'NA',ha='center',fontsize=6,color='#777777')
    ax.scatter([],[],color='#3b6ea8',s=19,label='PTYG (21 strains)');ax.scatter([],[],color='#b85b45',s=20,marker='s',label='M9 + succinate (19 strains)')
    ax.set_ylim(0,105);ax.set_ylabel('Unique BBB-associated genes');taxa_ticks(ax,'x',ORDER);ax.legend(loc='lower right',frameon=False,fontsize=6.5)
    title(ax,'A','Unique genes across all recorded BBB categories')
    ax=axs[1];bottom=np.zeros(len(shared))
    for col,label,c in [('shared','Shared between media','#717f8c'),('PTYG_only','PTYG only','#3b6ea8'),('M9_only','M9 + succinate only','#b85b45')]:
        ax.bar(range(len(shared)),overlap[col],bottom=bottom,color=c,width=.75,label=label);bottom+=overlap[col].values
    ax.set_ylim(0,105);ax.set_ylabel('Unique genes in the union');taxa_ticks(ax,'x',shared)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,-.38),ncol=3,frameon=False,fontsize=6.5)
    title(ax,'B','Gene identity across media in 19 matched strains')
    for ax in axs:ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    page1=save(fig,'S13_gene_sets_and_medium_overlap')
    fig,axes=plt.subplots(2,1,figsize=(6.5,4.7))
    fig.subplots_adjust(left=.22,right=.88,bottom=.19,top=.93,hspace=.64)
    maximum=max(r['unique_genes_in_category'] for r in membership if pd.notna(r['unique_genes_in_category']))
    for ax,med,letter in zip(axes,['PTYG','M9'],['C','D']):
        vals=np.full((6,21),np.nan)
        for j,s in enumerate(ORDER):
            if s in bycat[med]:
                for i,c in enumerate(categories):
                    if c in bycat[med][s]:vals[i,j]=len(bycat[med][s][c])
        cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#d9d9d9')
        im=ax.imshow(vals,vmin=0,vmax=maximum,cmap=cmap,aspect='auto',interpolation='nearest')
        taxa_ticks(ax,'x',ORDER,6);ax.set_yticks(range(6));ax.set_yticklabels(['Amino acids','Nucleotides','Cofactors','Fatty acids / lipids','Cell wall','Polyamines'],fontsize=6.5)
        ax.tick_params(length=0);clean_axis(ax);title(ax,letter,'PTYG' if med=='PTYG' else 'M9 + succinate')
    cb=fig.colorbar(im,cax=fig.add_axes([.92,.36,.018,.42]));cb.set_label('Unique genes per BBB class',fontsize=7);cb.ax.tick_params(labelsize=6)
    fig.text(.54,.11,'Grey: category not recorded or strain not included',ha='center',fontsize=7)
    group_legend(fig,.006)
    page2=save(fig,'S13_bbb_memberships')
    figures.append(dict(number=13,title='BBB-associated gene sets recorded in two media',pages=[page1,page2],sources=list(files.values()),
      caption='Figure S13. Biomass-building-block-associated gene sets recorded in two media. (A) Unique genes present in the union of all recorded biomass-building-block (BBB) categories for each strain in PTYG and M9 with succinate. Lines connect the same strain; mip and cohm have no M9 entry. (B) Disjoint sets of genes shared between media, recorded only in PTYG, or recorded only in M9 with succinate for the 19 strains represented in both media. (C,D) Unique genes assigned to individual BBB classes in PTYG and M9, respectively. A gene may occur in more than one class, so category counts are not additive totals. Grey cells identify unrecorded categories or unavailable strain–medium combinations, rather than zero genes. The retained PTYG classification has no amino-acid column; absence of that category must not be interpreted as absence of amino-acid-related requirements. Differences between media therefore include differences in recorded category coverage. Counts describe the archived BBB-associated deletion-result sets and do not constitute a separate experimental validation or a direct count of whole-biomass growth-essential genes. Label colours identify taxonomic groups.',
      checks=dict(ptyg_strains=21,m9_strains=19,ptyg_unique_range=[48,65],ptyg_median=58,m9_unique_range=[67,96],m9_median=85,overlap_partition_verified=True,missing_ptyg_amino_acid_category=True)))


def bbb_routes():
    source=[];allrows=[];matrices={}
    for med in ['M9','PTYG']:
        a=frame(NIS/med/'statistics_biosynthesis_a.xlsx');b=frame(NIS/med/'statistics_biosynthesis_b.xlsx')
        assert a.shape==b.shape and len(a)==49
        strains=[s for s in ORDER if s in a.columns];a=a[strains];b=b[strains]
        assert np.all(a.values>=0) and np.all(b.values>=0) and np.array_equal(a.values==0,b.values==0)
        matrices[med]=(a,b)
        for label in a.index:
            for s in strains:
                allrows.append(dict(medium=med,strain=s,bbb=label,enumerated_minimum_size_routes=int(a.loc[label,s]),route_reactions=int(b.loc[label,s]),route_recorded=bool(a.loc[label,s]>0)))
        source.extend([str((NIS/med/f'statistics_biosynthesis_{suffix}.xlsx').relative_to(ROOT)) for suffix in ['a','b']])
    detail=pd.DataFrame(allrows);source_table(detail,'S15_bbb_route_complexity.tsv')
    # Preserve the complete individual-BBB pairwise data rather than plotting 98 matrices.
    pairrows=[];seen={}
    for med in ['M9','PTYG']:
        paths=[p for p in (NIS/med).rglob('similarity_new_*.xlsx') if 'conflicted' not in p.name and 'old' not in p.parts]
        assert len(paths)==49
        for p in paths:
            bid=p.stem.replace('similarity_new_','');assert (med,bid) not in seen;seen[med,bid]=True
            d=frame(p);assert d.shape==(19,19) if med=='M9' else d.shape==(21,21)
            for s in d.index:
                for t in d.columns:pairrows.append(dict(medium=med,bbb_id=bid,strain_a=s,strain_b=t,recorded_similarity=d.loc[s,t],source=str(p.relative_to(ROOT))))
    source_table(pd.DataFrame(pairrows),'S15_complete_bbb_pairwise_similarity.tsv')
    labels=list(matrices['M9'][0].index)
    replacements={'5,10-Methylenetetrahydrofolate':'5,10-Methylene-THF','10-Formyltetrahydrofolate':'10-Formyl-THF','5-Methyltetrahydrofolate':'5-Methyl-THF','S-Adenosyl-L-methionine':'S-Adenosylmethionine','1,2-Diacyl-sn-glycerol':'1,2-Diacylglycerol'}
    display=[replacements.get(x,x) for x in labels]
    pages=[]
    for metric,index,stem,cmap,label in [('size',1,'S15_bbb_route_sizes','viridis','Reactions in smallest extracted route'),('alternatives',0,'S15_bbb_route_alternatives','magma','Enumerated routes of minimum recorded size')]:
        fig=plt.figure(figsize=(6.5,6.7));gs=fig.add_gridspec(1,2,left=.27,right=.875,bottom=.12,top=.925,width_ratios=[19,21],wspace=.055)
        for j,med in enumerate(['M9','PTYG']):
            d=matrices[med][index];v=d.values.astype(float);v[v==0]=np.nan
            ax=fig.add_subplot(gs[0,j]);cm=plt.get_cmap(cmap).copy();cm.set_bad('#dedede')
            kw={'norm':LogNorm(vmin=1,vmax=1000)} if index==0 else {'vmin':0,'vmax':65}
            im=ax.imshow(v,aspect='auto',cmap=cm,interpolation='nearest',**kw)
            taxa_ticks(ax,'x',list(d.columns),5.5);ax.set_yticks(range(49));ax.set_yticklabels(display if j==0 else ['']*49,fontsize=5.8);ax.tick_params(length=0,pad=2);clean_axis(ax)
            ax.set_title('M9 + succinate (19 strains)' if med=='M9' else 'PTYG (21 strains)',fontsize=8,pad=8)
            for y in [19.5,26.5,41.5]:ax.axhline(y,color='white',lw=.8)
        ca=fig.add_axes([.905,.26,.016,.42]);cb=fig.colorbar(im,cax=ca);cb.ax.tick_params(labelsize=6)
        if index==0:cb.set_ticks([1,10,100,1000]);cb.set_ticklabels(['1','10','100','1000'])
        cb.set_label(label,fontsize=7,labelpad=6)
        fig.text(.02,.97,'A' if index==1 else 'B',fontweight='bold',fontsize=11)
        fig.text(.32,.058,'Grey: no extracted route recorded',ha='center',fontsize=7)
        fig.text(.75,.058,'Rows grouped by BBB class',ha='center',fontsize=6.5,color='#555555')
        pages.append(save(fig,stem))
    figures.append(dict(number=15,title='BBB-specific biosynthetic route complexity and multiplicity',pages=pages,sources=source+['Code/Analysis/NIS/compareBiosynthesis.m','Code/Analysis/NIS/output/M9/PathwaySimilarity','Code/Analysis/NIS/output/M9/CostSimilarity','Code/Analysis/NIS/output/PTYG/PathwaySimilarity'],
      caption='Figure S15. Biomass-building-block-specific biosynthetic route complexity and multiplicity. (A) Number of reactions in the smallest extracted route for each of 49 biomass building blocks (BBBs) and each strain under M9 with succinate (19 strains) or PTYG (21 strains). (B) Number of enumerated routes having that minimum recorded size, on a logarithmic colour scale. Horizontal separators distinguish amino acids, nucleotides, cofactors and lipids. Grey cells denote no extracted route recorded, which may reflect direct availability from the medium or the extraction configuration and does not establish biosynthetic incapacity. Alternative counts describe retained enumeration results, not a proven exhaustive set; 11 strain–BBB–medium entries contain 1,000 recorded alternatives. Rows retain archived BBB identifiers where a resolved name was unavailable. THF, tetrahydrofolate. Strains follow the same taxonomic order as the other supplements; label colours identify taxonomic groups. Complete BBB-specific pairwise similarity matrices are provided as source data.',
      checks=dict(bbbs=49,m9_strains=19,ptyg_strains=21,metric_records=len(detail),entries_with_1000_alternatives=int((detail.enumerated_minimum_size_routes==1000).sum()),individual_bbb_matrices=98,pairwise_rows=len(pairrows))))
    audit.append('S15: retained alternative counts are not asserted exhaustive; the exact historical stopping ceiling has not been independently established.')


def getpayload(path, names=None):
    d=loadmat(path,variable_names=names,simplify_cells=True)
    if names is not None:return d
    return next(v for k,v in d.items() if not k.startswith('_'))


def eq(a,b):
    if issparse(a):a=a.toarray()
    if issparse(b):b=b.toarray()
    return np.array_equal(a,b)


def network_signature(model):
    h=hashlib.sha256()
    for k in ['rxns','mets']:
        h.update(('\n'.join(map(str,model[k]))+'\n').encode())
    for k in ['S','lb','ub']:
        v=model[k].toarray() if issparse(model[k]) else model[k]
        a=np.asarray(v,dtype='<f8')
        h.update(str(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()


def core_network():
    base=NIS/'Core'
    paths=[Path(p) for p in subprocess.check_output(['rg','--files',str(base)],text=True).splitlines() if p.endswith('.mat') and 'conflicted' not in p]
    lineage=[];bounddiff=[];sets={};counts=[];matrix_rows=[];maxerr={};configs=[];signatures={}
    for depth in range(3):
        ss={}
        for p in paths:
            if f'_D{depth}_core/Models/' not in str(p):continue
            if not p.stem.startswith(p.parent.name+'_'):continue
            ident=str(p).split(f'_D{depth}')[0].split('/')[-1];s=short(ident)
            model=getpayload(p);ss[s]=set(model['core_model']['rxns'])
            counts.append(dict(depth=depth,strain=s,core_reactions=len(ss[s])))
            signatures[depth,s]=network_signature(model['OriginalGEM'])
            props=model['ReductionProperties']
            assert int(props['D'])==depth
            configs.append(dict(strain=s,depth=depth,L=int(props['L']),add_ETC=props['AddETCAsSubsystem'],aerobic_anaerobic=props['AerobicAnaerobic'],original_GEM_signature=signatures[depth,s],source=str(p.relative_to(ROOT))))
        assert set(ss)==set(ORDER)
        sets[depth]=ss
        recorded=frame(base/f'comparison_D{depth}'/f'D{depth}_rxns_Overlap.xlsx').loc[ORDER,ORDER]
        recalc=np.array([[len(ss[s]&ss[t])/len(ss[s]|ss[t]) for t in ORDER] for s in ORDER])
        maxerr[depth]=float(np.max(np.abs(recorded.values-recalc)))
        assert maxerr[depth]<1e-10, ('Core metric not verified',depth,maxerr[depth])
        print(f'S16: verified D{depth} Jaccard matrix for 21 saved cores',flush=True)
        for i,s in enumerate(ORDER):
            for j,t in enumerate(ORDER):matrix_rows.append(dict(depth=depth,strain_a=s,strain_b=t,jaccard=recalc[i,j]))
    for p in paths:
        if '_D0_core/WorkSpaces/' not in str(p):continue
        ident=str(p).split('_D0')[0].split('/')[-1];s=short(ident)
        orig=getpayload(p,['GEMmodel_orig'])['GEMmodel_orig'];refpath=EXTERNAL_DATA_ROOT/'Code/Analysis/NIS/GEMs/PTYG'/f'{ident}.mat';ref=getpayload(refpath)
        row={'strain':s,'workspace':str(p.relative_to(ROOT)),'reference':str(refpath.relative_to(ROOT))}
        for k in ['rxns','mets','S','lb','ub']:row[k+'_identical']=eq(orig[k],ref[k])
        for k in ['lb','ub']:
            for ix in np.flatnonzero(orig[k]!=ref[k]):
                bounddiff.append(dict(strain=s,reaction=orig['rxns'][ix],bound=k,core_input=float(orig[k][ix]),archived_ptyg_input=float(ref[k][ix])))
        lineage.append(row)
        print(f'S16: checked input lineage for {s}',flush=True)
    assert len(lineage)==21 and all(all(r[k+'_identical'] for k in ['rxns','mets','S']) for r in lineage), lineage
    assert {r['strain'] for r in bounddiff}=={'racd','cohm'}
    assert all(abs(r['core_input'])==100 and abs(r['archived_ptyg_input'])==1000 for r in bounddiff)
    source_table(pd.DataFrame(bounddiff),'S16_core_input_bound_differences.tsv')
    assert all(signatures[d,s]==signatures[0,s] for d in range(3) for s in ORDER)
    source_table(pd.DataFrame(configs),'S16_saved_core_configurations.tsv')
    source_table(pd.DataFrame(lineage),'S16_core_input_lineage.tsv');source_table(pd.DataFrame(matrix_rows),'S16_core_reaction_similarity.tsv');source_table(pd.DataFrame(counts),'S16_core_reaction_counts.tsv')
    subs=frame(base/'comparison_D0/statistics_rxns_initial_subsystems.xlsx')[ORDER]
    assert subs.shape==(6,21)
    subrows=[dict(subsystem=p,strain=s,reactions=int(subs.loc[p,s])) for p in subs.index for s in ORDER]
    source_table(pd.DataFrame(subrows),'S16_initial_subsystem_reactions.tsv')
    fig=plt.figure(figsize=(6.5,6.35));gs=fig.add_gridspec(2,2,left=.1,right=.96,bottom=.15,top=.94,height_ratios=[1,1.8],width_ratios=[1.35,.75],wspace=.45,hspace=.65)
    ax=fig.add_subplot(gs[0,:]);im=ax.imshow(subs.values,aspect='auto',cmap='YlGnBu',vmin=0,vmax=60)
    names=['Glycolysis / gluconeogenesis','Pentose phosphate','TCA cycle','Pyruvate metabolism','Glyoxylate / dicarboxylates','Electron transport chain']
    ax.set_yticks(range(6));ax.set_yticklabels(names,fontsize=6);taxa_ticks(ax,'x',ORDER,6);clean_axis(ax)
    for i in range(6):
        for j in range(21):ax.text(j,i,str(int(subs.iloc[i,j])),ha='center',va='center',fontsize=5,color='white' if subs.iloc[i,j]>38 else '#222222')
    title(ax,'A','Reaction content of the six initial core subsystems')
    # A generous left margin for pathway names; colour bar is unnecessary with cell counts.
    pos=ax.get_position();ax.set_position([.255,pos.y0,.69,pos.height])
    ax=fig.add_subplot(gs[1,0]);z=np.array([[len(sets[0][s]&sets[0][t])/len(sets[0][s]|sets[0][t]) for t in ORDER] for s in ORDER]);np.fill_diagonal(z,np.nan)
    cm=plt.get_cmap('viridis').copy();cm.set_bad('white');im=ax.imshow(z,vmin=0,vmax=1,cmap=cm)
    taxa_ticks(ax,'x',ORDER,5.8);taxa_ticks(ax,'y',ORDER,5.8);clean_axis(ax);title(ax,'B','D0 core reaction similarity')
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.025);cb.ax.tick_params(labelsize=6);cb.set_label('Jaccard similarity',fontsize=7)
    ax=fig.add_subplot(gs[1,1]);vals=[]
    for depth in range(3):
        ss=sets[depth];vals.append([len(ss[s]&ss[t])/len(ss[s]|ss[t]) for i,s in enumerate(ORDER) for t in ORDER[i+1:]])
    vio=ax.violinplot(vals,positions=[0,1,2],showextrema=False,widths=.8)
    for pc,c in zip(vio['bodies'],['#7496af','#4c7598','#214d73']):pc.set_facecolor(c);pc.set_alpha(.7)
    ax.boxplot(vals,positions=[0,1,2],widths=.16,showfliers=False,patch_artist=True,boxprops={'facecolor':'white','edgecolor':'#333333'},medianprops={'color':'#222222','linewidth':1},whiskerprops={'linewidth':.6},capprops={'linewidth':.6})
    ax.set_xticks([0,1,2]);ax.set_xticklabels(['D0','D1','D2']);ax.set_ylim(0,1.02);ax.set_ylabel('Reaction-set Jaccard similarity');ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    title(ax,'C','Core expansion depth')
    ax.text(.5,.03,'210 pairs per depth',transform=ax.transAxes,ha='center',fontsize=6)
    group_legend(fig,.017)
    stem=save(fig,'S16_central_metabolic_cores')
    figures.append(dict(number=16,title='Central-metabolic subsystem composition and core similarity',pages=[stem],
      sources=['Code/Analysis/NIS/output/Core','Code/Analysis/NIS/GEMs/PTYG','Code/Analysis/NIS/compareCore.m'],
      caption='Figure S16. Central-metabolic subsystem composition and core similarity. (A) Reaction counts in six initial subsystems for 21 strains; values are printed in each cell. Subsystem assignments can overlap, so row counts are not additive whole-model totals. (B) Pairwise Jaccard similarity of reaction identifiers in the D0 core networks, recalculated from the saved core reaction sets. White cells omit self-comparisons. (C) Distributions of the same similarity measure over 210 unordered strain pairs at each core-expansion depth (D0–D2). Violin widths show the distribution; boxes indicate the interquartile range, centre lines the median and whiskers the most extreme observations within 1.5 interquartile ranges. D0 denotes the initial subsystem core (L = 0); D1 and D2 include progressively broader connecting paths, with both searches using L = 3. Pairwise values share strains and are shown descriptively without independent-pair significance tests. All 21 D0 input models share the reaction and metabolite identifiers and stoichiometry of the corresponding archived PTYG NIS inputs. Absolute flux-bound limits of 100 were used for one reaction in racd and 21 reactions in cohm where the reference inputs use 1,000; all other compared bounds match. These differences are documented in the source data. Label colours identify taxonomic groups.',
      checks=dict(strains=21,subsystems=6,depths=[0,1,2],pairs_per_depth=210,metric_verified='Jaccard from saved reaction sets',max_matrix_difference=maxerr,all_21_input_structures_identical=True,original_GEM_fields_consistent_across_depths=True,flux_bounds_identical_strains=19,strains_with_bound_cap_differences=['racd','cohm'],bound_fields_differing=len(bounddiff))))
    audit.append('S16: all 21 input network structures match archived PTYG NIS inputs; racd (1 reaction) and cohm (21 reactions) retain absolute bound caps of 100 versus 1,000 in the reference inputs. Figure describes the archived core sets; no reconstruction or solver rerun was performed.')


if __name__=='__main__':
    DATA.mkdir(parents=True,exist_ok=True)
    nis_ptyg();essentiality();bbb_routes();core_network()
    manifest('computational',figures,audit)
    print(json.dumps({'figures':[f['number'] for f in figures],'pages':sum(len(f['pages']) for f in figures),'audit':audit}))
