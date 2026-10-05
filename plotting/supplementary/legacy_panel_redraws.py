"""Re-layout two V2 figures without rerunning metabolic models.

Panel values are checked against manuscript-aligned pooled summary workbooks.
The original profile order and separate panel denominators are retained. The
metabolomics matrix follows ExperimentalValidation.ipynb cell 1 and corrects
V2's two misplaced secretion markers (buol/vpe, Pyridoxine -> Glycerol).
"""
from pathlib import Path
import os, sys, json, hashlib, re
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs'/'supplementary_legacy_redraws'
RAW=ROOT/'data'/'graphing'/'DiMEs'
os.environ['MPLCONFIGDIR']=str(OUT/'mpl_cache')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
sys.path.insert(0,str(ROOT))
from publication_style import apply_publication_style
apply_publication_style(font_family='Arial',extra_rc={'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7,'savefig.bbox':None,'svg.fonttype':'none','pdf.fonttype':42})
for p in ['figures','source_data']:(OUT/p).mkdir(parents=True,exist_ok=True)
GROUPS={**dict.fromkeys(['bbet','cauj','dri','mlot','pzu','tro'],'Alphaproteobacteria'),**dict.fromkeys(['buol','vpe'],'Betaproteobacteria'),**dict.fromkeys(['lgu','lur','pmor','pscf','racd'],'Gammaproteobacteria'),**dict.fromkeys(['cuvk','cxy','mip','rfa'],'Actinobacteria'),**dict.fromkeys(['cpi','flk','mgk'],'Bacteroidetes'),'cohm':'Firmicutes'}
GROUP_COLORS={'Actinobacteria':'#1b9e77','Alphaproteobacteria':'#d95f02','Betaproteobacteria':'#7570b3','Gammaproteobacteria':'#e7298a','Firmicutes':'#66a61e','Bacteroidetes':'#e6ab02'}
DISPLAY={'Inorganics':'Inorganic compounds','AminoAcids':'Amino acids','Glycosides':'Glycosides','Amides':'Amides','FattyAcids':'Fatty acids','HeterocyclicCompounds':'Heterocyclic compounds','SulfurCompounds':'Sulfur compounds','Peptides':'Peptides','Vitamins&Cofactors':'Vitamins & cofactors','Alcohols':'Alcohols','CarboxylicAcids':'Carboxylic acids','Aldehydes':'Aldehydes','Amines':'Amines','Amidines':'Amidines','Sugars':'Sugars','Ketones':'Ketones','Hydrocarbons':'Other carbon compounds','Hydroxides':'Cyanic acid','CarboxylicEsters':'Carboxylic esters','CarbocyclicCompounds':'Carbocyclic compounds','Phenols':'Phenols'}
META_LABEL={'Vitamins&Cofactors':'Vitamins &\ncofactors','Hydrocarbons':'Other carbon\ncompounds','HeterocyclicCompounds':'Heterocyclic\ncompounds','CarboxylicAcids':'Carboxylic acids','AminoAcids':'Amino acids'}
SOURCES=[]
AUDIT={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def record(p,locator):
 SOURCES.append({'path':str(p.relative_to(ROOT)),'sha256':sha(p),'locator':locator})
def save(fig,stem):
 for ax in fig.axes:
  for artist in ax.collections:artist.set_rasterized(False)
 for ext in ['svg','pdf','png']:fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=400,bbox_inches=None,facecolor='white')
 plt.close(fig)
 s=(OUT/'figures'/f'{stem}.svg').read_text()
 assert '<image' not in s
 assert '<text' in s
 return {'stem':stem,'width_inches':float(fig.get_size_inches()[0]),'height_inches':float(fig.get_size_inches()[1]),'svg_live_text':True,'svg_embedded_rasters':0}

# Four exact tables used in the original V2 profile figure. In particular,
# notebook cells 10/11 impose a later shared class order, which V2 did not use.
# The original independent-order algorithm is in cells 0/1/5/6.
profile_specs=[('A','Secretion — basal','products','1st','SecretionProfile1stIterUnionMatrix.xlsx'),('B','Secretion — conditioned','products','2nd','MetabolicProfiles/SecretionProfile2ndIterUnionRaw.xlsx'),('C','Uptake — basal','substrates','1st','UptakeProfile1stIterUnionMatrix.xlsx'),('D','Uptake — conditioned','substrates','2nd','UptakeProfile2ndIterUnionMatrix.xlsx')]
graph=pd.read_excel(RAW/'CompoundClassification.xlsx',sheet_name='Graphing',usecols=['KEGG ID','Name','Classification']).dropna().drop_duplicates('KEGG ID',keep='first')
class_map=dict(zip(graph['KEGG ID'].astype(str).str.strip(),graph['Classification']))
record(RAW/'CompoundClassification.xlsx','Graphing')
record(ROOT/'plotting/main/DiMEs/Uptake_SecretionProfile.ipynb','Cells 0, 1, 5 and 6; cell 10/11 class-order revision not used by V2')
profiles=[]
for letter,title,role,iteration,fn in profile_specs:
 p=RAW/fn
 matrix=pd.read_excel(p,sheet_name='ClassMatrix',index_col=0).drop(columns=['Breadth Score','Normalized Breadth Score'])
 presence=pd.read_excel(p,sheet_name='PresenceMatrix',index_col=0)
 mapping=pd.read_excel(p,sheet_name='MappingUsed')
 spath=RAW/f'summaryAll{iteration}Iteration.xlsx'
 sheets=pd.read_excel(spath,sheet_name=None)
 rawsets={s[:-8]:set(d[role].dropna().astype(str).str.strip()) for s,d in sheets.items() if s.endswith('_summary')}
 ids=sorted(set.union(*rawsets.values())-{'C00001'})
 ids=[c for c in ids if c in class_map]
 assert set(ids)==set(presence.columns),(letter,'compound mismatch')
 checked=0;denom=[]
 for cls in matrix.columns:
  cis=[cid for cid in ids if class_map[cid]==cls]
  denom.append({'panel':letter,'class_original':cls,'class_display':DISPLAY.get(cls,cls),'denominator':len(cis),'compound_ids':';'.join(cis)})
  for sp in matrix.index:
   val=sum(c in rawsets[sp] for c in cis)/len(cis)
   assert np.isclose(matrix.loc[sp,cls],val,rtol=0,atol=1e-12),(letter,sp,cls,val,matrix.loc[sp,cls])
   checked+=1
 for sp in presence.index:
  for cid in presence.columns: assert int(presence.loc[sp,cid])==int(cid in rawsets[sp])
 matrix.to_csv(OUT/'source_data'/f'profile_{letter}_fractions.csv')
 pd.DataFrame(denom).to_csv(OUT/'source_data'/f'profile_{letter}_denominators.csv',index=False)
 profiles.append({'letter':letter,'title':title,'matrix':matrix})
 AUDIT[letter]={'value_cells_checked_against_pooled_summary':checked,'strain_count':matrix.shape[0],'class_count':matrix.shape[1],'compound_count':len(ids),'original_strain_order':list(matrix.index),'original_class_order':list(matrix.columns),'all_match':True}
 record(p,'ClassMatrix, PresenceMatrix, MappingUsed');record(spath,'All *_summary sheets, '+role)

figures=[]
for pi,pair in enumerate([profiles[:2],profiles[2:]],1):
 fig=plt.figure(figsize=(6.5,6.3))
 for k,p in enumerate(pair):
  bottom_y,top_y=(3.40,5.70) if k==0 else (.73,2.91)
  ax=fig.add_axes([.267,bottom_y/6.3,.705,(top_y-bottom_y)/6.3])
  m=p['matrix'];ny,nx=m.shape[1],m.shape[0]
  mesh=ax.pcolormesh(np.arange(nx+1),np.arange(ny+1),m.T.values,cmap='Oranges',vmin=0,vmax=1,linewidth=.25,edgecolors='#ffffff',rasterized=False)
  ax.set_xlim(0,nx);ax.set_ylim(ny,0)
  ax.set_xticks(np.arange(nx)+.5);ax.set_xticklabels(m.index,rotation=90,fontsize=8)
  ax.set_yticks(np.arange(ny)+.5);ax.set_yticklabels([DISPLAY.get(c,c) for c in m.columns],fontsize=7.5)
  ax.tick_params(axis='both',length=0,pad=2)
  for spine in ax.spines.values():spine.set_visible(False)
  for tick,sp in zip(ax.get_xticklabels(),m.index):tick.set_color(GROUP_COLORS[GROUPS[sp]])
  ax.set_title(p['title'],loc='left',fontsize=10,pad=7)
  fig.text(.027,(top_y+.13)/6.3,p['letter'],fontsize=12,fontweight='bold',va='center')
 cax=fig.add_axes([.28,6.07/6.3,.52,.009])
 cb=fig.colorbar(mesh,cax=cax,orientation='horizontal',ticks=[0,.25,.5,.75,1]);cb.ax.tick_params(labelsize=7,length=2,pad=1)
 cb.ax.xaxis.set_label_position('top');cb.set_label('Fraction of metabolites in each group',fontsize=8,labelpad=4)
 handles=[Patch(color=c,label=g) for g,c in GROUP_COLORS.items()]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.005),ncol=3,frameon=False,fontsize=7,handlelength=1,handleheight=.65,columnspacing=1.0)
 figures.append(save(fig,f'S12_v2_profiles_sheet{pi}'))

# Execute only the data-reading/grouping portion of the source notebook, with
# its paths redirected. No plot or output-writing notebook code is executed.
nbpath=ROOT/'plotting/main/DiMEs/ExperimentalValidation.ipynb'
source=''.join(json.loads(nbpath.read_text())['cells'][1]['source'])
source=source[:source.index('# ----------------------------\n# Plot single matrix')]
source=source.replace("matplotlib.use('cairo')","matplotlib.use('Agg')")
cwd=Path.cwd();os.chdir(RAW)
g={};exec(compile(source,str(nbpath)+' cell 1','exec'),g);os.chdir(cwd)
apply_publication_style(font_family='Arial',extra_rc={'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,'savefig.bbox':None,'svg.fonttype':'none','pdf.fonttype':42})
record(nbpath,'Cell 1, pooled summaries; data preparation and matching rules')
record(RAW/'ExperimentalValidation.xlsx','ExperimentalUptakes and ExperiemntalSecretions')
strains=g['species_labels']; rows=[];records=[]
for role,rowids in [('secretion',list(reversed(g['sec_row_ids']))),('uptake',list(reversed(g['upt_row_ids'])))]:
 for mid in rowids:
  row={'role':role,'id':mid,'name':g['name_map'][mid],'class':g['class_map'][mid],'cells':[]}
  for j,sp in enumerate(strains):
   sheet=g['iter2_sheets'][j]
   if role=='uptake':
    in1=mid in g['iter1_subs'].get(sheet,set());in2=mid in g['iter2_subs'][sheet];ring=False
   else:
    idx=g['exp_sec_ids'].index(mid);star=g['star_flags'][idx];precs=g['exp_sec_prec'][idx]
    id1=mid in g['iter1_prods'].get(sheet,set());id2=mid in g['iter2_prods'][sheet]
    in1=id1 or (star and any(p in g['iter1_prods'].get(sheet,set()) for p in precs))
    in2=id2 or (star and any(p in g['iter2_prods'][sheet] for p in precs))
    flags=[]
    if in1:flags.append(id1)
    if in2:flags.append(id2)
    ring=bool(flags) and not any(flags)
   state='both' if in1 and in2 else 'basal_only' if in1 else 'conditioned_only' if in2 else 'absent'
   cell={'state':state,'ring':bool(ring)};row['cells'].append(cell)
   records.append({'role':role,'metabolite_id':mid,'metabolite_name':row['name'],'class_original':row['class'],'strain':sp,'basal':bool(in1),'conditioned':bool(in2),'state':state,'related_compound_only':bool(ring)})
  rows.append(row)
assert len(rows)==90 and len(strains)==21
met=pd.DataFrame(records);met.to_csv(OUT/'source_data'/'metabolomics_prediction_matrix.csv',index=False)
AUDIT['metabolomics']={'rows':90,'uptake_rows':73,'secretion_rows':17,'strains':21,'cells':1890,'uptake_cells':1533,'secretion_cells':357,'predicted_cells':int((met.state!='absent').sum()),'state_counts':met.groupby(['role','state']).size().to_dict().__str__(),'row_order':'Exact V2 top-to-bottom order (reverse of source notebook row arrays); secretion block above uptake.','correction':'Two misplaced blue secretion circles for buol/vpe are on Glycerol in the pooled source, but were on Pyridoxine in the V2 embedded image. Redraw follows source. All 1,533 uptake cells and the other 353 secretion cells had matched in the prior V2 audit.','prior_audit':'Historical V4 bitmap audit (not distributed in this clean repository)'}
for sp in ['buol','vpe']:
 a=met[(met.role=='secretion')&(met.strain==sp)&(met.metabolite_name=='Glycerol')];b=met[(met.role=='secretion')&(met.strain==sp)&(met.metabolite_name=='Pyridoxine')]
 assert a.state.iloc[0]=='both' and b.state.iloc[0]=='absent'

# Break at an original class boundary: upper 45 / lower 45 rows. Text remains
# 7.5 pt at the Word insertion width of 6.5 inches. Every prediction is vector.
COLORS={'both':'#0072B2','basal_only':'#E69F00','conditioned_only':'#CC79A7'}
for pi,rs in enumerate([rows[:45],rows[45:]],1):
 n=len(rs);h=6.3
 fig=plt.figure(figsize=(6.5,h));bottom=.80/h;top=.42/h
 ax=fig.add_axes([.326,bottom,.514,1-bottom-top])
 ax.set_xlim(-.5,20.5);ax.set_ylim(n-.5,-.5)
 ax.set_xticks(range(21));ax.set_xticklabels([sp+'*' if sp=='mip' else sp for sp in strains],rotation=90,fontsize=7.5)
 ax.set_yticks(range(n));ax.set_yticklabels([r['name'] for r in rs],fontsize=7.5)
 ax.tick_params(axis='both',length=0,pad=2)
 for sp in ax.spines.values():sp.set_visible(False)
 ax.grid(axis='x',color='#cccccc',linestyle=':',linewidth=.4)
 bounds=[];start=0
 for i in range(1,n+1):
  if i==n or (rs[i]['class'],rs[i]['role'])!=(rs[start]['class'],rs[start]['role']):
   bounds.append((start,i-1,rs[start]['class'],rs[start]['role']));start=i
 for k,(st,en,cls,role) in enumerate(bounds):
  if k%2==0:ax.axhspan(st-.5,en+.5,color='#f2f2f2',zorder=-2,lw=0)
  ax.plot([20.7,20.85,20.85,20.7],[st-.4,st-.4,en+.4,en+.4],color='#888888',lw=.6,clip_on=False)
  label=META_LABEL.get(cls,DISPLAY.get(cls,cls))
  if pi==2 and st==0 and cls=='Glycosides':label='Glycosides\n(continued)'
  ax.text(21.02,(st+en)/2,label,fontsize=7,va='center',ha='left',clip_on=False)
 for i,r in enumerate(rs):
  if i>0 and r['role']!=rs[i-1]['role']:ax.axhline(i-.5,color='#777777',lw=.8)
  for j,c in enumerate(r['cells']):
   if c['state']=='absent':continue
   col=COLORS[c['state']]
   ax.scatter(j,i,marker='s' if r['role']=='uptake' else 'o',s=13,facecolors='none' if c['ring'] else col,edgecolors=col if c['ring'] else '#333333',linewidths=.4,zorder=3)
 fig.text(.015,1-.12/h,'A' if pi==1 else 'B',fontweight='bold',fontsize=12,va='top')
 fig.text(.065,1-.12/h,'Predictions for experimentally listed metabolites'+(' — continued' if pi==2 else ''),fontsize=9,va='top')
 handles=[Patch(color=c,label=l) for c,l in zip(COLORS.values(),['Both environments','Basal only','Conditioned only'])]
 handles += [Line2D([],[],marker='s',ms=3.5,mfc='white',mec='#333333',ls='',label='Uptake'),Line2D([],[],marker='o',ms=4,mfc='white',mec='#333333',ls='',label='Secretion'),Line2D([],[],marker='o',ms=4,mfc='none',mec=COLORS['both'],ls='',label='Related-compound match')]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.15/h),ncol=3,frameon=False,fontsize=7,handlelength=1,handletextpad=.5,columnspacing=1.2)
 fig.text(.5,.025/h,'* Basal mip results unavailable',ha='center',va='bottom',fontsize=7)
 figures.append(save(fig,f'S15_v2_metabolomics_sheet{pi}'))

manifest={'figures':figures,'audit':AUDIT,'sources':SOURCES,'profile_caption_addendum':'The four original panel-specific denominators and strain/group orders are retained. Panels A–B are shown on the first sheet and C–D on the second; strain-label colours identify taxonomic groups. The source panels match notebook cells 0/1/5/6, before the shared class-order revision in cells 10/11.','metabolomics_caption_addendum':'The matrix continues from sheet A to sheet B, with the original top-to-bottom row order retained; the glycoside group spans both sheets. Basal predictions are unavailable for mip; its conditioned markers do not demonstrate a gain relative to basal. Related-compound-only matches are hollow circles. Two secretion markers (buol and vpe) were corrected from the V2 Pyridoxine row to Glycerol, in agreement with the pooled source summaries.','missing_vs_zero':'Basal mip is absent from the profile panels A/C. In the comparison matrix, missing basal mip records follow the original notebook set-membership convention; a lack of a marker is not an experimental absence.'}
(OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'figures':figures,'profile_checks':{k:v['value_cells_checked_against_pooled_summary'] for k,v in AUDIT.items() if k in 'ABCD'},'metabolomics':AUDIT['metabolomics']},indent=2))
