"""Build V4 S8-S11 from archived, manuscript-aligned DiMEs summaries.

Run with /opt/anaconda3/envs/cobra/bin/python; no solver is used.
The source workbooks and MATLAB archives are read only. All calculations
and eligible producer-metabolite-consumer triples are exported as TSV.
"""
from collections import Counter, defaultdict
from pathlib import Path
import re
import textwrap
import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from common import ROOT, EXTERNAL_DATA_ROOT, OUT, FIG, DATA, ORDER, COLORS, style, save, manifest, source_table, sha256

BASE = ROOT / 'data' / 'graphing' / 'DiMEs'
CODE = EXTERNAL_DATA_ROOT / 'Code' / 'Analysis' / 'DiMEs'
WATER = 'C00001'
MAIN_CLASSES = {'Sugars', 'CarboxylicAcids', 'Vitamins&Cofactors'}
RED, BLUE, GRAY = '#9E4A4F', '#4B7FAD', '#999999'
CLASS_LABELS = {'AminoAcids':'Amino acids', 'CarboxylicAcids':'Carboxylic acids',
    'Vitamins&Cofactors':'Vitamins & cofactors', 'SulfurCompounds':'Sulfur compounds',
    'HeterocyclicCompounds':'Heterocyclic compounds', 'FattyAcids':'Fatty acids',
    'CarboxylicEsters':'Carboxylic esters', 'CarbocyclicCompounds':'Carbocyclic compounds'}
DISPLAY_GROUP_LABELS = {'Hydrocarbons':'Other carbon compounds', 'Hydroxides':'Cyanic acid'}


def label_class(x):
    return DISPLAY_GROUP_LABELS.get(x, CLASS_LABELS.get(x, x))


def source_group_order(x):
    """Preserve archived group ordering and panel letters after display corrections."""
    return CLASS_LABELS.get(x, x)


def load_summary(path):
    sheets = pd.read_excel(path, sheet_name=None)
    result = {r:{} for r in ('products','substrates')}
    for sheet, df in sheets.items():
        if not sheet.endswith('_summary'):
            continue
        strain = sheet[:-8]
        assert strain in ORDER, strain
        for role in result:
            result[role][strain] = set(df[role].dropna().astype(str).str.strip())
    return result


def met_map(profiles):
    out = defaultdict(set)
    for strain, compounds in profiles.items():
        for met in compounds:
            out[met].add(strain)
    return dict(out)


def union(profiles):
    return set().union(*profiles.values())


def color_ticks(ax, axis='x'):
    ticks = ax.get_xticklabels() if axis == 'x' else ax.get_yticklabels()
    for t in ticks:
        if t.get_text() in COLORS:
            t.set_color(COLORS[t.get_text()])


def panel_label(fig, letter, x, y):
    fig.text(x, y, letter, fontsize=12, fontweight='bold', va='top')


def clean_axes(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)


def normalize_ids(values):
    return {re.search(r'C\d{5}', str(x)).group() for x in np.asarray(values).ravel()}


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    style()
    source_paths = {('Basal','pooled'):BASE/'summaryAll1stIteration.xlsx',
       ('Conditioned','pooled'):BASE/'summaryAll2ndIteration.xlsx',
       ('Basal','optimal'):BASE/'summaryAll100_1stIteration.xlsx',
       ('Conditioned','optimal'):BASE/'summaryAll100_2ndIteration.xlsx'}
    summaries = {key:load_summary(path) for key,path in source_paths.items()}
    basal, cond = summaries[('Basal','pooled')], summaries[('Conditioned','pooled')]
    class_file = BASE/'CompoundClassification.xlsx'
    mapping = pd.read_excel(class_file, sheet_name='Graphing').dropna(subset=['KEGG ID','Classification'])
    mapping = mapping.drop_duplicates('KEGG ID', keep='first')
    classes = dict(zip(mapping['KEGG ID'], mapping['Classification']))
    names = dict(zip(mapping['KEGG ID'], mapping['Name']))
    all_ids = set().union(*(union(p) for data in summaries.values() for p in data.values()))
    assert all_ids <= set(classes), sorted(all_ids-set(classes))
    source_table(mapping.rename(columns={'KEGG ID':'metabolite','Name':'name','Classification':'class'}),
                 'dimes_compound_mapping.tsv')
    presence_rows = []
    for (condition, selection), data in summaries.items():
        for role, profiles in data.items():
            for strain, mets in profiles.items():
                for met in sorted(mets):
                    presence_rows.append({'condition':condition,'selection':selection,'role':role,
                        'strain':strain,'metabolite':met,'name':names[met],'class':classes[met]})
    source_table(pd.DataFrame(presence_rows), 'dimes_all_profile_memberships.tsv')
    audit = []
    figures = []
    display_crosswalk = [{'source_group':old, 'display_group':new,
        'n_catalog_metabolites':sum(c==old for c in classes.values()),
        'metabolites':';'.join(sorted(m for m,c in classes.items() if c==old)),
        'membership_changed':False}
        for old,new in DISPLAY_GROUP_LABELS.items()]
    assert [r['n_catalog_metabolites'] for r in display_crosswalk] == [9,1]
    source_table(pd.DataFrame(display_crosswalk), 'S08_S10_group_display_labels.tsv')
    audit.append({'group_display_label_corrections':display_crosswalk,
        'reason':'Archived Hydrocarbons includes oxygenated and sulfur-containing compounds; '
                 'archived Hydroxides contains only cyanic acid. Only display labels changed; '
                 'source group memberships, denominators, order and exchange links are preserved.'})

    # S8: use a common classified-catalog denominator across all four panels.
    # This deliberately differs from the notebook's role/environment-specific
    # presence-column denominator, so colors are directly comparable.
    coverage_rows, denominators = [], []
    coverages = {}
    class_catalog = {c:{m for m in classes if classes[m]==c and m!=WATER}
                     for c in sorted(set(classes.values()),key=source_group_order)}
    assert set().union(*class_catalog.values()) == (
        union(basal['substrates'])|union(basal['products'])|
        union(cond['substrates'])|union(cond['products']))-{WATER}
    for condition, data in [('Basal',basal),('Conditioned',cond)]:
        for role, profiles in data.items():
            for c, mids in class_catalog.items():
                denominators.append({'condition':condition,'role':role,'class':c,'denominator':len(mids),
                    'metabolites':';'.join(sorted(mids))})
                for s in ORDER:
                    value = len(profiles[s]&mids)/len(mids) if s in profiles else np.nan
                    coverage_rows.append({'condition':condition,'role':role,'strain':s,'class':c,
                        'n_metabolites':len(profiles[s]&mids) if s in profiles else np.nan,
                        'denominator':len(mids),'fraction':value})
                    coverages[(condition,role,s,c)] = value
    source_table(pd.DataFrame(coverage_rows), 'S08_class_coverage.tsv')
    source_table(pd.DataFrame(denominators), 'S08_class_denominators.tsv')

    pool_files = [CODE/'EXCMetabolomics.mat',CODE/'MetabolomicsPool.mat',CODE/'EnviromentPool.mat']
    mapped, pool_b, pool_c = [normalize_ids(loadmat(p, simplify_cells=True)['EXCMetabo']) for p in pool_files]
    products = union(basal['products'])
    assert (len(mapped),len(pool_b),len(pool_c),len(products))==(120,145,315,255)
    assert mapped <= pool_b and pool_c == pool_b|products
    pool_rows=[]
    for met in sorted(pool_c):
        pool_rows.append({'metabolite':met,'name':names.get(met,met),'mapped_source':met in mapped,
          'basal_input':met in pool_b,'basal_product':met in products,'conditioned_input':True,
          'added_to_basal':met in pool_c-pool_b})
    source_table(pd.DataFrame(pool_rows), 'S08_resource_pool_membership.tsv')
    source_table(pd.DataFrame([{'step':'Mapped source list','n_metabolites':120},
        {'step':'Basal input','n_metabolites':145},{'step':'Basal candidate products','n_metabolites':255},
        {'step':'Conditioned input','n_metabolites':315},{'step':'Additions to mapped list','n_metabolites':25},
        {'step':'Additions to basal input','n_metabolites':170}]), 'S08_resource_pool_counts.tsv')
    pages=[]
    for page,role,color in [(1,'substrates',RED),(2,'products',BLUE)]:
        fig=plt.figure(figsize=(6.5,6.5))
        row_classes=list(class_catalog)
        cmap=LinearSegmentedColormap.from_list('coverage',['#FFFFFF',color])
        cmap.set_bad('#D8D8D8')
        for k,(condition,data) in enumerate([('Basal',basal),('Conditioned',cond)]):
            y0=.605 if k==0 else .22
            ax=fig.add_axes([.29,y0,.605,.30])
            arr=np.array([[coverages.get((condition,role,s,c),np.nan) for s in ORDER] for c in row_classes])
            im=ax.imshow(np.ma.masked_invalid(arr),aspect='auto',cmap=cmap,vmin=0,vmax=1,interpolation='nearest')
            ax.set_yticks(range(len(row_classes)))
            ns={(r['condition'],r['role'],r['class']):r['denominator'] for r in denominators}
            ax.set_yticklabels([f'{label_class(c)} ({ns.get((condition,role,c),0)})' for c in row_classes],fontsize=7)
            ax.set_xticks(range(len(ORDER)));ax.set_xticklabels(ORDER,rotation=90,fontsize=7)
            color_ticks(ax);clean_axes(ax)
            ax.set_title(f'{condition}: {"uptake" if role=="substrates" else "secretion"} group coverage',loc='left',pad=8)
            panel_label(fig,chr(65+(page-1)*2+k),.025,y0+.332)
        cax=fig.add_axes([.92,.23,.022,.675]);cb=fig.colorbar(im,cax=cax,ticks=[0,.25,.5,.75,1])
        cb.ax.tick_params(labelsize=7);cb.ax.set_title('Fraction',fontsize=7,pad=8)
        fig.text(.29,.135,'Parentheses: common group denominator. Gray: unavailable basal mip.',fontsize=7)
        fig.text(.5,.065,'Fixed strain order; mip has no basal profile. Water is excluded from group fractions.',
                 ha='center',fontsize=7)
        pages.append(save(fig,f'S08_class_coverage_{page:02d}'))
    figures.append({'number':8,'title':'Predicted uptake and secretion across metabolite classes.',
       'caption':'(A,B) Feasible uptake profiles in the basal and conditioned environments. (C,D) Corresponding secretion profiles. '
       'Each cell is the fraction of metabolites in the indicated class that occurs in at least one retained DiMEs profile for that strain. '
       'The common denominator, shown in parentheses, is the number of compounds in that class in the classified catalog, identical across '
       'A-D. This catalog is the union of all four profiles (267 compounds after excluding water, C00001). White indicates zero class coverage; '
       'the unavailable basal mip profile is gray. Strain order and the 0-1 color scale are fixed; strain-label colors denote taxonomic groups. '
       'Profiles combine alternatives from 10-100% of optimal biomass yield and describe predicted capabilities, not measured exchange.',
       'pages':pages,'sources':[str(p.relative_to(ROOT)) for p in [source_paths[('Basal','pooled')],source_paths[('Conditioned','pooled')],class_file,*pool_files]],
       'page_descriptions':{'S08_class_coverage_01':'A, basal uptake; B, conditioned uptake. Common class denominators exclude water.',
          'S08_class_coverage_02':'C, basal secretion; D, conditioned secretion. Same class denominators and strain order; no panel E.'},
       'checks':{'basal_strains':20,'conditioned_strains':21,'pool_sizes':[120,145,255,315],'pool_union_exact':True,
        'water_excluded_from_class_fractions':True,'common_catalog_denominator':267,
        'source_tables':['S08_class_coverage.tsv','S08_class_denominators.tsv','S08_resource_pool_membership.tsv']}})

    # S9: distinct questions on separate sheets. Exact sharing counts replace
    # mixed percentages, and supplier counts retain the recipient perspective.
    bins=['1 strain','2-5 strains','6-10 strains','>10 strains']
    bin_colors=['#D8D8D8','#9DB8CB','#5D8DA8','#245875']
    sharing_detail=[];sharing_counts=[]
    specs=[('Basal','substrates'),('Conditioned','substrates'),('Basal','products'),('Conditioned','products')]
    for condition,role in specs:
        mm=met_map(summaries[(condition,'pooled')][role]);counts=Counter()
        for met,ss in sorted(mm.items()):
            if met==WATER:continue
            n=len(ss);idx=0 if n==1 else 1 if n<=5 else 2 if n<=10 else 3
            counts[bins[idx]]+=1
            sharing_detail.append({'condition':condition,'role':role,'metabolite':met,'name':names[met],
                'n_strains':n,'bin':bins[idx],'strains':';'.join(s for s in ORDER if s in ss)})
        for label in bins:sharing_counts.append({'condition':condition,'role':role,'bin':label,'n_metabolites':counts[label],
            'total_metabolites':sum(counts.values()),'percent':100*counts[label]/sum(counts.values())})
    source_table(pd.DataFrame(sharing_detail), 'S09_metabolite_sharing.tsv')
    source_table(pd.DataFrame(sharing_counts), 'S09_sharing_bins.tsv')
    basal_producers=met_map(basal['products']);donor_rows=[];donor_counts=[]
    for recipient in ORDER:
        cnt=Counter()
        for met in sorted(cond['substrates'][recipient]-pool_b):
            donors=basal_producers.get(met,set())-{recipient}
            category='0 other' if not donors else '1 other' if len(donors)==1 else '2+ others'
            cnt[category]+=1
            donor_rows.append({'recipient':recipient,'metabolite':met,'name':names[met],
               'n_other_donors':len(donors),'bin':category,'other_donors':';'.join(s for s in ORDER if s in donors),
               'recipient_is_basal_producer':recipient in basal_producers.get(met,set())})
        donor_counts.append({'strain':recipient,**{c:cnt[c] for c in ['0 other','1 other','2+ others']},'total':sum(cnt.values())})
    ddf=pd.DataFrame(donor_counts);details=pd.DataFrame(donor_rows)
    assert ddf[['0 other','1 other','2+ others']].sum().tolist()==[5,65,1826]
    source_table(details,'S09_other_donor_memberships.tsv');source_table(ddf,'S09_other_donor_counts.tsv')
    exact_counts=[]
    for condition,role in specs:
        vals=Counter(r['n_strains'] for r in sharing_detail if r['condition']==condition and r['role']==role)
        for n in range(1,22):exact_counts.append({'condition':condition,'role':role,'n_strains':n,'distinct_metabolites':vals[n]})
    exact_counts=pd.DataFrame(exact_counts);source_table(exact_counts,'S09_exact_sharing_counts.tsv')
    fig=plt.figure(figsize=(6.5,5.2))
    ymax=int(np.ceil(exact_counts.distinct_metabolites.max()*1.25/10)*10)
    for i,(condition,role) in enumerate(specs):
        row,col=divmod(i,2);ax=fig.add_axes([.105+col*.49,.575-row*.43,.355,.29])
        dat=exact_counts[(exact_counts.condition==condition)&(exact_counts.role==role)]
        color=RED if role=='substrates' else BLUE
        ax.bar(dat.n_strains,dat.distinct_metabolites,width=.82,color=color,linewidth=0)
        ax.set(xlim=(.3,21.7),ylim=(0,ymax),xlabel='Number of strains',ylabel='Number of distinct metabolites')
        ax.set_xticks([1,5,10,15,20]);ax.tick_params(labelsize=8)
        ax.set_title(f'{condition} {"uptake" if role=="substrates" else "secretion"}',loc='left',fontsize=10,pad=8)
        ax.text(.03,.95,f'{int(dat.distinct_metabolites.sum())} metabolites; {20 if condition=="Basal" else 21} strains',
                transform=ax.transAxes,va='top',fontsize=8)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        panel_label(fig,chr(65+i),.025+col*.49,.94-row*.43)
    fig.text(.5,.985,'How many strains can take up or secrete each metabolite?',ha='center',va='top',fontsize=10)
    fig.text(.5,.026,'Each metabolite is counted once per panel. Water is excluded.',ha='center',fontsize=8)
    sharing_stem=save(fig,'S09_sharing_distribution')

    fig=plt.figure(figsize=(6.5,6.0))
    panel_label(fig,'E',.025,.985)
    fig.text(.51,.98,'Metabolites available only after conditioning',ha='center',va='top',fontsize=11)
    fig.text(.51,.928,'Other strains predicted to secrete each metabolite in basal medium',ha='center',va='top',fontsize=9)
    ax2=fig.add_axes([.15,.15,.79,.66]);left=np.zeros(len(ORDER))
    supplier_cols=['#555555','#E69F00','#8FB5C9']
    for cat,col in zip(['0 other','1 other','2+ others'],supplier_cols):
        vals=ddf.set_index('strain').loc[ORDER,cat].values
        ax2.barh(range(len(ORDER)),vals,left=left,color=col,height=.73,linewidth=0)
        left+=vals
    for i,total in enumerate(left):ax2.text(total+1,i,str(int(total)),va='center',fontsize=8)
    ax2.set_yticks(range(len(ORDER)));ax2.set_yticklabels(ORDER,fontsize=8);ax2.invert_yaxis();color_ticks(ax2,'y')
    ax2.set_xlim(0,137);ax2.set_xlabel('Number of metabolites predicted to be taken up',fontsize=9)
    ax2.set_ylabel('Recipient strain',fontsize=9);ax2.tick_params(axis='x',labelsize=8)
    ax2.spines['top'].set_visible(False);ax2.spines['right'].set_visible(False)
    ax2.legend([Patch(facecolor=c) for c in supplier_cols],['0 suppliers','1 supplier','2+ suppliers'],
               loc='lower center',bbox_to_anchor=(.5,1.015),ncol=3,frameon=False,fontsize=9)
    fig.text(.5,.060,'Excludes the 145 basal-medium inputs and the recipient itself as a supplier.',ha='center',fontsize=8)
    fig.text(.5,.029,'The five metabolites with 0 other suppliers are self-only predictions for racd.',ha='center',fontsize=8)
    donor_stem=save(fig,'S09_other_potential_suppliers')
    for ext in ['png','pdf','svg']:
        old=FIG/f'S09_sharing_and_other_donors.{ext}'
        if old.exists():old.unlink()
    figures.append({'number':9,'title':'Shared metabolic capabilities and potential metabolite suppliers.',
      'caption':'(A-D) Numbers of distinct metabolites that can be taken up or secreted by a given number of strains, '
      'under basal or conditioned conditions. Each metabolite is counted once in each panel; water is excluded. '
      '(E) For each recipient, newly available metabolites are grouped by whether zero, one, or at least two other strains '
      'could secrete them in basal medium. The 145 basal inputs and the recipient itself are excluded from this supplier count. '
      'Zero suppliers denotes five self-only predictions for racd. Predictions combine alternatives across 10-100% of optimal biomass yield '
      'and do not establish metabolite transfer or growth dependence.',
      'pages':[sharing_stem,donor_stem],
      'page_descriptions':{sharing_stem:'A, basal uptake; B, conditioned uptake; C, basal secretion; D, conditioned secretion. Bars count distinct metabolites at each exact number of strains.',
         donor_stem:'E, potential suppliers of metabolites available only after conditioning. Each bar is one recipient; segments count metabolites with 0, 1, or 2+ other suppliers.'},
      'sources':[str(source_paths[(c,'pooled')].relative_to(ROOT)) for c in ['Basal','Conditioned']]+[str(pool_files[1].relative_to(ROOT))],
      'checks':{'other_donor_bins':[5,65,1826],'recipient_metabolite_memberships':1896,'water_excluded_from_sharing':True,
         'sharing_totals':{c+' '+r:sum(x['n_metabolites'] for x in sharing_counts if x['condition']==c and x['role']==r) for c,r in specs}}})

    # S10: exact historical dual-role exclusion; retain every eligible
    # metabolite, unlike the archived top-ten display cap.
    bprod=met_map(basal['products']);csub=met_map(cond['substrates']);bsub=union(basal['substrates'])
    edge_rows=[];eligibility=[]
    for met in sorted(set(bprod)|set(csub)):
        p0=bprod.get(met,set());c0=csub.get(met,set());dual=p0&c0;p=p0-dual;c=c0-dual
        eligible=met not in bsub and bool(p) and bool(c)
        eligibility.append({'metabolite':met,'name':names[met],'class':classes[met],
          'basal_substrate':met in bsub,'eligible':eligible,'original_producers':';'.join(sorted(p0)),
          'original_consumers':';'.join(sorted(c0)),'dual_role_removed':';'.join(sorted(dual)),
          'retained_producers':';'.join(s for s in ORDER if s in p),
          'retained_consumers':';'.join(s for s in ORDER if s in c)})
        if eligible:
            for producer in ORDER:
                if producer not in p:continue
                for consumer in ORDER:
                    if consumer in c:
                        assert producer!=consumer
                        edge_rows.append({'class':classes[met],'metabolite':met,'name':names[met],
                                          'producer':producer,'consumer':consumer})
    edges=pd.DataFrame(edge_rows);source_table(edges,'S10_all_eligible_exchange_triples.tsv')
    source_table(pd.DataFrame(eligibility),'S10_metabolite_eligibility.tsv')
    extra_classes=sorted(set(edges['class'])-MAIN_CLASSES,key=source_group_order);assert len(extra_classes)==14
    pages=[];page_meta=[];incidence=[];row_maps=[];page_descriptions={}
    display_names=dict(names)
    # The classification workbook gives the same generic label to two
    # different CMP isomers. Use the archived model's explicit 3'-CMP name.
    name_check_file=CODE/'Models/ThermobbetPTYGVitconstrained.mat'
    name_model=loadmat(name_check_file,simplify_cells=True)['ttmodel_bbet']
    name_index=list(name_model['mets']).index('C05822_c')
    assert name_model['metNames'][name_index]=="3'-CMP"
    display_names['C05822']='Cytidine 3′-monophosphate'
    class_letters={c:chr(65+i) for i,c in enumerate(extra_classes)}
    class_rows={}
    for cls in extra_classes:
        e=edges[edges['class']==cls]
        mets=sorted(set(e.metabolite),key=lambda m:(str(display_names[m]).lower(),m))
        class_rows[cls]=[]
        for m in mets:
            suppliers=[s for s in ORDER if s in set(e.loc[e.metabolite==m,'producer'])]
            recipients=[s for s in ORDER if s in set(e.loc[e.metabolite==m,'consumer'])]
            left_text=textwrap.fill(', '.join(suppliers),width=25,break_long_words=False)
            right_text=textwrap.fill(', '.join(recipients),width=25,break_long_words=False)
            met_text=textwrap.fill(str(display_names[m]),width=26,break_long_words=False)
            nlines=max(t.count('\n')+1 for t in [left_text,right_text,met_text])
            row_height=max(.44,.16*nlines+.14)
            class_rows[cls].append({'class':cls,'metabolite':m,'left':left_text,'middle':met_text,
                                   'right':right_text,'height':row_height,'triples':len(suppliers)*len(recipients)})
            row_maps.append({'class':cls,'panel':class_letters[cls],'metabolite':m,'source_name':names[m],
                'display_name':display_names[m],'producers':';'.join(suppliers),'consumers':';'.join(recipients),
                'n_producers':len(suppliers),'n_consumers':len(recipients),'n_triples':len(suppliers)*len(recipients)})
            for role,ss in [('producer',suppliers),('consumer',recipients)]:
                for s in ORDER:incidence.append({'class':cls,'metabolite':m,'role':role,'strain':s,'present':s in ss})
    row_map_df=pd.DataFrame(row_maps)
    assert len(row_map_df)==54 and row_map_df.n_triples.sum()==804
    source_table(row_map_df,'S10_readable_row_maps.tsv')
    source_table(pd.DataFrame(incidence),'S10_extra_class_incidence.tsv')

    # Paginate at readable type sizes, repeating class headers after a break.
    sheets=[];blocks=[];used=0.;header_height=.32;capacity=4.95
    for cls in extra_classes:
        for i,row in enumerate(class_rows[cls]):
            need_header=(i==0 or not blocks)
            need=row['height']+(header_height if need_header else 0)
            if blocks and used+need>capacity:
                sheets.append(blocks);blocks=[];used=0.;need_header=True
            if need_header:
                blocks.append({'header':True,'class':cls,'continued':i>0,'height':header_height});used+=header_height
            blocks.append(row);used+=row['height']
    if blocks:sheets.append(blocks)
    all_drawn=[]
    for page,blocks in enumerate(sheets,1):
        content_height=sum(b['height'] for b in blocks)
        height=content_height+1.13
        fig=plt.figure(figsize=(6.5,height));ax=fig.add_axes([0,0,1,1])
        ax.set(xlim=(0,6.5),ylim=(0,height));ax.axis('off')
        ax.text(.14,height-.12,f'Candidate metabolite exchanges · sheet {page} of {len(sheets)}',
                fontsize=10,va='top',weight='bold')
        header_y=height-.50
        for x,t,c in [(1.04,'Predicted secretion\nin basal medium',BLUE),(3.25,'Metabolite','#333333'),
                      (5.47,'Predicted uptake\nafter conditioning',RED)]:
            ax.text(x,header_y,t,ha='center',va='center',fontsize=9,color=c,weight='bold')
        for left,right in [(2.04,2.27),(4.24,4.47)]:
            ax.annotate('',xy=(right,header_y),xytext=(left,header_y),arrowprops={'arrowstyle':'->','color':'#555555','lw':.8})
        y=height-.85;nr=0
        for row in blocks:
            h=row['height']
            if row.get('header'):
                ax.add_patch(Rectangle((.10,y-h),6.30,h,facecolor='#E8ECEF',edgecolor='none'))
                suffix=' (continued)' if row['continued'] else ''
                ax.text(.18,y-h/2,f'{class_letters[row["class"]]}  {label_class(row["class"])}{suffix}',
                        fontsize=9,weight='bold',va='center')
            else:
                if nr%2==0:ax.add_patch(Rectangle((.10,y-h),6.30,h,facecolor='#F7F8F9',edgecolor='none'))
                yc=y-h/2
                ax.text(1.04,yc,row['left'],ha='center',va='center',fontsize=8.5,color=BLUE,linespacing=1.18)
                ax.text(3.25,yc,row['middle'],ha='center',va='center',fontsize=9,color='#222222',linespacing=1.18)
                ax.text(5.47,yc,row['right'],ha='center',va='center',fontsize=8.5,color=RED,linespacing=1.18)
                for left,right in [(2.04,2.27),(4.24,4.47)]:
                    ax.annotate('',xy=(right,yc),xytext=(left,yc),arrowprops={'arrowstyle':'->','color':'#777777','lw':.65})
                nr+=1;all_drawn.append(row['metabolite'])
            y-=h
        ax.text(3.25,.11,'Read each row from left to right. Lists show strain codes; arrows show predicted possibilities.',
                ha='center',va='center',fontsize=7.5)
        stem=save(fig,f'S10_exchange_classes_{page:02d}');pages.append(stem)
        sheet_groups=list(dict.fromkeys(b['class'] for b in blocks))
        page_descriptions[stem]='; '.join(f'{class_letters[c]}, {label_class(c)}' for c in sheet_groups)+'. Read producer strains → metabolite → consumer strains.'
        for cls in sheet_groups:
            rr=[b for b in blocks if not b.get('header') and b['class']==cls]
            page_meta.append({'page':page,'panel':class_letters[cls],'class':cls,'n_metabolites':len(rr),
                              'n_triples':sum(r['triples'] for r in rr),'stem':stem})
    assert len(all_drawn)==54 and len(set(all_drawn))==54
    source_table(pd.DataFrame(page_meta),'S10_page_index.tsv')
    for old_page in range(len(sheets)+1,15):
        for ext in ['png','pdf','svg']:
            old=FIG/f'S10_exchange_classes_{old_page:02d}.{ext}'
            if old.exists():old.unlink()
    figures.append({'number':10,'title':'Predicted producers and consumers of additional metabolite classes.',
      'caption':'Each row connects strains predicted to secrete a metabolite in basal medium (left, blue) to strains predicted to '
      'take it up after conditioning (right, red). Every listed producer can be paired with every listed consumer in that row. '
      'Panels A-N cover 54 metabolites in 14 additional classes; sugars, carboxylic acids, and vitamins and cofactors are shown in the main figure. '
      'Metabolites already present among basal uptake profiles are excluded. For each metabolite, strains predicted in both roles are '
      'removed from both lists, following the original map-selection rule. The maps preserve 804 candidate links across pooled alternatives '
      'and identify hypotheses for exchange, not observed transfer. Full metabolite identifiers and link lists accompany the figure.',
      'pages':pages,'page_descriptions':page_descriptions,
      'sources':[str(source_paths[(c,'pooled')].relative_to(ROOT)) for c in ['Basal','Conditioned']]+[
          str(class_file.relative_to(ROOT)),'plotting/main/DiMEs/Uptake_SecretionProfile.ipynb',str(name_check_file.relative_to(ROOT))],
      'checks':{'extra_classes':extra_classes,'page_index':page_meta,'all_classes_triples':len(edges),
          'extra_classes_metabolites':int(edges[~edges['class'].isin(MAIN_CLASSES)].metabolite.nunique()),
          'extra_classes_triples':804,'all54_metabolites_drawn_once':True,
          'dual_role_exclusion_verified':True,'self_edges':int((edges.producer==edges.consumer).sum()),'top_n_truncation':False}})
    audit.append({'display_name_clarification':{'metabolite':'C05822','source_label':names['C05822'],
       'figure_label':display_names['C05822'],'evidence':'Archived bbet model labels C05822_c as 3-prime CMP.',
       'source':str(name_check_file.relative_to(ROOT)),'sha256':sha256(name_check_file)}})

    # S11: use archived optimum and pooled summaries only. Raw per-yield
    # workspaces are intentionally not joined because their provenance differs.
    comparison=[]
    for condition in ['Basal','Conditioned']:
        for role in ['substrates','products']:
            o=summaries[(condition,'optimal')][role];u=summaries[(condition,'pooled')][role]
            for s in ORDER:
                if s not in u:continue
                assert o[s]<=u[s],(condition,role,s)
                comparison.append({'condition':condition,'role':role,'strain':s,
                 'optimal_count':len(o[s]),'pooled_count':len(u[s]),'additional_count':len(u[s]-o[s]),
                 'optimal_fraction_of_pooled':len(o[s])/len(u[s]) if u[s] else np.nan,
                 'additional_metabolites':';'.join(sorted(u[s]-o[s]))})
    comparison=pd.DataFrame(comparison);source_table(comparison,'S11_optimal_vs_pooled.tsv')
    fig=plt.figure(figsize=(6.5,6.1))
    for i,(condition,role) in enumerate([(c,r) for c in ['Basal','Conditioned'] for r in ['substrates','products']]):
        col=i%2;row=i//2;ax=fig.add_axes([.105+col*.48,.555-row*.455,.365,.345])
        color=RED if role=='substrates' else BLUE
        dat=comparison[(comparison.condition==condition)&(comparison.role==role)].set_index('strain')
        for j,s in enumerate(ORDER):
            if s not in dat.index:
                ax.axhspan(j-.45,j+.45,color='#E8E8E8',zorder=0);continue
            a,b=dat.loc[s,['optimal_count','pooled_count']]
            ax.plot([a,b],[j,j],color='#B3B3B3',lw=.7,zorder=1)
            ax.plot(b,j,'o',ms=3.1,color=color,zorder=2)
            ax.plot(a,j,'o',ms=3.1,mfc='white',mec='#333333',mew=.7,zorder=3)
        ax.set_yticks(range(len(ORDER)));ax.set_yticklabels(ORDER,fontsize=7);color_ticks(ax,'y');ax.set_ylim(len(ORDER)-.4,-.6)
        mx=comparison[comparison.role==role].pooled_count.max();ax.set_xlim(-3,mx*1.04)
        ax.set_title(f'{condition}: {"uptake" if role=="substrates" else "secretion"}',loc='left',fontsize=9,pad=7)
        ax.tick_params(axis='y',length=0,pad=3);ax.tick_params(axis='x',labelsize=7)
        ax.set_xlabel('Number of feasible metabolites',fontsize=8)
        panel_label(fig,chr(65+i),.02+col*.48,.943-row*.455)
    fig.legend([Line2D([],[],marker='o',mfc='white',mec='#333333',color='none',ms=4),
                Line2D([],[],marker='o',color='#555555',lw=0,ms=4)],
               ['Optimal yield (100%)','Pooled yields (10-100%)'],ncol=2,loc='upper center',bbox_to_anchor=(.5,1),frameon=False,fontsize=8)
    fig.text(.5,.019,'Paired model predictions; grey marks the unavailable basal mip profile.',ha='center',fontsize=7)
    stem=save(fig,'S11_optimal_vs_pooled')
    # Directional overlap uses the same 20 strains and 380 off-diagonal
    # directed pairs in every condition/selection; no inferential test is used.
    common=[s for s in ORDER if all(s in d['substrates'] for d in summaries.values())]
    assert len(common)==20 and 'mip' not in common
    noi_rows=[]
    for (condition,selection),data in summaries.items():
        p=data['substrates']
        for a in common:
            for b in common:
                if a==b:continue
                shared=len(p[a]&p[b]);n=len(p[a]);value=shared/n
                noi_rows.append({'condition':condition,'selection':selection,'focal_strain':a,'other_strain':b,
                   'intersection_count':shared,'focal_uptake_count':n,'directional_NOI':value,
                   'NOI_at_least_0_75':value>=.75})
    noi=pd.DataFrame(noi_rows);source_table(noi,'S11_directional_NOI.tsv')
    noi_wide=noi.pivot(index=['focal_strain','other_strain'],columns=['condition','selection'],values='directional_NOI')
    threshold=noi.groupby(['selection','condition'],sort=False).NOI_at_least_0_75.sum().astype(int)
    assert [threshold[x] for x in [('optimal','Basal'),('optimal','Conditioned'),
                                    ('pooled','Basal'),('pooled','Conditioned')]]==[132,89,343,192]
    source_table(threshold.rename('n_pairs').reset_index().assign(threshold=.75,total_directed_pairs=380),
                 'S11_NOI_threshold_counts.tsv')
    fig=plt.figure(figsize=(6.5,5.7))
    for k,condition in enumerate(['Basal','Conditioned']):
        ax=fig.add_axes([.10+k*.48,.525,.365,.37])
        ax.plot([0,1],[0,1],color='#AAAAAA',lw=.8,zorder=0)
        ax.scatter(noi_wide[(condition,'optimal')],noi_wide[(condition,'pooled')],
                   s=10,color=BLUE if k==0 else RED,alpha=.4,edgecolors='none')
        ax.set(xlim=(-.025,1.025),ylim=(-.025,1.025),xlabel='Optimal-yield directional NOI',
               ylabel='Pooled-yield directional NOI')
        ax.set_xticks([0,.25,.5,.75,1]);ax.set_yticks([0,.25,.5,.75,1])
        ax.set_title(condition,loc='left',pad=8)
        panel_label(fig,chr(69+k),.02+k*.48,.965)
    ax=fig.add_axes([.10,.10,.82,.285])
    x=np.arange(2)
    for selection,color,marker in [('optimal','#444444','o'),('pooled',BLUE,'s')]:
        vals=[threshold[(selection,c)] for c in ['Basal','Conditioned']]
        ax.plot(x,vals,color=color,marker=marker,ms=4,lw=1.2,
                label='Optimal yield (100%)' if selection=='optimal' else 'Pooled yields (10-100%)')
        for xx,v in zip(x,vals):ax.annotate(str(v),(xx,v),xytext=(7,5),textcoords='offset points',fontsize=8,color=color)
    ax.set_xticks(x);ax.set_xticklabels(['Basal','Conditioned']);ax.set_xlim(-.12,1.12);ax.set_ylim(0,400)
    ax.set_yticks([0,100,200,300,380]);ax.set_ylabel(r'Directed pairs with NOI $\geq$ 0.75')
    ax.legend(frameon=False,loc='upper right',fontsize=8)
    panel_label(fig,'G',.02,.415)
    fig.text(.5,.022,'Matched 20 strains; 380 directed pairs in every comparison. Water is retained.',ha='center',fontsize=7)
    noi_stem=save(fig,'S11_niche_overlap_sensitivity')
    figures.append({'number':11,'title':'Effect of biomass yield on predicted exchange and niche overlap.',
      'caption':'Paired counts compare the archived optimal-yield (100%; open circles) and pooled-yield (10-100%; filled circles) DiMEs '
      'summaries for uptake and secretion in basal (A,B) and conditioned (C,D) environments. Lines connect the two summaries for each strain. '
      'Each count is the number of distinct metabolites appearing in at least one retained profile; water is retained, matching the archived '
      'breadth summaries. All optimal-yield metabolite sets are subsets of their corresponding pooled sets. Basal mip has no feasible '
      'profile and is shown in gray. (E,F) Directional niche overlap, NOI(i,j) = |U_i ∩ U_j|/|U_i|, compares optimal and pooled uptake '
      'sets for the same 380 directed pairs among the 20 strains present in both environments; mip and self-pairs are excluded. '
      'Each point is one directed pair, with transparency showing overplotting; gray lines mark equality. (G) Numbers of pairs with '
      'NOI ≥ 0.75 fall from 132 to 89 between basal and conditioned environments at optimal yield and from 343 to 192 for pooled yields. '
      'Water is also retained in E-G, matching the archived NOI calculation. These paired comparisons describe model capabilities; '
      'directed pairs are not independent biological replicates. Pooled alternatives need not be realized simultaneously. '
      'The figure does not establish exhaustive enumeration, solution probabilities, or iteration convergence.',
      'pages':[stem,noi_stem],'sources':[str(p.relative_to(ROOT)) for p in source_paths.values()]+['plotting/main/DiMEs/NOI.ipynb'],
      'page_descriptions':{stem:'A, basal uptake; B, basal secretion; C, conditioned uptake; D, conditioned secretion. Open symbols denote 100% yield; filled symbols combine 10-100% yield predictions.',
         noi_stem:'E, basal niche overlap; F, conditioned niche overlap; G, number of directed pairs with NOI at least 0.75. All comparisons use the same 20 strains and 380 directed pairs.'},
      'checks':{'paired_strain_role_conditions':len(comparison),'optimal_sets_are_subsets':True,
                'water_retained':True,'raw_per_yield_files_used':False,'NOI_common_strains':20,'NOI_directed_pairs':380,
                'NOI_threshold':.75,'NOI_threshold_optimal':[132,89],'NOI_threshold_pooled':[343,192]}})
    audit.append({'S8_normalization':'A common 267-compound catalog replaces the role/environment-specific denominator in archived '
       'Uptake_SecretionProfile.ipynb cells 10 and 11. The catalog equals the union of all four pooled profiles, excluding water.'})
    audit.append({'note':'Per-yield provenance conflicts are documented rather than used to recompute manuscript-aligned pooled summaries.',
      'conflicts':['basal cuvk 70/80% raw-vs-pooled profiles','conditioned rfa raw-vs-pooled profiles'],
      'missing':['basal flk 80%','conditioned mlot 20% Active/CSV','conditioned mip per-yield archives'],
      'enumeration_caution':'600 alternatives is a cap, not proof of exhaustive coverage.'})
    audit.append({'input_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in [*source_paths.values(),class_file,*pool_files]}})
    manifest('dimes',figures,audit)
    print('Created',[(f['number'],len(f['pages'])) for f in figures])


if __name__=='__main__':
    main()
