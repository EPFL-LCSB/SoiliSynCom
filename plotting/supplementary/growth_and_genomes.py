"""Plot growth and genome panels underlying current Supplementary V5.

No growth fitting, imputation, outlier exclusion, or model simulation is performed.
"""
from datetime import time, timedelta, datetime
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from openpyxl import load_workbook
from openpyxl.utils.datetime import to_excel
from common import ROOT, EXTERNAL_DATA_ROOT, DATA, ORDER, COLORS, GROUPS, GROUP_COLORS, save, manifest, source_table, sha256

IDS = {'6146':'vpe','6147':'flk','6150':'mip','6152':'mgk','6155':'cuvk',
       '6477':'rfa','6478':'bbet','6483':'cxy','6500':'lur','6508':'pzu',
       '6519':'mlot','6525':'lgu','6533':'pscf','6537':'cauj','6551':'pmor',
       '6593':'dri','6695':'cohm','6697':'tro','6699':'buol','6700':'racd',
       'spr':'cpi','spread':'cpi','Pv':'psvr'}
INV_IDS = {v:k for k,v in IDS.items() if k not in ['spread','Pv']}
GENUS = {'bbet':'Bradyrhizobium','cauj':'Caulobacter','dri':'Devosia','mlot':'Mesorhizobium',
 'pzu':'Phenylobacterium','tro':'Tardiphaga','buol':'Burkholderia','vpe':'Variovorax',
 'lgu':'Lysobacter','lur':'Luteibacter','pmor':'Pseudomonas 1','pscf':'Pseudomonas 2',
 'racd':'Rahnella','cuvk':'Curtobacterium','cxy':'Cellulomonas','mip':'Microbacterium',
 'rfa':'Rhodococcus','cpi':'Chitinophaga','flk':'Flavobacterium','mgk':'Mucilaginibacter','cohm':'Cohnella'}
SUCC = ROOT/'tables/graphing/ExperimentalGrowthData/Senka_OD600_succinategrowth_pooleddata.xlsx'
PTYG = ROOT/"tables/graphing/ExperimentalGrowthData/PTYG_OD600_21strain_Senka&Tania.xlsx"
RESCUE = ROOT/'tables/graphing/ExperimentalValidation/MicroCohnella_MinRichMedia_0924.xlsx'
GENOMES = EXTERNAL_DATA_ROOT/'Code/DraftModelReconstruction/Genomes'
BATCH_STYLE = {'part1':'-', 'part2':(0,(5,2)), 'selected_strains_repeat':(0,(1.2,1.7))}
BATCH_NAME = {'part1':'Original batch 1','part2':'Original batch 2','selected_strains_repeat':'Repeat batch'}


def elapsed_hours(x):
    if isinstance(x,timedelta):
        return x.total_seconds()/3600
    if isinstance(x,time):
        return x.hour+x.minute/60+x.second/3600+x.microsecond/3.6e9
    if isinstance(x,datetime):
        return to_excel(x)*24
    return float(x)*24


def clean_axes(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(length=2.6,pad=2)


def atlas(frame, number, medium, batches, xmax, ymax):
    pages=[]
    # Twelve panels on page one; nine on page two. Panel ordering is shared throughout V4.
    for page, strains in enumerate([ORDER[:12],ORDER[12:]],1):
        nrows=math.ceil(len(strains)/3)
        fig,axes=plt.subplots(nrows,3,figsize=(6.5,6.15 if nrows==4 else 4.95),squeeze=False)
        fig.subplots_adjust(left=.092,right=.98,bottom=.135,top=.90,wspace=.30,hspace=.48)
        for i,(ax,strain) in enumerate(zip(axes.flat,strains)):
            sub=frame[(frame.model_id==strain)&(frame.role=='inoculated')]
            for batch in batches:
                s=sub[sub.sheet==batch]
                if s.empty: continue
                for _,rep in s.groupby('sample_id',sort=False):
                    ax.plot(rep.time_h,rep.corrected_od600,color=COLORS[strain],
                            linestyle=BATCH_STYLE.get(batch,'-'),lw=.48,alpha=.42)
                mean=s.groupby('time_h',sort=True).corrected_od600.mean()
                ax.plot(mean.index,mean.values,color=COLORS[strain],
                        linestyle=BATCH_STYLE.get(batch,'-'),lw=1.25)
            ax.set_title(f'{strain}  |  {GENUS[strain]}',loc='left',fontsize=8,pad=4)
            ax.set_xlim(0,xmax);ax.set_ylim(-.025,ymax)
            ax.set_xticks([0,24,48,72,90] if number==1 else [0,24,48,72,96,116])
            ax.set_yticks([0,.5,1.0] if number==1 else [0,.5,1.0])
            ax.axhline(0,color='#b8b8b8',lw=.4,zorder=0)
            clean_axes(ax)
        fig.suptitle(f'{medium} growth trajectories'+(' (continued)' if page==2 else ''),fontsize=11,y=.973)
        fig.text(.016,.51,'Blank-corrected OD$_{600}$',rotation=90,va='center',fontsize=9)
        fig.text(.53,.074,'Time (h)',ha='center',fontsize=9)
        if number==1:
            handles=[Line2D([0],[0],color='#333333',lw=1.25,ls=BATCH_STYLE[b],label=BATCH_NAME[b]) for b in batches]
            fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.52,.002),ncol=3,frameon=False)
        else:
            fig.text(.52,.027,'Thin lines: individual wells (n = 4); thick line: mean',ha='center',fontsize=7)
        pages.append(save(fig,f'S{number:02d}_{"succinate" if number==1 else "PTYG"}_atlas_p{page}'))
    return pages


def succinate():
    w=load_workbook(SUCC,read_only=True,data_only=True)
    rows=[]; checks={}
    for sheet in BATCH_STYLE:
        a=list(w[sheet].values); headers=list(a[0]); values=np.array([r[1:] for r in a[1:]],float)
        t=np.array([elapsed_hours(r[0]) for r in a[1:]])
        assert np.all(np.diff(t)>0) and np.isfinite(values).all()
        blankcols=[i for i,h in enumerate(headers[1:]) if h.startswith('MM_')]
        assert len(blankcols)==4
        blank=values[:,blankcols].mean(axis=1)
        for j,sample in enumerate(headers[1:]):
            isolate,rep=sample.rsplit('_',1); role='blank' if isolate=='MM' else 'inoculated'
            model=IDS.get(isolate,'')
            for k,th in enumerate(t):
                rows.append({'source_file':str(SUCC.relative_to(ROOT)),'sheet':sheet,
                 'source_row':k+2,'source_column':j+2,'sample_id':sample,'isolate_id':isolate,
                 'model_id':model,'role':role,'replicate':int(rep),'time_h':th,
                 'raw_od600':values[k,j],'blank_mean_od600':blank[k],
                 'corrected_od600':values[k,j]-blank[k] if role=='inoculated' else np.nan,
                 'plotted':model in ORDER})
        checks[sheet]={'timepoints':len(t),'start_h':float(t[0]),'end_h':float(t[-1]),
                       'blank_wells':4,'inoculated_wells':len(headers)-5}
    d=pd.DataFrame(rows)
    assert set(d[d.plotted].model_id)==set(ORDER)
    assert d[d.plotted].groupby(['sheet','model_id']).sample_id.nunique().eq(4).all()
    source_table(d,'S01_succinate_all_measurements.tsv')
    means=d[d.plotted].groupby(['sheet','model_id','time_h'],as_index=False).agg(
        mean_corrected_od600=('corrected_od600','mean'),sd_corrected_od600=('corrected_od600','std'),n=('corrected_od600','count'))
    source_table(means,'S01_succinate_means.tsv')
    checks['sha256']=sha256(SUCC)
    return {'number':1,'title':'Replicate growth trajectories on succinate',
      'caption':'Optical density at 600 nm (OD600) of the 21 SynCom isolates in succinate medium. '
      'Thin lines show individual inoculated wells and thick lines show the mean of four wells within each experimental batch. '
      'At every time point, the mean of the four uninoculated medium wells from the same batch was subtracted from each inoculated well. '
      'Original batches 1 and 2 are shown by solid and long-dashed lines, respectively; the repeat batch is dotted. '
      'Ten isolates were measured again in the repeat batch; pmor was present in both original batches. Batches are displayed separately. '
      'All recorded inoculated wells are retained, without smoothing or baseline normalization. '
      'The original batches extend to 72.07 and 69.52 h and the repeat batch to 90.32 h; lines stop at the final measurement. '
      'Axes use the same limits across isolates. Isolates are ordered by their taxonomic groups, with matching group colors. '
      'Model codes are followed by the genus labels used in the original culture records. '
      'The two pages together show all 21 isolates.',
      'pages':atlas(d,1,'Succinate',list(BATCH_STYLE),91,1.10),'sources':[str(SUCC.relative_to(ROOT)),
       'Code/DraftModelReconstruction/Genomes/GenomeInformation_Senka_summary.xlsx'],'checks':checks}


def ptyg():
    w=load_workbook(PTYG,read_only=True,data_only=True)
    a=list(w['Sheet1'].values);headers=list(a[0]);values=np.array([r[1:] for r in a[1:]],float)
    t=np.array([elapsed_hours(r[0]) for r in a[1:]])
    assert values.shape==(232,96) and np.isfinite(values).all() and np.all(np.diff(t)>0)
    wellmap={r[4]:r[5] for r in list(w['Sheet2'].values)[1:] if r[4]}
    rows=[];corrected=[]
    for j,well in enumerate(headers[1:]):
        label=wellmap[well]; isblank=label=='#N/A'
        model='' if isblank else IDS[label.rsplit('_',1)[1]]
        blankrow='E' if well[0] in 'ABCDE' else well[0]
        blankwells=[f'{blankrow}{i}' for i in [10,11,12]]
        blank=values[:,[headers.index(b)-1 for b in blankwells]].mean(axis=1)
        y=values[:,j]-blank
        if not isblank: corrected.append(y)
        for k,th in enumerate(t):
            rows.append({'source_file':str(PTYG.relative_to(ROOT)),'sheet':'Sheet1',
             'source_row':k+2,'source_column':j+2,'sample_id':well,'source_label':label,
             'model_id':model,'role':'blank' if isblank else 'inoculated','time_h':th,
             'raw_od600':values[k,j],'blank_wells':','.join(blankwells),'blank_mean_od600':blank[k],
             'corrected_od600':y[k] if not isblank else np.nan,'plotted':not isblank})
    # Independent confirmation against all existing corrected numeric cells.
    existing=np.array(list(w['Data'].values)[1:],float)
    diff=float(np.max(np.abs(np.array(corrected).T-existing[:,1:])))
    assert diff<1e-12 and np.allclose(t,existing[:,0],atol=1e-12)
    d=pd.DataFrame(rows);assert d[d.plotted].groupby('model_id').sample_id.nunique().eq(4).all()
    assert set(d[d.plotted].model_id)==set(ORDER)
    source_table(d,'S02_PTYG_all_measurements.tsv')
    means=d[d.plotted].groupby(['model_id','time_h'],as_index=False).agg(
        mean_corrected_od600=('corrected_od600','mean'),sd_corrected_od600=('corrected_od600','std'),n=('corrected_od600','count'))
    source_table(means,'S02_PTYG_means.tsv')
    return {'number':2,'title':'Replicate growth trajectories in PTYG',
      'caption':'Optical density at 600 nm (OD600) of the 21 SynCom isolates in PTYG medium. '
      'Thin lines represent four individual inoculated wells per isolate and thick lines their mean. '
      'Blank correction follows the recorded plate scheme: wells in rows A–D and E1–E9 are corrected using the mean of E10–E12; '
      'F1–F9, G1–G9 and H1–H9 use the three uninoculated wells in columns 10–12 of their respective rows. '
      'All wells, including wells at the plate edge, are retained without smoothing, outlier exclusion or baseline normalization. '
      'Measurements span 0.50–116.00 h. Axes use the same limits across isolates; ordering and colors follow the same taxonomic groups '
      'as in Supplementary Fig. S1. Model codes are followed by the genus labels used in the original culture records. '
      'The two pages together show all 21 isolates.',
      'pages':atlas(d,2,'PTYG',['Sheet1'],116,1.35),'sources':[str(PTYG.relative_to(ROOT)),
       'Code/DraftModelReconstruction/Genomes/GenomeInformation_Senka_summary.xlsx'],
      'checks':{'timepoints':len(t),'isolates':21,'inoculated_wells':84,'blank_wells':12,
       'start_h':float(t[0]),'end_h':float(t[-1]),'max_difference_from_recorded_Data_sheet':diff,'sha256':sha256(PTYG)}}


def supplementation():
    w=load_workbook(RESCUE,read_only=True,data_only=True);a=list(w['Sheet1'].values)
    headers=list(a[0]);v=np.array(a[1:],float);t=v[:,0]
    assert np.array_equal(t,np.arange(165))
    rows=[]
    for model,species in [('mip','Microbacterium'),('cohm','Cohnella')]:
        for medium in ['MinMedia','RichMedia']:
            bi=[i for i,h in enumerate(headers) if h.startswith(f'{medium}_blank_{species}_')]
            assert len(bi)==4 and np.isfinite(v[:,bi]).all()
            blank=v[:,bi].mean(axis=1)
            for role in ['blank','inoc']:
                indices=[i for i,h in enumerate(headers) if h.startswith(f'{medium}_{role}_{species}_')]
                assert len(indices)==4
                for j in indices:
                    corrected=v[:,j]-blank;delta=corrected-corrected[0]
                    for k,th in enumerate(t):
                        rows.append({'source_file':str(RESCUE.relative_to(ROOT)),'sheet':'Sheet1',
                         'source_row':k+2,'source_column':j+1,'sample_id':headers[j],
                         'model_id':model,'species_label':species,'medium':medium,'role':role,
                         'replicate':int(headers[j].rsplit('_',1)[1]),'time_h':th,'raw_od600':v[k,j],
                         'blank_mean_od600':blank[k],'corrected_od600':corrected[k],
                         'baseline_corrected_od600':corrected[0],
                         'delta_od600':delta[k] if role=='inoc' else np.nan,'plotted':role=='inoc'})
    d=pd.DataFrame(rows); inoc=d[d.plotted]
    missing=inoc[inoc.raw_od600.isna()]
    assert len(missing)==1 and missing.iloc[0].sample_id=='RichMedia_inoc_Cohnella_4' and missing.iloc[0].time_h==164
    endpoint=inoc[inoc.time_h==150].copy()
    stats=endpoint.groupby(['model_id','medium'],sort=False).delta_od600.agg(['mean','std','count'])
    assert np.isclose(stats.loc[('mip','MinMedia'),'mean'],.028) and np.isclose(stats.loc[('mip','MinMedia'),'std'],.0034156502553)
    source_table(d,'S07_supplementation_all_measurements.tsv')
    source_table(endpoint,'S07_supplementation_150h_replicates.tsv')
    source_table(stats.reset_index().rename(columns={'mean':'mean_delta_od600','std':'sd_delta_od600','count':'n'}),'S07_supplementation_150h_summary.tsv')
    means=inoc.groupby(['model_id','medium','time_h'],as_index=False).agg(mean_delta_od600=('delta_od600','mean'),sd_delta_od600=('delta_od600','std'),n=('delta_od600','count'))
    source_table(means,'S07_supplementation_means.tsv')
    medcol={'MinMedia':'#4B7FAD','RichMedia':'#9E4A4F'};medstyle={'MinMedia':'-','RichMedia':(0,(5,2))}
    fig,axes=plt.subplots(1,2,figsize=(6.5,3.2))
    fig.subplots_adjust(left=.12,right=.975,top=.82,bottom=.18,wspace=.30)
    for col,model in enumerate(['mip','cohm']):
        ax=axes[col]
        for medium in medcol:
            sub=inoc[(inoc.model_id==model)&(inoc.medium==medium)]
            for _,rep in sub.groupby('sample_id',sort=False):
                ax.plot(rep.time_h,rep.delta_od600,c=medcol[medium],ls=medstyle[medium],lw=.55,alpha=.4)
            mean=sub.groupby('time_h').delta_od600.mean()
            ax.plot(mean.index,mean.values,c=medcol[medium],ls=medstyle[medium],lw=1.35)
        ax.axvline(150,color='#999999',ls=':',lw=.75);ax.axhline(0,color='#aaaaaa',lw=.6)
        ax.set(xlim=(0,165),ylim=(-.041,.071),xticks=[0,50,100,150],xlabel='Time (h)')
        ax.set_title(f'{GENUS[model]} ({model})',fontsize=10,pad=7)
        ax.set_ylabel('Change in OD$_{600}$')
        ax.text(-.15,1.10,'AB'[col],transform=ax.transAxes,fontweight='bold',fontsize=11)
        clean_axes(ax)
    fig.legend(handles=[Line2D([0],[0],c=medcol[m],ls=medstyle[m],lw=1.4,label=m) for m in medcol],
               loc='upper center',bbox_to_anchor=(.56,.997),ncol=2,frameon=False)
    stem=save(fig,'S07_supplementation_validation')
    return {'number':7,'title':'Measured growth in the two supplementation conditions',
     'caption':'Growth of Microbacterium (mip) and Cohnella (cohm) in the two conditions recorded as MinMedia and RichMedia. '
     '(A, B) Individual inoculated-well trajectories (thin lines) and their means (thick lines). At each time point, the mean of four '
     'matching uninoculated medium wells was subtracted, followed by subtraction of each inoculated well’s corrected value at time zero. '
     'Thus, the ordinate is the change in blank-corrected OD600 from the initial measurement. Solid blue and dashed red lines identify '
     'MinMedia and RichMedia, respectively; the vertical dotted line marks 150 h. '
     'Four inoculated wells were measured per condition. The individual endpoint values and mean ± sample standard deviation at 150 h are retained in the source tables. All recorded measurements are retained, '
     'including negative changes. One cohm RichMedia measurement is missing at 164 h, where the mean uses the three available wells. '
     'The medium names reproduce the experimental record. No condition is identified as an unsupplemented control in the experimental record.',
     'pages':[stem],'sources':[str(RESCUE.relative_to(ROOT))],
     'checks':{'timepoints':165,'blank_wells':16,'inoculated_wells':16,'n_at_150h_each_condition':4,
       'missing_measurement':{'sample_id':'RichMedia_inoc_Cohnella_4','time_h':164,'source_cell':'AH166'},
       'summary_150h':stats.reset_index().to_dict('records'),'sha256':sha256(RESCUE)}}


def read_fasta(path):
    entries=[];header=None;seq=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if header is not None:entries.append((header,''.join(seq).upper()))
            header=line[1:];seq=[]
        elif line.strip():seq.append(line.strip())
    if header is not None:entries.append((header,''.join(seq).upper()))
    return entries


def assemblies():
    contigs=[];summary=[];sources=[]
    for model in ORDER:
        p=GENOMES/f'{model}.fasta';sources.append(str(p.relative_to(ROOT)))
        entries=read_fasta(p);lens=sorted([len(s) for _,s in entries],reverse=True)
        total=sum(lens);n50=next(length for length,acc in zip(lens,np.cumsum(lens)) if acc>=total/2)
        acgt=sum(sum(s.count(b) for b in 'ACGT') for _,s in entries)
        gc=sum(s.count('G')+s.count('C') for _,s in entries);ns=sum(s.count('N') for _,s in entries)
        assert acgt+ns==total
        summary.append({'model_id':model,'source_fasta':str(p.relative_to(ROOT)),'total_bp':total,
          'contig_count':len(lens),'n50_bp':n50,'gc_percent':100*gc/acgt,'n_bases':ns,'sha256':sha256(p)})
        for rank,(header,seq) in enumerate(sorted(entries,key=lambda x:len(x[1]),reverse=True),1):
            contigs.append({'model_id':model,'source_fasta':str(p.relative_to(ROOT)),
             'contig_rank_by_length':rank,'contig_header':header,'length_bp':len(seq),
             'gc_percent':100*(seq.count('G')+seq.count('C'))/len(seq),'n_bases':seq.count('N')})
    summ=pd.DataFrame(summary);con=pd.DataFrame(contigs)
    assert summ.contig_count.value_counts().to_dict()=={1:13,2:6,3:2}
    assert len(con)==31 and summ.n_bases.sum()==0
    source_table(summ,'S14_final_assembly_summary.tsv');source_table(con,'S14_final_assembly_contigs.tsv')
    fig,axes=plt.subplots(1,3,figsize=(6.5,6.3),sharey=True,gridspec_kw={'width_ratios':[3.4,1.0,1.4]})
    fig.subplots_adjust(left=.12,right=.965,bottom=.14,top=.90,wspace=.29)
    y=np.arange(len(ORDER))
    for i,model in enumerate(ORDER):
        left=0
        for row in con[con.model_id==model].itertuples():
            rank=row.contig_rank_by_length; length=row.length_bp/1e6
            axes[0].barh(i,length,left=left,height=.64,color=COLORS[model],
              alpha=[1,.63,.35][rank-1],edgecolor='white',linewidth=.55,hatch=['','///','...'][rank-1])
            left+=length
        sr=summ[summ.model_id==model].iloc[0]
        axes[0].scatter(sr.n50_bp/1e6,i,marker='D',s=16,color='#202020',
            edgecolor='white',linewidth=.35,zorder=5)
        axes[1].scatter(sr.contig_count,i,color=COLORS[model],s=24,edgecolor='#333333',linewidth=.3,zorder=3)
        axes[2].scatter(sr.gc_percent,i,color=COLORS[model],s=24,edgecolor='#333333',linewidth=.3,zorder=3)
    axes[0].set_yticks(y,ORDER);axes[0].invert_yaxis()
    axes[0].set(xlim=(0,9.2),xticks=[0,3,6,9],xlabel='Assembly length (Mb)')
    axes[1].set(xlim=(.65,3.35),xticks=[1,2,3],xlabel='Contigs')
    axes[2].set(xlim=(30,76),xticks=[35,50,65],xlabel='GC (%)')
    for ax,letter,title in zip(axes,'ABC',['Contig lengths','Contiguity','GC content']):
        ax.set_title(title,fontsize=9,pad=11)
        ax.text(-.13,1.055,letter,transform=ax.transAxes,fontweight='bold',fontsize=11)
        ax.set_ylim(20.7,-.7);clean_axes(ax)
        ax.grid(axis='x',color='#e5e5e5',lw=.4)
    for ax in axes[1:]: ax.tick_params(axis='y',left=False,labelleft=False);ax.spines['left'].set_visible(False)
    contig_handles=[Patch(facecolor='#555555',alpha=a,hatch=h,label=f'Contig {r}') for r,a,h in [(1,1,''),(2,.63,'///'),(3,.35,'...')]]
    contig_handles.append(Line2D([0],[0],marker='D',markersize=4,color='#202020',linestyle='none',label='N50'))
    fig.legend(handles=contig_handles,
         loc='lower center',bbox_to_anchor=(.51,.016),ncol=4,frameon=False)
    fig.text(.51,.963,'21 assemblies · 31 contigs',ha='center',fontsize=10)
    stem=save(fig,'S14_genome_assembly_structure')
    return {'number':14,'title':'Structure and GC content of the genome assemblies used for reconstruction',
     'caption':'Genome assembly characteristics calculated from the final nucleotide FASTA files for the 21 SynCom isolates. '
     '(A) Total assembly length, partitioned into individual contigs ordered from longest to shortest. Contigs are drawn to scale; '
     'hatching distinguishes the second and third contigs where present. Black diamonds mark N50, the contig length at which '
     'contigs ordered from longest to shortest reach at least half the assembly length. (B) Number of contigs per isolate. '
     '(C) GC content as a percentage of A, C, G and T bases. The collection contains 13 single-contig, six two-contig and two '
     'three-contig assemblies, comprising 31 contigs in total. No N bases occur in these sequence files. '
     'Strain ordering and colors follow the taxonomic grouping used throughout the supplementary figures. '
     'These are assembly-structure metrics and do not constitute independent estimates of genome completeness or contamination.',
     'pages':[stem],'sources':sources,'checks':{'assemblies':21,'contigs':31,'contig_count_distribution':{1:13,2:6,3:2},'n_bases':0}}


if __name__=='__main__':
    source_table(pd.DataFrame([{'model_id':m,'original_culture_genus_label':GENUS[m],
        'original_isolate_id':INV_IDS[m],'plot_group':GROUPS[m],'plot_color':COLORS[m]} for m in ORDER]),
        'experimental_strain_key.tsv')
    figures=[succinate(),ptyg(),supplementation(),assemblies()]
    manifest('experimental',figures,audit=[
      'S1: The original part1 workbook includes P. veronii (Pv), which is outside the 21-isolate collection. Its raw values are preserved in the source TSV with plotted=False; batch-specific blank wells are also preserved.',
      'S1: pmor has two original batches and ten other isolates have original/repeat batches; no batches or distinct replicate sets were pooled.',
      'S2: The blank correction exactly reproduces all 19,488 cells in the existing Data sheet (within floating-point precision). Plate-edge wells are retained; the source README warns about possible evaporation. The old thesis Caulobacter n=2 selection was not used.',
      'S7: MinMedia/RichMedia labels could not be linked unambiguously to recipes. No assertion that RichMedia is PTYG or that either condition tests individual-vitamin auxotrophy is made.',
      'S7: One previously unreported empty cell AH166 is the final 164-h cohm RichMedia replicate-4 reading. It remains missing; no interpolation. All four 150-h readings are present.',
      'S14: Metrics derive from final FASTAs, not mixed historical QUAST snapshots; particularly lur differs from historical PB03. No independent completeness/contamination or circularity claim is made.',
      'S1/S2: Genus subtitles reproduce culture-record labels, while main Table1 GTDB labels include revised names for buol, pmor, pscf and rfa. Stable model codes provide the crosswalk; color grouping follows the project common.py convention.',
      'Source sample IDs refer to well labels as recorded. Where physical plate coordinates were not supplied (S1/S7), they were not inferred. Replicates are described as wells, without assuming biological independence.'
    ])
    print(json.dumps({'figures':[f['number'] for f in figures],'pages':sum(len(f['pages']) for f in figures),'status':'complete'},indent=2))
