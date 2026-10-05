"""Shared figure conventions for the current Supplementary V5 assets."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
EXTERNAL_DATA_ROOT = ROOT / 'external_data'
OUT = ROOT / 'outputs' / 'supplementary'
FIG = OUT / 'figures'
DATA = OUT / 'source_data'
MANIFEST = OUT / 'manifests'
HIST = EXTERNAL_DATA_ROOT
sys.path.insert(0, str(ROOT))
from publication_style import apply_publication_style, setup_axes, add_panel_letter, PUBLICATION_COLORS

ORDER = ['bbet','cauj','dri','mlot','pzu','tro','buol','vpe','lgu','lur','pmor','pscf','racd','cuvk','cxy','mip','rfa','cpi','flk','mgk','cohm']
GROUPS = {**dict.fromkeys(['bbet','cauj','dri','mlot','pzu','tro'],'Alphaproteobacteria'),
          **dict.fromkeys(['buol','vpe'],'Betaproteobacteria'),
          **dict.fromkeys(['lgu','lur','pmor','pscf','racd'],'Gammaproteobacteria'),
          **dict.fromkeys(['cuvk','cxy','mip','rfa'],'Actinobacteria'),
          **dict.fromkeys(['cpi','flk','mgk'],'Bacteroidetes'),'cohm':'Firmicutes'}
GROUP_COLORS = {'Actinobacteria':'#1b9e77','Alphaproteobacteria':'#d95f02',
                'Betaproteobacteria':'#7570b3','Gammaproteobacteria':'#e7298a',
                'Firmicutes':'#66a61e','Bacteroidetes':'#e6ab02'}
COLORS = {s:GROUP_COLORS[GROUPS[s]] for s in ORDER}

V5_STEM_MAP = {
    'S03a_genome_and_curation': 'S06a_genome_and_curation',
    'S03b_substrate_prediction_errors': 'S06b_substrate_prediction_errors',
    'S04_quality_and_thermodynamics': 'S08_quality_and_thermodynamics',
    'S05_gapfill_evidence_scores': 'S07_gapfill_evidence_scores',
    'S06_nutrient_alternative_cardinality': 'S09_nutrient_alternative_cardinality',
    'S07_supplementation_validation': 'S10_supplementation_validation',
    'S11_optimal_vs_pooled': 'S14_optimal_vs_pooled',
    'S11_niche_overlap_sensitivity': 'S14_niche_overlap_sensitivity',
    'S12_ptyg_biosynthetic_similarity': 'S18_ptyg_biosynthetic_similarity',
    'S13_gene_sets_and_medium_overlap': 'S20_gene_sets_and_medium_overlap',
    'S13_bbb_memberships': 'S20_bbb_memberships',
    'S14_genome_assembly_structure': 'S03_genome_assembly_structure',
    'S15_bbb_route_sizes': 'S19_bbb_route_sizes',
    'S15_bbb_route_alternatives': 'S19_bbb_route_alternatives',
}

V5_NUMBER_MAP = {3: 6, 4: 8, 5: 7, 6: 9, 7: 10, 8: 11,
                 11: 14, 12: 18, 13: 20, 14: 3, 15: 19}

def v5_stem(stem):
    if stem.startswith('S08_class_coverage_'):
        return stem.replace('S08_', 'S11_', 1)
    if stem.startswith(('S09_sharing_', 'S09_other_', 'S10_exchange_classes_')):
        return 'historical_' + stem
    return V5_STEM_MAP.get(stem, stem)

def style():
    apply_publication_style(font_family='Arial', extra_rc={'font.size':8,'axes.labelsize':8,
       'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
       'legend.fontsize':7,'savefig.bbox':None,'svg.fonttype':'none','pdf.fonttype':42})

def save(fig, stem):
    """Write V5-named, live-text SVGs with editable vector marks."""
    stem = v5_stem(stem)
    FIG.mkdir(parents=True, exist_ok=True)
    # imshow normally embeds a bitmap even in an SVG. Replace each numerical
    # heatmap with a mesh so Illustrator can edit the cells as well as labels.
    for ax in fig.axes:
        for im in list(ax.images):
            values=im.get_array()
            if values.ndim != 2:
                raise ValueError('Unexpected image rather than a numerical heatmap')
            left,right,bottom,top=im.get_extent()
            rows,cols=values.shape
            xx=np.linspace(left,right,cols+1)
            yy=np.linspace(top,bottom,rows+1) if im.origin=='upper' else np.linspace(bottom,top,rows+1)
            xlim,ylim=ax.get_xlim(),ax.get_ylim()
            ax.pcolormesh(xx,yy,values,cmap=im.get_cmap(),norm=im.norm,
                          shading='flat',edgecolors='none',antialiased=False,
                          alpha=im.get_alpha(),zorder=im.get_zorder(),rasterized=False)
            im.remove()
            ax.set_xlim(xlim);ax.set_ylim(ylim)
        for artist in ax.collections:
            artist.set_rasterized(False)
    for ext in ['pdf','svg','png']:
        fig.savefig(str(FIG / (stem+'.'+ext)), dpi=400, bbox_inches=None, facecolor='white')
    plt.close(fig)
    return stem

def manifest(name, figures, audit=None):
    MANIFEST.mkdir(parents=True, exist_ok=True)
    for figure in figures:
        old_number = figure.get('number')
        if old_number in V5_NUMBER_MAP:
            figure['source_number'] = old_number
            figure['number'] = V5_NUMBER_MAP[old_number]
        elif old_number in (9, 10):
            figure['status'] = 'historical panel not retained in V5'
    for figure in figures:
        old_number = figure.get('number')
        if old_number in V5_NUMBER_MAP:
            figure['source_number'] = old_number
            figure['number'] = V5_NUMBER_MAP[old_number]
        elif old_number in (9, 10):
            figure['status'] = 'historical panel not retained in V5'
    caption_file=MANIFEST/'captions_v4.json'
    if caption_file.exists():
        captions=json.loads(caption_file.read_text())
        for figure in figures:
            figure.update(captions[str(figure['number'])])
    (MANIFEST/(name+'.json')).write_text(json.dumps({'figures':figures,'audit':audit or []},indent=2))

def source_table(frame, filename):
    DATA.mkdir(parents=True,exist_ok=True)
    frame.to_csv(DATA/filename,index=False,sep='\t')

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

style()
