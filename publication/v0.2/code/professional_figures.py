"""Publication artwork from unchanged CSV evidence; no generative image model."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parent
BLUE, ORANGE, BLACK = '#0072B2', '#D55E00', '#202020'
COLORS = [BLUE, ORANGE]
NAMES = ['Interface - core', 'Structured extra - core']
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,
    'axes.titlesize':10,'axes.titleweight':'bold','axes.linewidth':.85,
    'xtick.labelsize':8,'ytick.labelsize':9,'lines.linewidth':1.2,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'svg.hashsalt':'membranecal-academic-v2','text.color':BLACK,
    'axes.labelcolor':BLACK,'axes.edgecolor':BLACK})


def rows(name):
    with (ROOT/'tables'/name).open(encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f))


def save(fig, stem):
    for ext in ['pdf','svg','png','tiff']:
        opts = {'metadata':{'CreationDate':None,'ModDate':None}} if ext=='pdf' else {}
        if ext=='svg': opts={'metadata':{'Date':None}}
        if ext=='tiff': opts={'pil_kwargs':{'compression':'tiff_lzw'}}
        fig.savefig(ROOT/f'figures/{stem}.{ext}',dpi=600,facecolor='white',**opts)
    from PIL import Image
    for ext in ['png','tiff']:
        path = ROOT/f'figures/{stem}.{ext}'
        with Image.open(path) as image:
            rgb = image.convert('RGB')
        rgb.save(path, dpi=(600,600), **({'compression':'tiff_lzw'} if ext=='tiff' else {}))
    plt.close(fig)


def clean(ax):
    for side in ['top','right','left']: ax.spines[side].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=8)
    ax.tick_params(axis='x',length=3,width=.85)


def flow():
    fig,ax=plt.subplots(figsize=(7.2,7.0))
    fig.subplots_adjust(left=.015,right=.985,top=.99,bottom=.015)
    ax.set(xlim=(0,1),ylim=(0,1)); ax.axis('off')
    def node(x,y,w,h,title,detail=''):
        ax.add_patch(Rectangle((x,y),w,h,fill=False,lw=1.05,edgecolor=BLACK))
        ax.text(x+w/2,y+h*.65 if detail else y+h/2,title,ha='center',va='center',weight='bold',fontsize=9)
        if detail: ax.text(x+w/2,y+h*.27,detail,ha='center',va='center',fontsize=8.2,linespacing=1.35)
    def arrow(x1,y1,x2,y2):
        ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops={'arrowstyle':'-|>','lw':1.05,'color':BLACK,'mutation_scale':9})
    node(.04,.90,.48,.075,'8,915 OPM index records')
    arrow(.28,.90,.28,.66)
    ax.plot([.28,.60],[.79,.79],color=BLACK,lw=1.05)
    ax.text(.62,.85,'Sequential exclusions',weight='bold',fontsize=9)
    ax.text(.62,.815,'65  Development-entry overlap\n3,540  Type / resolution availability\n1,497  Resolution / source boundary',va='top',fontsize=8.2,linespacing=1.6)
    node(.04,.57,.48,.09,'800 candidate entries','Unique-code and per-family selection; cap 800')
    arrow(.28,.57,.28,.44)
    ax.plot([.28,.60],[.505,.505],color=BLACK,lw=1.05)
    ax.text(.62,.54,'Preparation exclusions',weight='bold',fontsize=9)
    ax.text(.62,.505,'286  Ineligible entries\n5  Source / parser / QC failures',va='top',fontsize=8.2,linespacing=1.6)
    node(.04,.35,.48,.09,'509 eligible entries','622 eligible protein cases')
    arrow(.28,.35,.28,.22)
    ax.plot([.28,.60],[.285,.285],color=BLACK,lw=1.05)
    ax.text(.62,.305,'Canonical representative selection',weight='bold',fontsize=8.5)
    ax.text(.62,.275,'68 duplicate protein cases removed',fontsize=8.2)
    node(.04,.125,.48,.095,'554 canonical proteins','473 entries; 223,714 paired-observed residues')
    arrow(.52,.173,.60,.173)
    ax.text(.62,.225,'Primary matching',weight='bold',fontsize=9.5)
    ax.text(.62,.195,'Interface - core: 437 proteins\n45.27% effective common mass',va='top',fontsize=8.2,linespacing=1.35)
    ax.text(.62,.135,'Structured extra - core: 229 proteins\n23.81% effective common mass',va='top',fontsize=8.2,linespacing=1.35)
    ax.text(.04,.045,'Counts use the units named at each stage. Candidate selection is capped, not exhaustive.',fontsize=8)
    save(fig,'01_cohort_flow')


def primary():
    data=rows('primary_effects.csv')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.5),gridspec_kw={'width_ratios':[1,1.45]})
    fig.subplots_adjust(left=.22,right=.98,bottom=.22,top=.79,wspace=.26)
    for k,ax in enumerate(axes):
        clean(ax); ax.axvline(0,color='#707070',lw=.85)
        for i,r in enumerate(data):
            x,lo,hi=[float(r[v]) for v in ['estimate_points','lower_points','upper_points']]
            ax.errorbar(x,1-i,xerr=[[x-lo],[hi-x]],fmt='o',color=COLORS[i],capsize=3,ms=5)
            if k: ax.text(x,1-i+.18,f'{x:+.3f} [{lo:+.3f}, {hi:+.3f}]',ha='center',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','pad':1})
        ax.set_ylim(-.4,1.55); ax.set_yticks([1,0],NAMES if k==0 else ['', ''])
        ax.set_xlabel('Agreement - confidence contrast\n(percentage points)')
    axes[0].set_xlim(-6,6); axes[0].set_xticks([-5,0,5])
    for x in [-5,5]: axes[0].axvline(x,color='#555555',ls=(0,(3,3)),lw=1.0)
    axes[1].set_xlim(-1.15,1.4)
    axes[0].set_title('a   Investigator-chosen\n     ±5-point decision scale',loc='left',pad=10)
    axes[1].set_title('b   Detail of the same estimates',loc='left',pad=16)
    fig.text(.22,.965,'97.5% marginal block-bootstrap intervals',fontsize=9,va='top')
    save(fig,'02_primary_effects')


def support():
    data=rows('support_balance.csv')
    fig,axes=plt.subplots(2,1,figsize=(7.2,5.0))
    fig.subplots_adjust(left=.25,right=.92,top=.85,bottom=.12,hspace=1.08)
    for ax in axes: clean(ax)
    a,b=axes
    for i,r in enumerate(data):
        y=1-i
        for field,off,marker in [('raw_retention',.12,'o'),('effective_retention',-.12,'D')]:
            x=float(r[field])*100
            a.plot(x,y+off,marker=marker,color=COLORS[i],ms=6 if marker!='o' else 5,mfc='white' if marker!='o' else COLORS[i],mew=1.5 if marker!='o' else 1.0)
            a.text(x+2,y+off,f'{x:.2f}%',va='center',fontsize=8)
        for field,off,marker in [('group_balanced_confidence_difference',.12,'o'),('maximum_absolute_protein_confidence_difference',-.12,'^')]:
            x=float(r[field])*100
            b.plot(x,y+off,marker=marker,color=COLORS[i],ms=6 if marker!='o' else 5,mfc='white' if marker!='o' else COLORS[i],mew=1.5 if marker!='o' else 1.0)
            b.text(x-.025 if x<0 else x+.025,y+off,f'{x:+.3f}' if marker=='o' else f'{x:.3f}',ha='right' if x<0 else 'left',va='center',fontsize=8)
    a.set(xlim=(0,100),ylim=(-.45,1.45),xlabel='Eligible regional residues retained (%)')
    a.set_yticks([1,0],[f'{NAMES[i]}\nN = {int(r["available"]):,}' for i,r in enumerate(data)])
    a.set_title('a   Common support',loc='left',pad=26)
    a.plot([],[],'o',color=BLACK,label='Residues in retained cells',ms=4)
    a.plot([],[],'D',mfc='white',color=BLACK,label='Effective common mass',ms=5,mew=1.5)
    a.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.02,1.01),ncol=2,fontsize=8,handletextpad=.4,columnspacing=1)
    b.axvline(0,color='#707070',lw=.85)
    b.set(xlim=(-.2,.85),ylim=(-.45,1.45),xlabel='Residual confidence difference (pLDDT points)')
    b.set_yticks([1,0],NAMES)
    b.set_title('b   Residual confidence balance',loc='left',pad=26)
    b.plot([],[],'o',color=BLACK,label='Signed group-balanced mean',ms=4)
    b.plot([],[],'^',mfc='white',color=BLACK,label='Maximum absolute protein value',ms=5,mew=1.5)
    b.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.02,1.01),ncol=2,fontsize=7.5,handletextpad=.4,columnspacing=.8)
    save(fig,'03_support_balance')


def sensitivities():
    data=rows('sensitivity_effects.csv')
    contrasts=list(dict.fromkeys(r['contrast'] for r in data))
    fig,axes=plt.subplots(1,2,figsize=(7.2,5.7))
    fig.subplots_adjust(left=.27,right=.965,bottom=.19,top=.89,wspace=.27)
    for j,(ax,contrast) in enumerate(zip(axes,contrasts)):
        selected=[r for r in data if r['contrast']==contrast]
        assert len(selected)==12
        clean(ax); ax.axvline(0,color='#707070',lw=.85)
        for i,r in enumerate(selected):
            x,lo,hi=[float(r[v]) for v in ['estimate_points','lower_points','upper_points']]
            ax.errorbar(x,11-i,xerr=[[x-lo],[hi-x]],fmt='D' if i==0 else 'o',color=COLORS[j],capsize=2,ms=5.5 if i==0 else 4)
            ax.text(2.65,11-i,r['proteins'],ha='right',va='center',fontsize=8)
        ax.set(xlim=(-1.3,2.85),ylim=(-.6,12),xlabel='Contrast (percentage points)')
        ax.set_xticks([-1,0,1,2]); ax.spines['bottom'].set_bounds(-1,2)
        ax.set_yticks(range(11,-1,-1),[r['label'] for r in selected] if j==0 else ['']*12)
        if j==0: ax.get_yticklabels()[0].set_fontweight('bold')
        ax.set_title(('a   ' if j==0 else 'b   ')+NAMES[j],loc='left',pad=16)
        ax.text(2.65,11.7,'n',ha='right',fontsize=8,style='italic')
    fig.text(.27,.035,'97.5% marginal intervals; n = matched proteins.\nNo simultaneous guarantee across sensitivity scenarios.',fontsize=8,linespacing=1.4)
    save(fig,'04_sensitivities')


def cases():
    data=rows('descriptive_cases.csv')
    fig,ax=plt.subplots(figsize=(7.2,4.4)); fig.subplots_adjust(left=.24,right=.96,bottom=.20,top=.91)
    clean(ax); ax.axvline(0,color='#707070',lw=.85)
    for i,r in enumerate(data):
        x=float(r['mean_error_points']); color=BLUE if x>=0 else ORANGE
        ax.plot([0,x],[7-i,7-i],color=color,lw=1.2)
        ax.plot(x,7-i,'o',color=color,ms=5,mew=1)
        ax.annotate(f'{x:+.2f}',(x,7-i),xytext=(6 if x>=0 else -6,0),textcoords='offset points',ha='left' if x>=0 else 'right',va='center',fontsize=8)
    ax.set_yticks(range(7,-1,-1),[f'{r["accession"]} / {r["entry_id"]}' for r in data])
    ax.set(xlim=(-46,22),ylim=(-.6,7.6),xlabel='All-observed-residue mean (agreement - confidence)\n(percentage points)')
    fig.text(.24,.97,'Descriptive unmatched protein means',weight='bold',fontsize=10,va='top')
    save(fig,'05_descriptive_cases')


def main():
    flow(); primary(); support(); sensitivities(); cases()
    manifest=[{'path':p.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'origin':'Matplotlib vector artwork from unchanged retained evidence'} for p in sorted((ROOT/'figures').iterdir()) if p.suffix in ['.pdf','.svg','.png','.tiff']]
    (ROOT/'evidence/figure_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8',newline='\n')
    print('Five figures exported as PDF/SVG and 600 dpi PNG/TIFF')


if __name__=='__main__': main()
