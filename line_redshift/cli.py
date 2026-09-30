"""Command-line reporting, plotting and optional local desktop notification."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from dataclasses import asdict
from .core import read_lines, read_pairs, match_pair, anchor_hypotheses, predicted_lines, Hypothesis, frame_status
from . import __version__


def desktop_notify(message):
    """Local OS notification only; no shell, email, webhooks or remote messages."""
    if sys.platform == 'darwin' and shutil.which('osascript'):
        script = 'on run argv\ndisplay notification (item 1 of argv) with title (item 2 of argv)\nend run'
        command = ['osascript', '-e', script, message, 'Line-redshift']
    elif shutil.which('notify-send'):
        command = ['notify-send', 'Line-redshift', message]
    else:
        print('Desktop notification unavailable; see terminal/report.', file=sys.stderr)
        return
    try:
        subprocess.run(command, check=True, capture_output=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        print('Desktop notification failed; scientific outputs were saved.', file=sys.stderr)


def load_spectra(manifest):
    import numpy as np
    spectra = {}
    if manifest is None:
        return spectra
    with open(manifest, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not {'pair_id', 'path', 'flux_unit'}.issubset(reader.fieldnames or []):
            raise ValueError('spectra manifest requires pair_id,path,flux_unit')
        for row in reader:
            path = Path(row['path'])
            if not path.is_absolute():
                path = Path(manifest).resolve().parent / path
            a = np.loadtxt(path, comments='#', ndmin=2)
            if a.shape[1] != 2 or a.shape[0] < 2 or not np.isfinite(a).all() or (a[:,0] <= 0).any():
                raise ValueError(f'{path.name}: require two finite numeric columns, frequency in GHz and flux')
            order = np.argsort(a[:,0])
            if np.any(np.diff(a[order,0]) <= 0):
                raise ValueError(f'{path.name}: spectrum frequencies must be distinct')
            spectra.setdefault(row['pair_id'], []).append((a[order,0], a[order,1], row['flux_unit'], path.name))
    for key, segments in spectra.items():
        if len({r[2] for r in segments}) != 1:
            raise ValueError(f'{key}: convert spectra to one common flux unit before plotting')
    return spectra


def plot_ladder(path, title, observed, hypotheses, lines, coverage, spectra, max_rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    shown = hypotheses[:max_rows]
    fig, (top, ladder) = plt.subplots(2, 1, figsize=(12, 4.0 + .85*len(shown)),
                                     gridspec_kw={'height_ratios':[2, max(1, len(shown)*.85)]}, sharex=True)
    frequencies = list(observed)
    for x,y,unit,label in spectra:
        top.step(x,y,where='mid',lw=.8,label=label)
        frequencies.extend([float(x.min()),float(x.max())])
    for lo,hi in coverage:
        frequencies.extend([lo,hi])
        for ax in [top,ladder]:ax.axvspan(lo,hi,color='lightgray',alpha=.3)
    lo,hi=min(frequencies),max(frequencies)
    margin=max((hi-lo)*.05, .01*(lo+hi)/2)
    lo,hi=lo-margin,hi+margin
    if spectra:
        top.set_ylabel(spectra[0][2]);top.legend(fontsize=7,loc='upper right')
    else:
        top.set_ylim(0,1);top.set_yticks([]);top.text(.5,.5,'Input line frequencies (no spectrum supplied)',ha='center',transform=top.transAxes)
    for i,f in enumerate(observed):
        for ax in [top,ladder]:ax.axvline(f,color='tab:red',ls='--',lw=.8)
        top.text(f,.95,f'{f:.5f} GHz',rotation=90,va='top',ha='right',transform=top.get_xaxis_transform(),fontsize=8)
    for i,h in enumerate(shown):
        selected={h.get('line_a'),h.get('line_b'),h.get('anchor_line')}
        visible=[]
        for line in lines:
            f=line.rest_ghz/(1+h['z'])
            if lo<=f<=hi:
                color='tab:blue' if line.name in selected else 'gray'
                ladder.plot([f,f],[i-.10,i+.10],color=color,lw=1.4)
                visible.append((f,line.name))
        clusters=[]
        for item in visible:
            if clusters and item[0]-clusters[-1][-1][0] < .015*(hi-lo):clusters[-1].append(item)
            else:clusters.append([item])
        for cluster in clusters:
            f=sum(v[0] for v in cluster)/len(cluster)
            color='tab:blue' if any(v[1] in selected for v in cluster) else 'gray'
            ladder.text(f,i+.13,'\n'.join(v[1] for v in cluster),rotation=90,fontsize=6,color=color,ha='center',va='bottom')
    labels=[]
    for h in shown:
        suffix=(f"  dv={h['velocity_mismatch_kms']:.1f} km/s" if 'velocity_mismatch_kms' in h else f"  {h['anchor_line']}")
        labels.append(f"z={h['z']:.5f}"+suffix)
    ladder.set_yticks(range(len(shown)));ladder.set_yticklabels(labels,fontsize=8)
    ladder.set_ylim(-.35,max(.9,len(shown)+.05));ladder.set_xlim(lo,hi);ladder.set_xlabel('Observed frequency (GHz)')
    if not shown:ladder.text(.5,.5,'No compatible assignment within the specified catalogue, z range and tolerance',ha='center',transform=ladder.transAxes,fontsize=9)
    top.set_title(title+'\nFrequency hypotheses only; no automatic redshift confirmation',fontsize=11)
    ladder.grid(axis='x',alpha=.2)
    fig.subplots_adjust(left=.24,right=.98,bottom=.09,top=.89,hspace=.16)
    fig.savefig(path,dpi=160,bbox_inches='tight');fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight');plt.close(fig)


def plot_spectrum_zoom(path, ident, observed, hypotheses, spectra):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(spectra), 1, figsize=(10, 3.2*len(spectra)), squeeze=False)
    for ax, (x,y,unit,label) in zip(axes[:,0],spectra):
        ax.step(x,y,where='mid',color='black',lw=.9)
        ax.axhline(0,color='gray',lw=.6)
        for f in observed:
            if min(x)<=f<=max(x):ax.axvline(f,color='tab:red',ls='--',label='Input frequency')
        for h,color in zip(hypotheses[:3],['tab:blue','tab:orange','tab:green']):
            for key,linekey in [('predicted_a_ghz','line_a'),('predicted_b_ghz','line_b')]:
                f=h.get(key)
                if f is not None and min(x)<=f<=max(x):
                    ax.axvline(f,color=color,lw=.9,label=f"{h[linekey]} z={h['z']:.5f}")
        ax.set(title=label,xlabel='Observed frequency (GHz)',ylabel=unit,xlim=(min(x),max(x)))
        handles,labels=ax.get_legend_handles_labels()
        if handles:ax.legend(fontsize=7,loc='best')
    fig.suptitle(ident+' | individual spectrum segments; first three assignments only',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(path,dpi=150);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)


def parser():
    p=argparse.ArgumentParser(description='Legacy lines.freq redshift ladders plus pair matching and notifications.')
    p.add_argument('--lines',type=Path,required=True,help='Whitespace frequency/name catalogue (legacy default MHz)')
    p.add_argument('--line-unit',choices=['Hz','kHz','MHz','GHz'],default='MHz')
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--pairs',type=Path,help='CSV: pair_id,frequency_a_ghz,frequency_b_ghz,frame_a,frame_b')
    mode.add_argument('--anchor-ghz',type=float,help='Original single-anchor redshift-ladder mode')
    p.add_argument('--tolerance-kms',type=float,default=100.,help='Pair-wise centroid velocity mismatch; not a sigma or per-line tolerance')
    p.add_argument('--z-min',type=float,default=0.)
    p.add_argument('--z-max',type=float,default=10.)
    p.add_argument('--allow-mixed-frames',action='store_true',help='Unverified exploratory matches only; NO frame conversion is performed')
    p.add_argument('--coverage',action='append',default=[],metavar='LOW:HIGH',help='Actual coverage interval in GHz; repeat for gaps')
    p.add_argument('--spectra-manifest',type=Path,help='CSV pair_id,path,flux_unit; 2-column GHz/flux ASCII segments')
    p.add_argument('--out',type=Path,default=Path('results'))
    p.add_argument('--plot-limit',type=int,default=12,help='Maximum rows per plot; all matches remain in CSV/JSON')
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--notify',action='store_true',help='Optional LOCAL desktop completion notification')
    p.add_argument('--version',action='version',version=__version__)
    return p


def run(args):
    if not (math.isfinite(args.z_min) and math.isfinite(args.z_max) and 0<=args.z_min<=args.z_max):
        raise ValueError('require finite 0 <= z_min <= z_max')
    if not math.isfinite(args.tolerance_kms) or args.tolerance_kms<=0 or args.plot_limit<1:
        raise ValueError('tolerance and plot-limit must be positive')
    coverage=[]
    for text in args.coverage:
        values=text.split(':')
        if len(values)!=2:raise ValueError('coverage format is LOW:HIGH in GHz')
        lo,hi=map(float,values)
        if not (math.isfinite(lo) and math.isfinite(hi) and 0<lo<hi):raise ValueError('require finite 0 < coverage LOW < HIGH')
        coverage.append((lo,hi))
    lines=read_lines(args.lines,args.line_unit)
    spectra=load_spectra(args.spectra_manifest)
    rows=[];groups=[];input_records=[];frames={}
    if args.pairs:
        pairs=read_pairs(args.pairs)
        unknown_ids=set(spectra)-{p.pair_id for p in pairs}
        if unknown_ids:raise ValueError('spectra manifest has unknown pair IDs: '+', '.join(sorted(unknown_ids)))
        # Validate all pairs before writing a partial output.
        for pair in pairs:
            input_records.append(asdict(pair));frames[pair.pair_id]=frame_status(pair,args.allow_mixed_frames)
            if frames[pair.pair_id].startswith('unverified'):print(f'WARNING {pair.pair_id}: {frames[pair.pair_id]}',file=sys.stderr)
            hits=[asdict(h) for h in match_pair(pair,lines,args.tolerance_kms,args.z_min,args.z_max,args.allow_mixed_frames)]
            groups.append((pair.pair_id,[pair.frequency_a_ghz,pair.frequency_b_ghz],hits))
            rows.extend(hits)
    else:
        if set(spectra)-{'anchor'}:raise ValueError('anchor-mode spectra must have pair_id=anchor')
        hits=anchor_hypotheses(args.anchor_ghz,lines,args.z_min,args.z_max)
        groups=[('anchor',[args.anchor_ghz],hits)]
        rows=[dict(pair_id='anchor',**h,status='single_line_ambiguous_hypothesis') for h in hits]
    args.out.mkdir(parents=True,exist_ok=True)
    summary=[];predictions=[];reports=['# Frequency-consistent redshift hypotheses','',
        'Every row is a hypothesis, not a secure redshift or a detection probability. Sorting by velocity mismatch is not probabilistic ranking.', '',
        f'Pair-wise tolerance: {args.tolerance_kms:g} km/s. Redshift interval: [{args.z_min:g}, {args.z_max:g}].',
        'No line-strength prior, missing-line exclusion, frame conversion, or centroid fitting is performed. The default tolerance is illustrative; calibrate it for your measurements.', '']
    for index,(ident,observed,hits) in enumerate(groups,1):
        basename=f'pair_{index:03d}' if args.pairs else 'anchor'
        summary.append(dict(pair_id=ident,frame_status=frames.get(ident,'single-anchor mode'),hypotheses=len(hits),figure=basename+'.png' if not args.no_plots else None))
        print(f'{ident}: {len(hits)} hypothesis(es)')
        for h in hits[:3]:
            print('  z={:.6f}  {}'.format(h['z'],h.get('line_a',h.get('anchor_line',''))+' / '+h.get('line_b','anchor')))
        reports.extend([f'## {ident}','',f'{len(hits)} compatible assignment(s).',''])
        if hits:
            reports.extend(['| Assignment | z | mismatch (km/s) | frame status |','|---|---:|---:|---|'])
            for h in hits:
                assignment=h.get('line_a',h.get('anchor_line',''))+' / '+h.get('line_b','anchor')
                reports.append(f"| {assignment} | {h['z']:.6f} | {h.get('velocity_mismatch_kms',0):.3f} | {h.get('frame_status','single-anchor mode')} |")
                predictions.append(dict(pair_id=ident,z=h['z'],assignment=assignment,lines=predicted_lines(h['z'],lines,coverage)))
            reports.append('')
        else:reports.extend(['No compatible assignment under these settings. This does not rule out a real source or transitions absent from the catalogue.',''])
        if not args.no_plots:
            plot_ladder(args.out/(basename+'.png'),ident,observed,hits,lines,coverage,spectra.get(ident,[]),args.plot_limit)
            reports.extend([f'![{ident}]({basename}.png)',''])
            if spectra.get(ident):
                plot_spectrum_zoom(args.out/(basename+'_spectra.png'),ident,observed,hits,spectra[ident])
                reports.extend([f'![{ident} spectrum segments]({basename}_spectra.png)',''])
    fields=list(dict.fromkeys(k for r in rows for k in r)) or list(Hypothesis.__dataclass_fields__)
    with open(args.out/'hypotheses.csv','w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    payload=dict(version=__version__,inputs=input_records,anchor_ghz=args.anchor_ghz,parameters=dict(tolerance_kms=args.tolerance_kms,z_min=args.z_min,z_max=args.z_max,line_unit=args.line_unit,allow_mixed_frames=args.allow_mixed_frames,coverage_ghz=coverage,plot_limit=args.plot_limit),line_catalogue_sha256=hashlib.sha256(args.lines.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.pairs.read_bytes()).hexdigest() if args.pairs else None,summary=summary,hypotheses=rows,predictions=predictions)
    (args.out/'hypotheses.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (args.out/'report.md').write_text('\n'.join(reports)+'\n',encoding='utf-8')
    message=f'{sum(bool(h) for _,_,h in groups)}/{len(groups)} inputs have frequency-consistent hypotheses ({len(rows)} assignments).'
    print(message);print('Saved:',args.out)
    if args.notify:desktop_notify(message)
    return 0


def main(argv=None):
    p=parser();args=p.parse_args(argv)
    try:return run(args)
    except (ValueError,OSError) as exc:
        p.error(str(exc))
