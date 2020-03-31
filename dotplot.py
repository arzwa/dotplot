#!/usr/env/python3
# Author: Arthur Zwaenepoel (2020)
import pandas as pd
import matplotlib.pyplot as plt
import logging
from matplotlib import transforms

__author__ = "Arthur Zwaenepoel <arzwa@psb.vib-ugent.be>"
__version__ = "0.1.0"

def parse_sgff(fname):
    gff = pd.read_csv(fname, header=None, sep="\t")
    return process_sgff(gff)

def process_sgff(gff):
    gff.columns  = ["chr", "gene", "start", "stop"]
    gff.index    = gff['gene']
    gff['stop']  = gff['stop'].apply(lambda x: str(x).split(',')[0])
    gff['start'] = gff['start'].apply(lambda x: str(x).split(',')[0])
    gff['stop']  = pd.to_numeric(gff['stop'])
    gff['start'] = pd.to_numeric(gff['start'])

    # sort genomic elements by length and compute coordinates
    elems = gff.groupby(['chr']).apply(
        lambda x: max(x['stop']) - min(x['start'])).sort_values(
            ascending=False)
    curr_offset = 0
    gff['coord1'] = -1
    gff['coord2'] = -1
    gff['elemlen'] = -1.
    gff['elembound'] = -1
    for elem in elems.index:
        rows = gff[gff['chr'] == elem].index
        gff.loc[rows,'elemlen'] = elems.loc[elem]
        gff.loc[rows,'coord1'] = gff.loc[rows,'start'] + curr_offset
        gff.loc[rows,'coord2'] = gff.loc[rows,'stop' ] + curr_offset
        gff.loc[rows,'elembound'] = curr_offset
        curr_offset += elems.loc[elem]
    return gff

def parse_gff(fname, feature="gene", attribute="ID"):
    gff = pd.read_csv(fname, header=None, sep="\t")
    gff = gff.loc[gff[2] == feature]
    gff[9] = gff[8].apply(lambda x: extract_from_gffinfo(x, attribute))
    sgff = gff.dropna()[[0,9,3,4]]
    return process_sgff(sgff)

def extract_from_gffinfo(x, key):
    d = {}
    for pair in x.split(";"):
        k,v=pair.split("=")
        d[k] = v
        if key in d.keys():
            return d[key]
        else:
            return np.nan

def parse_mcscan(fname):
    with open(fname, 'r') as f:
        content = f.read().split("## Alignment")[1:]
    anchors = []
    for aln in content:
        lines = aln.split("\n")
        for line in lines[1:-1]:
            l = line.split()
            anchors.append("__".join(sorted([l[2], l[3]])))
    return set(anchors)

def parse_pairs(fname):
    ds = []
    with open(fname, 'r') as f:
        for line in f:
            l = line.split()
            g1, g2 = l[0], l[1]
            if g1 == g2:
                continue
            pair = "__".join(sorted([g1,g2]))
            d = {"pair":pair, "g1": g1, "g2": g2}
            ds.append(d)
    df = pd.DataFrame.from_dict(ds)
    df.index = df['pair']
    return df

def parse_mcl(fname):
    ds = []
    with open(fname, 'r') as f:
        for line in f:
            l = line.split()
            for i in range(len(l)):
                for j in range(i+1,len(l)):
                    g1, g2 = l[i], l[j]
                    pair = "__".join(sorted([g1,g2]))
                    d = {"pair":pair, "g1": g1, "g2": g2}
                    ds.append(d)
    df = pd.DataFrame.from_dict(ds)
    df.index = df['pair']
    return df

def assemble_df(coords, hits, anchors=None, coord='coord1', minlen=0):
    coords_ = coords[coords['elemlen'] > minlen]
    hits['x'] = list(coords_.reindex(hits['g1'])['coord1'])
    hits['y'] = list(coords_.reindex(hits['g2'])['coord1'])
    hits['anchor'] = False
    if anchors:
        hits.loc[anchors, 'anchor'] = True
    l = len(hits)
    df = hits.dropna().drop('pair', axis=1)
    m = len(df)
    notfound = hits.loc[list(set(hits.index) - set(df.index))]
    logging.warning("Went from {} to {} after dropna()".format(l, m))
    df['a'] = df['x'].where(df.x < df.y, df['y'])
    df['b'] = df['y'].where(df.x < df.y, df['x'])
    df['color'] = 'k'
    df.loc[df['anchor'] == True, 'color'] = 'red'
    return df, coords_, notfound

def genomescatter_triangle(ax, df, coords,
        minlen=10**6, lw=0.1, s=0.8, alpha=0.2):
    base = plt.gca().transData
    rot = transforms.Affine2D().rotate_deg(-45)
    ax.scatter(df['a'], df['b'],
        alpha=alpha, s=s, c=df['color'], transform=rot+base);
    anchors = df.loc[df['anchor'] == True]
    ax.scatter(anchors['a'], anchors['b'],
        alpha=1, s=s, c=anchors['color'], transform=rot+base);
    bounds = list(coords['elembound'].unique()) + [max(coords['coord2'])]
    ymin, ymax = min(bounds), max(bounds)
    ax.vlines(bounds, ymin=ymin, ymax=ymax, lw=lw, transform=rot+base)
    ax.hlines(bounds, xmin=ymin, xmax=ymax, lw=lw, transform=rot+base)
    ax.set_ylim(0, (ymax**2/2)**0.5)
    ax.set_xlim(0, 2*(ymax**2/2)**0.5)
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.hlines([0], xmin=0, xmax=2*(ymax**2/2)**0.5, lw=lw)
    return ax
