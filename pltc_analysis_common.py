"""
pltc_analysis_common.py
Shared, deterministic image-analysis routines for PltC–CD31 colocalization in
CLARITY-cleared small intestine (Alexa Fluor 488-labelled PltC intoxication).

Channels (Imaris .ims / HDF5):  ch0 = DAPI (405), ch1 = CD31 (640/AF647), ch2 = PltC (488).

Requires: numpy, scipy, scikit-image, h5py  (Python 3.10).
The pipeline is fully deterministic: given the same input files it returns
bit-identical results (no random components).
"""
import numpy as np, h5py
from scipy import ndimage
from skimage.filters import frangi
from skimage.morphology import remove_small_objects, white_tophat, disk
from skimage.transform import rescale
from skimage.measure import label, regionprops

DOWNSCALE = 0.5     # in-plane downscale for speed (deterministic, anti-aliased)

def load_ims(path):
    """Return (dapi, cd31, pltc) single optical plane (float32) and pixel size (um)."""
    hf = h5py.File(path, 'r')
    b = hf['DataSet']['ResolutionLevel 0']['TimePoint 0']
    ch = [b['Channel %d' % i]['Data'][0].astype(np.float32) for i in range(3)]
    inf = hf['DataSetInfo']['Image'].attrs
    ext = lambda k: float(b''.join(inf[k]).decode())
    ps = (ext('ExtMax0') - ext('ExtMin0')) / int(b''.join(inf['X']).decode())
    hf.close()
    return ch[0], ch[1], ch[2], ps

def _ds(a):
    return rescale(a, DOWNSCALE, anti_aliasing=True, preserve_range=True).astype(np.float32)

def villus_interior(dapi, ps, rim_um=16):
    """Villus core = inside the DAPI+ epithelial ring, minus the ~rim_um epithelial rim."""
    d = ndimage.gaussian_filter(dapi, 2)
    solid = ndimage.binary_fill_holes(ndimage.binary_closing(d > np.percentile(d, 55), iterations=3))
    solid = remove_small_objects(solid, 4000)
    edt = ndimage.distance_transform_edt(solid)
    return edt > (rim_um / ps)

def whole_tissue(dapi, ps):
    """Whole DAPI+ tissue mask (no villus-interior restriction, no rim removal)."""
    d = ndimage.gaussian_filter(dapi, 2)
    solid = ndimage.binary_fill_holes(ndimage.binary_closing(d > np.percentile(d, 40), iterations=3))
    return remove_small_objects(solid, 2000)

def cd31_endothelium(cd31, ps):
    """Tubular CD31+ endothelial mask: vesselness filter + size/shape gating."""
    sm = ndimage.median_filter(cd31, 2)
    flat = sm - ndimage.gaussian_filter(sm, int(6 / ps)); flat[flat < 0] = 0
    v = frangi(flat / (flat.max() + 1e-6), sigmas=range(1, 6, 1), black_ridges=False)
    thr = v > np.percentile(v[v > 0], 80) if (v > 0).any() else v > 0
    thr = ndimage.binary_closing(remove_small_objects(thr, 40), iterations=1)
    lb = label(thr); keep = np.zeros_like(thr)
    for r in regionprops(lb):
        if r.area > 2500 and r.eccentricity < 0.8 and r.extent > 0.45:   # big round debris -> drop
            continue
        keep[lb == r.label] = True
    return flat * keep, keep

def metrics(dapi, cd31, pltc, ps, region='villus'):
    """
    Per-image PltC–CD31 colocalization metrics on the chosen analysis region.
      region='villus' -> villus interior (main analysis)
      region='whole'  -> whole DAPI+ tissue (robustness control)
    Returns dict, or None if the image fails QC (no vessel / too small).
    """
    core = villus_interior(dapi, ps) if region == 'villus' else whole_tissue(dapi, ps)
    if core.sum() < 800:
        return None
    cint, vmask = cd31_endothelium(cd31, ps)
    vcore = vmask & core
    if vcore.sum() < 20:
        return None
    near = ndimage.binary_dilation(vcore, iterations=max(1, int(round(1.2 / ps))))     # on-endothelium
    off = core & ~ndimage.binary_dilation(vcore, iterations=3)                          # off-vessel core
    if off.sum() < 50:
        return None
    r = max(1, int(round(2.2 / ps)))
    t488 = ndimage.median_filter(white_tophat(pltc, disk(r)), 1)
    t640 = ndimage.median_filter(white_tophat(cd31, disk(r)), 1)
    # remove broadband autofluorescence = objects bright in BOTH 488 and 640
    cob = (t488 > np.percentile(t488[core], 99.0)) & (t640 > np.percentile(t640[core], 98.0))
    tc = t488.copy(); tc[cob] = 0
    Timg = np.percentile(tc[off], 99.9)                 # per-image autofluorescence floor
    pltc_pos = (tc > Timg) & core
    sp = remove_small_objects((tc > Timg) & near, 3)
    n_punc = len([p for p in regionprops(label(sp)) if 2 <= p.area <= 400])
    ves_um2 = float(vcore.sum()) * ps * ps
    return dict(
        n_puncta      = n_punc,
        puncta_density= n_punc / ves_um2 * 1000.0,                         # per 1000 um^2 vessel
        enrichment    = tc[near].mean() / max(tc[off].mean(), 1e-6),       # on/off-vessel
        M1            = 100.0 * tc[core & near].sum() / max(tc[core].sum(), 1e-6),
        M2            = 100.0 * t640[core & pltc_pos].sum() / max(t640[core].sum(), 1e-6),
        coverage      = 100.0 * (pltc_pos & near).sum() / max(near.sum(), 1),
        vessel_area_um2 = ves_um2,
    )

def run_folder(base, region, out_csv):
    """Process the 2D Images folder (4 numbered condition subfolders) and write a per-image CSV."""
    import glob, os, re, csv
    COND = ['1. B6_NoToxin', '2. B6_Intoxication', '3. Tek-Cre_NoToxin', '4. Tek-Cre_Intoxication']
    rows = []
    for c in COND:
        for f in sorted(glob.glob(os.path.join(base, c, '*.ims'))):
            try:
                dapi, cd31, pltc, ps = load_ims(f)
            except Exception as e:
                print('skip (unreadable):', os.path.basename(f), e); continue
            res = metrics(_ds(dapi), _ds(cd31), _ds(pltc), ps * 2, region=region)
            if res:
                res['condition'] = c.split('. ')[1]; res['file'] = os.path.basename(f)
                rows.append(res)
    cols = ['condition', 'file', 'vessel_area_um2', 'n_puncta', 'puncta_density',
            'enrichment', 'M1', 'M2', 'coverage']
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({k: r[k] for k in cols})
    print('wrote %d rows -> %s' % (len(rows), out_csv))
    return rows
