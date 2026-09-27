"""Read-only native MCA semantic samples; no graph or tree representation."""
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.morphology import skeletonize
from .topbrain_dataset import discover, load_case
from .topbrain_qc import TopBrainError, crop, contact, physical_points, sha256, write_json

LABELS = {0:'UNKNOWN', 1:'M1', 2:'M2', 3:'M3'}
VERSION = 'native-semantic-points-1'


@dataclass
class Geometry:
    points: np.ndarray
    mass: np.ndarray
    root: np.ndarray
    direction: np.ndarray
    case_id: str
    side: str


@dataclass
class SemanticPoints:
    geometry: Geometry
    labels: np.ndarray
    native_points: np.ndarray
    native_labels: np.ndarray
    representative_ids: np.ndarray
    assignment: np.ndarray
    metadata: dict


def downsample(points, maximum=275):
    """Grid means select nearest actual sample; ties use original sample order.

    Geometry alone determines the cells and representatives, so hidden target
    labels cannot leak through sampling. A mixed cell takes its representative's
    actual native label, never an interpolated class. Every original point maps
    to exactly one representative and contributes its count to sample mass.
    """
    points=np.asarray(points,float)
    if len(points)<=maximum:return np.arange(len(points)),np.arange(len(points)),np.ones(len(points))/len(points)
    origin=points.min(axis=0);lo=0.;hi=np.linalg.norm(np.ptp(points,axis=0))
    for _ in range(40):
        step=(lo+hi)/2
        cells,inverse=np.unique(np.floor((points-origin)/step).astype(np.int64),axis=0,return_inverse=True)
        if len(cells)>maximum:lo=step
        else:hi=step
    cells,inverse=np.unique(np.floor((points-origin)/hi).astype(np.int64),axis=0,return_inverse=True)
    ids=[]
    for i in range(len(cells)):
        members=np.flatnonzero(inverse==i);center=points[members].mean(axis=0)
        ids.append(members[np.argmin(np.linalg.norm(points[members]-center,axis=1))])
    mass=np.bincount(inverse,minlength=len(ids)).astype(float);mass/=mass.sum()
    return np.asarray(ids),inverse,mass


def trunk_direction(points, root):
    distances=np.linalg.norm(points-root,axis=1)
    size=np.percentile(distances,90)
    near=(distances>=max(.5,.025*size))&(distances<=max(3.,.12*size))
    if near.sum()<3:near=np.argsort(distances)[1:min(12,len(points))]
    direction=(points[near]-root).mean(axis=0)
    if np.linalg.norm(direction)<1e-8:raise TopBrainError('ROOT_DIRECTION_AMBIGUOUS','No proximal outward direction')
    return direction/np.linalg.norm(direction)


def extract(case, side, maximum=275):
    values=[case.label_map.resolve(f'{side}-M{i}') for i in (1,2,3)]
    masks=[case.labels==v for v in values]
    if not all(m.any() for m in masks):raise TopBrainError('MISSING_NATIVE_SEGMENT',f'{case.paths.case_id}/{side}')
    union=np.logical_or.reduce(masks)
    interface=contact(union,case.mask(f'{side}-ICA'))
    patches,n=ndi.label(interface,structure=np.ones((3,3,3)))
    if not n:raise TopBrainError('MCA_ICA_ANCHOR_MISSING',f'{case.paths.case_id}/{side}')
    counts=np.bincount(patches.ravel())[1:];order=np.argsort(counts)[::-1]
    if n>1 and counts[order[1]]>=.2*counts[order[0]]:
        raise TopBrainError('MCA_ICA_ANCHOR_AMBIGUOUS',f'{case.paths.case_id}/{side}',patch_sizes=counts)
    interface=patches==order[0]+1
    slices,origin=crop(union);native_voxels=np.argwhere(skeletonize(np.pad(union[slices],1),method='lee'))+origin-1
    points=physical_points(native_voxels,case.affine_mm)
    native=np.array([values.index(int(v))+1 for v in case.labels[tuple(native_voxels.T)]],dtype=np.int8)
    anchor_points=physical_points(np.argwhere(interface),case.affine_mm)
    d=cKDTree(anchor_points).query(points)[0]
    plateau=np.flatnonzero(np.isclose(d,d.min(),atol=1e-8,rtol=0))
    root_id=plateau[np.argmin(np.linalg.norm(points[plateau]-anchor_points.mean(axis=0),axis=1))]
    root=points[root_id];direction=trunk_direction(points,root)
    ids,assignment,mass=downsample(points,maximum)
    geometry=Geometry(points[ids],mass,root,direction,case.paths.case_id,side)
    metadata=dict(version=VERSION,case_id=case.paths.case_id,side=side,modality='mr',
        affine=case.affine_mm,spacing=case.spacing_mm,source_hashes={k:case.provenance[k] for k in ['image_sha256','label_sha256']},
        source_paths=dict(image=str(case.paths.image),label=str(case.paths.label)),labelmap=case.label_map.report(),
        native_voxel_counts={f'M{i+1}':int(m.sum()) for i,m in enumerate(masks)},points_voxel=native_voxels,
        root_anchor='MCA union / ipsilateral ICA native contact; no internal class boundaries',
        root_voxel=native_voxels[root_id],root_world_mm=root,root_contact_distance_mm=float(d[root_id]),
        contact_patch_sizes=counts,downsampling='geometry-only grid; nearest actual sample to cell mean supplies native class',maximum_samples=maximum)
    return SemanticPoints(geometry,native[ids],points,native,ids,assignment,metadata)


def save(points, path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);g=points.geometry
    np.savez_compressed(path,points=g.points,mass=g.mass,root=g.root,direction=g.direction,labels=points.labels,
                        native_points=points.native_points,native_labels=points.native_labels,
                        representative_ids=points.representative_ids,assignment=points.assignment)
    write_json(path.with_suffix('.json'),points.metadata)


def load(path):
    a=np.load(path);meta=json.loads(Path(path).with_suffix('.json').read_text())
    g=Geometry(a['points'],a['mass'],a['root'],a['direction'],meta['case_id'],meta['side'])
    return SemanticPoints(g,a['labels'],a['native_points'],a['native_labels'],a['representative_ids'],a['assignment'],meta)


def build_pool(root, output, maximum=275):
    folder=Path(output)/'semantic_cache';folder.mkdir(parents=True,exist_ok=True)
    pool=[];excluded=[]
    for pair in discover(root,'2025'):
        if pair.modality!='mr':continue
        case=None
        for side in ['L','R']:
            path=folder/f'MRA{pair.case_id}_{side}.npz'
            if path.exists():
                pointset=load(path);meta=pointset.metadata
                if meta['version']==VERSION and meta['maximum_samples']==maximum and all(
                    sha256(getattr(pair,k))==meta['source_hashes'][k+'_sha256'] for k in ['image','label']):
                    pool.append(pointset);continue
            if case is None:case=load_case(root,pair.case_id,'mr',version='2025')
            try:pointset=extract(case,side,maximum)
            except TopBrainError as e:
                excluded.append(dict(case_id=pair.case_id,side=side,status=e.code,details=e.details));continue
            save(pointset,path);pool.append(pointset)
        print(f'Semantic cache MRA{pair.case_id}: {len(pool)} valid sides, {len(excluded)} exclusions',flush=True)
    write_json(Path(output)/'semantic_pool.json',dict(available=[dict(case_id=p.geometry.case_id,side=p.geometry.side,samples=len(p.labels)) for p in pool],excluded=excluded))
    return pool,excluded
