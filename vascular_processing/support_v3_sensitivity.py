"""Traceable full-rib width/gap variants; no physical calibration inference."""
import pickle
from pathlib import Path
from .support_v3_geometry import curtain,union
from .support_inputs import save,table


def generate(record,profile,out):
    folder=Path(record['folder']);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    variants=[]
    for kind,values in [('line_count',[1,2,3]),('roof_gap_mm',[.15,.20,.25])]:
        for value in values:
            lines=int(value) if kind=='line_count' else record['line_count']
            gap=float(value) if kind=='roof_gap_mm' else .20
            built=[]
            for r in ribs:
                built.append(curtain(r['path_xy'],r['top']+.20-gap,r['base'],profile,lines,r['roof_width_mm'],r['segment_length_mm'],3.)['mesh'])
            mesh=union(built);name=f'{kind}_{value:g}_support_only.stl';mesh.export(out/name)
            variants.append(dict(kind=kind,value=value,file=name,volume_mm3=float(mesh.volume),purpose='SENSITIVITY_ONLY_NOT_PRODUCTION',
                source_candidate=record['candidate_id'],width_mm=lines*profile['line_width'],roof_gap_mm=gap,
                actual_slice_validated=False,physical_status='NOT PHYSICALLY VALIDATED'))
    table(out/'variants.csv',variants);save(out/'manifest.json',dict(variants=variants,bridge_policy='NO_UNCALIBRATED_BRIDGE',
        bridge_2_4_6_8mm_reference=str(out.parents[1]/'v2/coupon'),physical_coupon_available=False))
    return variants
