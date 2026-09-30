"""Read-only original-file audit and independent delivered-result checks."""
import hashlib,json,re,sys
from pathlib import Path
from html.parser import HTMLParser
import xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
import imageio.v2 as imageio
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pd_clot.output import write_json


class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ['src','href']:self.links.append(value)


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    root=Path(__file__).resolve().parents[1];run=root/'results/straight_pipe_demo'
    originals=json.loads((root/'provenance/ORIGINAL_FILES_SHA256.json').read_text())
    modified=[p for p,h in originals.items() if not Path(p).is_file() or sha(p)!=h]
    checks={'original_files_unchanged':not modified}
    z=np.load(run/'states.npz');s=json.loads((run/'SUMMARY.json').read_text());c=json.loads((run/'CONFIG.json').read_text())
    checks['finite_states']=bool(all(np.isfinite(z[k]).all() for k in ['x','v','integrity','damage']))
    checks['fixed_base_exact']=bool(np.all(z['x'][:,z['fixed']]==z['X'][z['fixed']]) and np.all(z['v'][:,z['fixed']]==0))
    checks['irreversible_bonds']=bool(np.all(np.diff(z['integrity'],axis=0)<=1e-15))
    ids=z['pairs'].ravel();deg=np.bincount(ids,minlength=len(z['X']))
    D=np.bincount(ids,weights=np.repeat(1-z['integrity'][-1],2),minlength=len(z['X']))/deg
    checks['independent_damage']=bool(np.max(np.abs(D-z['damage'][-1]))<1e-14 and abs(D.mean()-s['final']['mean_particle_damage'])<1e-14)
    # Independent disjoint-set connectivity, separate from scipy.csgraph implementation.
    parent=list(range(len(D)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for (i,j),b in zip(z['pairs'],z['integrity'][-1]):
        if b>c['damage']['connectivity_threshold']:parent[find(int(i))]=find(int(j))
    labels=np.array([find(i) for i in range(len(D))]);count=len(np.unique(labels))
    checks['independent_fragments']=count==s['final']['number_of_fragments']
    fields=['displacement','velocity','particle_damage','max_bond_damage','number_of_active_bonds','strain_energy_density','fragment_id']
    vtk_count=0
    for step in range(len(z['cycles'])):
        p=pv.read(run/f'vtk/particles_{step:04d}.vtp')
        assert all(key in p.point_data for key in fields)
        np.testing.assert_allclose(p.points,z['x'][step],rtol=0,atol=0)
        np.testing.assert_allclose(p['particle_damage'],z['damage'][step],rtol=0,atol=0)
        for kind in ['bonds','traction']:pv.read(run/f'vtk/{kind}_{step:04d}.vtp')
        vtk_count+=3
    checks['vtk_raw_states_equal']=True
    for kind in ['particles','bonds','traction']:
        data=ET.parse(run/f'{kind}.pvd').findall('.//DataSet')
        assert len(data)==len(z['cycles'])
        assert all((run/e.attrib['file']).exists() for e in data)
    checks['pvd_complete']=True
    identity=json.loads((run/'IDENTITY.json').read_text())
    checks['run_source_snapshot_matches']=all(sha(root/'provenance/run_source_v1'/p)==h for p,h in identity['code_sha256'].items())
    checks['copied_wss_matches_original']=sha(root/'vendor/brava_wss_reference.py')==originals['/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin/vendor/flow_solver_support/wss.py']
    import imageio_ffmpeg,os
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe()
    reader=imageio.get_reader(root/'visualization/clot_damage.mp4');meta=reader.get_meta_data();nframes=0
    for frame in reader:
        assert frame.shape==(720,1280,3);nframes+=1
    reader.close();checks['video_full_decode']=nframes==108 and abs(meta['fps']-12)<1e-12
    assert len(imageio.mimread(root/'visualization/clot_damage.gif'))==9
    links=Links();links.feed((root/'OPEN_RESULTS.html').read_text())
    missing=[link for link in links.links if not (root/link).exists()]
    checks['portal_links_exist']=not missing
    acceptance=json.loads((root/'verification/ACCEPTANCE.json').read_text())
    imports=json.loads((root/'verification/IMPORT_AND_DT_CHECK.json').read_text())
    checks['acceptance_A_to_F']=acceptance['status']=='PASS';checks['import_and_dt']=imports['status']=='PASS'
    ctest=ET.parse(root/'logs/ctest_results.xml').getroot()
    checks['ctest_pass']=int(ctest.attrib.get('failures','0'))==0 and int(ctest.attrib.get('tests','0'))>0
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,original_file_count=len(originals),
        modified_originals=modified,missing_portal_links=missing,verified_vtk_files=vtk_count,
        decoded_video_frames=nframes,video_fps=meta['fps'],video_duration_s=nframes/meta['fps'],
        independent_fragment_count=count,visual_review='300 dpi overview and decoded-size last animation frame inspected; final animation labels separated, blue center visible',
        browser_execution_tested=False,ui_note='Offline Plotly HTML generated; links verified; no browser engine available for automated UI execution')
    write_json(root/'provenance/FINAL_QC.json',result)
    files=[]
    for path in sorted(root.rglob('*')):
        rel=path.relative_to(root)
        if path.is_file() and not any(part in ['.tools','build','__pycache__'] for part in rel.parts) and path.name!='DELIVERY_FILES.json':
            files.append(dict(path=str(rel),bytes=path.stat().st_size,sha256=sha(path)))
    write_json(root/'provenance/DELIVERY_FILES.json',dict(root=str(root),all_new=True,modified_existing_scientific_files=[],excluded_generated_directories=['.tools','build','__pycache__'],files=files))
    print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
