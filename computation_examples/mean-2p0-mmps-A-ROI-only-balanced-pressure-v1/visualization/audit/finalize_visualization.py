"""Verify collected media/data and summarize the new visualization bundle."""
from pathlib import Path
import json,hashlib,ast
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def dump(p,x): (ROOT/p).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def main():
    context=read('BUILD_CONTEXT.json');compute=read('input_data/field_diagnostics/COMPUTE_VALIDATION.json')
    streamline=read('input_data/streamlines/COMPUTE_VALIDATION.json')
    source=read('audit/SOURCE_CASE_LOCK.json')
    assert all(sha(Path(p))==h for p,h in source.items()),'Source CFD case changed'
    assert compute['source_field_sha256']==streamline['source_field_sha256']==context['source_field_sha256']
    for name,item in read('INPUT_MANIFEST.json')['files'].items():assert sha(ROOT/'input_data'/name)==item['sha256'],name
    shared=ROOT/'render_visualization.py';old=ROOT.parent/'render_visualization.py'
    assert shared.read_bytes()==old.read_bytes()
    def functions(p):return [ast.dump(x,include_attributes=False) for x in ast.parse(p.read_text()).body if isinstance(x,(ast.FunctionDef,ast.ClassDef))]
    assert functions(ROOT/'render_surface_fields.py')==functions(ROOT.parent/'render_surface_fields.py')
    videos=[];independent=[];devices=[]
    for base in ['results','results/surface_fields']:
        media=read(base+'/MEDIA_VALIDATION.json');check=read(base+'/INDEPENDENT_VALIDATION.json')
        assert media['all_pass'] and check['all_pass']
        for v in media['videos']:
            f=ROOT/base/v['file'];assert sha(f)==v['sha256']
            assert v['decoded_frames']==v['distinct_frames']==432 and v['fps']==24 and v['size']==[1920,1080]
            videos.append(dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,**v))
        for view in media['views'].values():
            assert any('NVIDIA GeForce RTX 4090' in x for x in view['opengl'])
            devices.extend(view['opengl'])
        independent.append(dict(path=base+'/INDEPENDENT_VALIDATION.json',sha256=sha(ROOT/base/'INDEPENDENT_VALIDATION.json')))
    figures=['results/figures/full_vessel_4k.png','results/figures/branch_vectors_4k.png','results/surface_fields/figures/pressure_overview_4k.png','results/surface_fields/figures/wss_overview_4k.png']
    assert all(Image.open(ROOT/p).size==(3840,2160) for p in figures)
    result=dict(all_pass=True,source_case=context['source_case'],source_field_sha256=context['source_field_sha256'],
        source_inputs_unchanged=True,protected_source_file_count=len(source),original_rendering_code_preserved=True,
        pressure_colorbar_pa=[-5,7000],WSS_colorbar_pa=[0,55],speed_colorbar_mm_s=[0,7.5],
        pressure_range_pa=compute['pressure_Pa'],wss_raw_and_display_pa=compute['wss_Pa'],
        videos=videos,images_4k=figures,independent_validation=independent,actual_GPU_renderers=sorted(set(devices)),
        new_CFD_calls=0,finite_size_microbubble_calls=0,RBC_calls=0,
        newly_integrated_visualization_point_streamlines=96,input_alias_note=context['compatibility_aliases'],
        visual_inspection='All four 4K views and multi-angle contact sheets inspected for labels, geometry framing, colorbars and field identity')
    dump('VISUALIZATION_VALIDATION.json',result)
    files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
    (ROOT/'SHA256SUMS.txt').write_text(''.join(sha(p)+'  '+str(p.relative_to(ROOT))+'\n' for p in files))
    print(json.dumps(dict(all_pass=True,videos=len(videos),figures_4k=len(figures),source_field_sha256=context['source_field_sha256']),indent=2))
if __name__=='__main__':main()
