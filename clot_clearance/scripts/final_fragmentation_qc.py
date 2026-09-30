"""Delivery audit: preserved inputs, source identity, VTK, media and local links."""
import hashlib,json,os,re
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote,urlsplit
import numpy as np
import pyvista as pv
import imageio.v2 as imageio
import imageio_ffmpeg
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'results/streaming_fragmentation_demo'
VIS=ROOT/'visualization/streaming_fragmentation_demo'
PROV=ROOT/'provenance/fragmentation_extension'

def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()

def jsonread(path):return json.loads(path.read_text())

def preserved(manifest,relative_root):
    entries=jsonread(manifest);bad=[]
    for name,expected in entries.items():
        path=relative_root/name
        if not path.is_file() or sha(path)!=expected:bad.append(name)
    return dict(status='PASS' if not bad else 'FAIL',file_count=len(entries),changed_or_missing=bad)

class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ['href','src','poster'] and value:self.links.append(value)

def main():
    checks={};details={}
    details['preexisting_clot_files']=preserved(PROV/'BASELINE_FILES_SHA256.json',ROOT)
    checks['preexisting_clot_files']=details['preexisting_clot_files']['status']=='PASS'
    manifest=ROOT/'provenance/ORIGINAL_FILES_SHA256.json'
    original=preserved(manifest,Path('/'));handoff=Path('/home/lzy/projects/BRAVA_CONTEXT_HANDOFF.md')
    unexpected=[p for p in original['changed_or_missing'] if p!=str(handoff)]
    checks['original_BraVa_and_source_files_except_shared_handoff']=not unexpected
    shared_append=False
    if str(handoff) in original['changed_or_missing']:
        body=handoff.read_bytes();expected=jsonread(manifest)[str(handoff)]
        # Preserve evidence of a concurrent append without accepting arbitrary edits.
        candidates=[i+1 for i,c in enumerate(body) if c==10]+[len(body)]
        matching=[i for i in candidates if hashlib.sha256(body[:i]).hexdigest()==expected]
        shared_append=bool(matching)
        details['shared_handoff_change']=dict(path=str(handoff),original_sha256=expected,current_sha256=sha(handoff),
            preserved_original_prefix_bytes=matching,appended_bytes=len(body)-max(matching) if matching else None,
            observed_mtime_ns=handoff.stat().st_mtime_ns,
            explanation='Concurrent appended mesh_particle_validation_v2 progress; this extension does not write the shared handoff.')
        (PROV/'SHARED_HANDOFF_OBSERVED.md').write_bytes(body)
        checks['shared_handoff_original_content_preserved']=shared_append
    details['original_BraVa_and_source_files']=original
    details['original_BraVa_and_source_files']['status']='PASS_WITH_SHARED_CONTEXT_APPEND' if shared_append and not unexpected else original['status']
    identity=jsonread(RUN/'IDENTITY.json');mismatch=[]
    for name,expected in identity['source_sha256'].items():
        for path in [ROOT/name,PROV/'source_v1'/name]:
            if not path.is_file() or sha(path)!=expected:mismatch.append(str(path))
    checks['run_source_and_snapshot_match']=not mismatch;details['source_mismatch']=mismatch
    config_file=ROOT/'configs/streaming_fragmentation_demo.json'
    canonical_hash=hashlib.sha256(json.dumps(jsonread(config_file),sort_keys=True).encode()).hexdigest()
    checks['selected_config_matches_run']=sha(config_file)==sha(RUN/'CONFIG.json') and canonical_hash==identity['config_sha256']
    checks['mild_config_preserved']=sha(ROOT/'configs/straight_pipe_damage_demo.json')==sha(ROOT/'configs/straight_pipe.json')
    checks['mild_result_alias_correct']=(ROOT/'results/straight_pipe_damage_demo').resolve()==(ROOT/'results/straight_pipe_demo').resolve()
    checks['fragmentation_result_alias_correct']=RUN.resolve()==ROOT/'runs/fragmentation_pilot_003'
    for name,path in [
        ('independent_audit',RUN/'INDEPENDENT_AUDIT.json'),
        ('transport_force_attribution',RUN/'TRANSPORT_FORCE_ATTRIBUTION.json'),
        ('mild_bitwise_reproduction',ROOT/'verification/fragmentation_extension/MILD_REPRODUCTION.json')]:
        checks[name]=jsonread(path)['status']=='PASS'
    testlog=(ROOT/'logs/fragmentation/ctest_final.xml').read_text()
    checks['all_26_tests_passed']='26 passed' in testlog and 'failures="0"' in testlog
    z=np.load(RUN/'states.npz');vtk_ok=True;vtk_count=0
    for k in range(len(z['cycles'])):
        p=pv.read(RUN/f'vtk/particles_{k:04d}.vtp');b=pv.read(RUN/f'vtk/bonds_{k:04d}.vtp');vtk_count+=2
        vtk_ok &= p.n_points==360 and b.n_cells==12592
        vtk_ok &= np.array_equal(p.points,z['x'][k]) and np.array_equal(b.points,z['x'][k])
        for field,stored in [('damage','damage'),('fragment_id','fragment_id'),('attached_to_base','attached'),('is_surface_particle','surface'),('velocity','v'),('deformation_rank','rank')]:
            vtk_ok &= field in p.point_data and np.array_equal(p[field],z[stored][k])
        vtk_ok &= np.array_equal(p['displacement'],z['x'][k]-z['X'])
        for field,stored in [('integrity','integrity'),('accumulated_D','bond_D'),('active','bond_active'),('recently_broken','recent_broken'),('measured_cycle_amplitude','amplitude')]:
            vtk_ok &= field in b.cell_data and np.array_equal(b[field],z[stored][k])
    checks['all_VTK_fields_match_saved_states']=bool(vtk_ok);details['vtk_files_read']=vtk_count
    render=jsonread(VIS/'RENDER_MANIFEST.json');clip_ok=True
    for axis,name in enumerate(['x','y','z']):
        lo,hi=render['camera'][name+'lim_mm'];values=z['x'][:,:,axis]*1e3
        clip_ok &= values.min()>=lo and values.max()<=hi
    checks['all_particle_positions_inside_fixed_camera_limits']=bool(clip_ok)
    checks['true_displacement_and_equal_units']=render['displacement_scale']==1 and render['equal_axis_units']
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe()
    decode=PROV/'decoded_event_frames';decode.mkdir(exist_ok=True)
    media=[]
    for focus in ['damage','fragments']:
        movie=VIS/f'fragmentation_{focus}.mp4';reader=imageio.get_reader(movie)
        meta=reader.get_meta_data();count=0;shape_ok=True;event_image_errors={}
        selected={i*6:i for i in [0,2,11,25]}
        for i,frame in enumerate(reader):
            count+=1;shape_ok &= frame.shape==(900,1600,3)
            if i in selected:
                state=selected[i]
                imageio.imwrite(decode/f'{focus}_state_{state:04d}.png',frame)
                reference=imageio.imread(VIS/f'frames/{focus}/frame_{state:04d}.png')[:,:,:3]
                event_image_errors[str(state)]=float(np.abs(frame.astype(float)-reference).mean())
        reader.close()
        with Image.open(VIS/f'fragmentation_{focus}.gif') as gif:
            gif_count=gif.n_frames
            gif_duration_ms=0
            for i in range(gif_count):gif.seek(i);gif.load();gif_duration_ms+=gif.info.get('duration',0)
        good=count==156 and shape_ok and meta['fps']==12 and gif_count==26 and gif_duration_ms==13000 and max(event_image_errors.values())<3
        checks[f'{focus}_movie_and_gif_full_decode']=bool(good)
        media.append(dict(name=focus,decoded_mp4_frames=count,fps=meta['fps'],size=list(meta['size']),duration_s=meta['duration'],gif_frames=gif_count,gif_duration_ms=gif_duration_ms,event_frame_mean_absolute_pixel_error_0to255=event_image_errors))
    details['media']=media
    parser=Links();parser.feed((ROOT/'OPEN_FRAGMENTATION.html').read_text());badlinks=[]
    for url in parser.links:
        parts=urlsplit(url)
        if not parts.scheme and parts.path and not (ROOT/unquote(parts.path)).exists() and ROOT/unquote(parts.path)!=PROV/'FINAL_QC.json':badlinks.append(url)
    checks['offline_portal_static_links_exist']=not badlinks;details['bad_links']=badlinks
    checks['all_slider_frames_exist']=all((VIS/f'frames/damage/frame_{i:04d}.png').is_file() for i in range(26))
    # Markdown links belong to newly added documentation only.
    missing_doc_links=[]
    for name in ['README_FRAGMENTATION.md','FRAGMENTATION_MODEL.md','FRAGMENTATION_REPORT_ZH.md']:
        for url in re.findall(r'\]\(([^)]+)\)',(ROOT/name).read_text()):
            if not urlsplit(url).scheme and not (ROOT/url.split('#')[0]).exists() and ROOT/url.split('#')[0]!=PROV/'FINAL_QC.json':missing_doc_links.append([name,url])
    checks['new_document_links_exist']=not missing_doc_links;details['missing_doc_links']=missing_doc_links
    result=dict(status=('PASS_WITH_SHARED_CONTEXT_APPEND' if shared_append else 'PASS') if all(checks.values()) else 'FAIL',checks=checks,details=details,
        visual_inspection=dict(source_frames_reviewed=[0,2,11,25],figures_reviewed=['four_panel_summary','time_histories'],
            decoded_frames_reviewed=['damage_state_0000','damage_state_0002','damage_state_0011','fragments_state_0025'],
            observation='Fixed view shows initially connected material, newly broken red bonds, first detached IDs and downstream separation at true coordinates; decoded event frames inspected and compared with source PNGs.'),
        browser_interaction_automated_tested=False,
        limitations='Verification only; 172 isolated particles and explicit rank-deficient continuation; no quantitative physical validation.')
    (PROV/'FINAL_QC.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not all(checks.values()):raise SystemExit(1)

if __name__=='__main__':main()
