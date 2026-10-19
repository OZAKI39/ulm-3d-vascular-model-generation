"""Validate delivered figures, video decoding and accepted-state frame provenance."""
from common import *
from PIL import Image
from html.parser import HTMLParser
from urllib.parse import urlparse,unquote
import csv,subprocess,datetime

class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        self.links.extend(v for k,v in attrs if k in ('href','src'))

def main():
    assert json.loads((O/'POST_FINISH_COMPLETE.json').read_text())['status']=='PASS'
    checks=json.loads((O/'DELIVERY_VALIDATION.json').read_text());assert checks['status']=='PASS'
    freeze=json.loads((O/'SOURCE_FREEZE.json').read_text())
    for p,h in freeze['implementation_sha256'].items():assert sha(R/p)==h,p
    for p,h in freeze['input_sha256'].items():assert sha(O/p)==h,p
    request=json.loads((O/'FORMAL_REQUEST.json').read_text())
    assert sha(O/'SOURCE_FREEZE.json')==request['source_freeze_sha256']
    pngs=sorted((F/'figures').glob('FIG_*.png'));assert len(pngs)==8
    figures=[]
    for index,p in enumerate(pngs,1):
        assert p.name.startswith(f'FIG_{index:02d}_')
        with Image.open(p) as im:im.load();width,height=im.size
        assert width>=2500 and height>=1600
        pdf=p.with_suffix('.pdf');raw=pdf.read_bytes()
        assert raw.startswith(b'%PDF-') and b'%%EOF' in raw[-100:]
        figures.append(dict(png=p.name,png_sha256=sha(p),pdf_sha256=sha(pdf),pixels=[width,height]))
    generation=json.loads((O/'VIDEO_GENERATION.json').read_text());assert generation['status']=='PASS'
    names=['VIDEO_01_SURFACE_COUPLING_BEFORE_AFTER','VIDEO_02_FRAGMENT_FLUID_TRACTION_AND_SURVIVAL']
    videos=[];source_cache={};topology_cache={}
    graph=module('artifact_membership_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
    ref=load_npz(O/'provenance/reference.npz')
    for name in names:
        p=F/'videos'/(name+'.mp4');metadata=next(v for v in generation['videos'] if v['file']==p.name)
        assert sha(p)==metadata['sha256']
        probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-show_format','-of','json',str(p)],text=True))
        video=next(v for v in probe['streams'] if v['codec_type']=='video')
        assert (video['width'],video['height'],video['avg_frame_rate'])==(1920,1080,'24/1')
        with p.with_name(name+'_FRAME_TIMES.csv').open() as f:rows=list(csv.DictReader(f))
        assert len(rows)==int(video['nb_frames'])==metadata['frames']
        previous=[-1.,-1.];seen=set()
        for i,row in enumerate(rows):
            assert int(row['frame'])==i and abs(float(row['playback_s'])-i/24)<1e-12
            assert row['state_interpolation']=='False'
            for j,side in enumerate(['left','right']):
                source=R/row[side+'_source'];actual=float(row[side+'_actual_time_s'])
                assert source.name!='failed_unaccepted.npz' and actual>=previous[j]
                assert source.is_relative_to(V2_O/'cases/biphasic/keys') or source.is_relative_to(O/'cases/biphasic/keys')
                if source not in source_cache:
                    with np.load(source) as s:
                        source_cache[source]=dict(path=str(source.relative_to(R)),time_s=float(s['time']),step=int(s['step']),sha256=sha(source))
                assert actual==source_cache[source]['time_s'];previous[j]=actual;seen.add(source)
            if name==names[1]:assert row['left_source']==row['right_source']
            ids=np.asarray(json.loads(row['focus_original_atom_ids']),np.int64)
            assert len(ids)==int(row['focus_nodes']) and np.array_equal(ids,np.unique(ids))
            assert hashlib.sha256(ids.tobytes()).hexdigest()==row['focus_member_ID_sha256']
            source=R/row['right_source']
            if source not in topology_cache:
                with np.load(source) as s:topology_cache[source]=graph.topology(ref,s['d'])[:2]
            labels,sizes=topology_cache[source]
            assert ids.min()>=1 and ids.max()<=int(ref['clot_count'])
            assert len(np.unique(labels[ids-1]))==1 and int(sizes[labels[ids[0]-1]])==len(ids)
            if row['hold']=='True':
                assert i>0 and row['left_source']==rows[i-1]['left_source'] and row['right_source']==rows[i-1]['right_source']
        log=O/(name+'_DECODE.log')
        with log.open('w') as stream:
            subprocess.run(['ffmpeg','-hide_banner','-v','error','-xerror','-i',str(p),'-map','0:v:0','-f','null','-'],stdout=stream,stderr=stream,check=True)
        assert not log.stat().st_size
        with Image.open(p.with_name(name+'_COVER.png')) as im:im.load();assert im.size==(1920,1080)
        videos.append(dict(file=p.name,sha256=sha(p),frames=len(rows),decoded_without_errors=True,source_keys=len(seen),frame_csv_sha256=sha(p.with_name(name+'_FRAME_TIMES.csv'))))
    parser=Links();parser.feed((F/'OPEN_RESULTS.html').read_text())
    for href in parser.links:
        parsed=urlparse(href)
        if not parsed.scheme and parsed.path:assert (F/unquote(parsed.path)).exists(),href
    dump(O/'ARTIFACT_VALIDATION.json',dict(status='PASS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_and_input_freeze_unchanged=True,figures=figures,videos=videos,actual_frame_sources=list(source_cache.values()),
        html_local_links_checked=len(parser.links),visual_review_required_separately=True))
    print('ARTIFACT_VALIDATION_PASS',len(figures),len(videos),len(source_cache),flush=True)

if __name__=='__main__':main()
