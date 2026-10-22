import os
os.environ['QT_QPA_PLATFORM']='offscreen'
from pathlib import Path
import sys,time,json,statistics,hashlib
ROOT=Path('/home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1')
sys.path.insert(0,str(ROOT/'src'))
from PySide6 import QtWidgets as W
from PySide6.QtTest import QTest
from ui import Window
from pzt_status import NAMES,ENDPOINTS
out=ROOT/'evidence/pzt_latency_20261019'
app=W.QApplication([])
cfg=json.loads((ROOT/'configs/default.json').read_text());cfg['source_mode']='LIVE_CAMERA'
w=Window(cfg,out/'ui_probe');w.show()
report=dict(query_only=True,manual=[],samples=[])
try:
 deadline=time.monotonic()+5
 while (w.pzt_panel.data.get('controller') in (None,'CLOSED') or not w.pzt_panel.data.get('refresh_supported')) and time.monotonic()<deadline:QTest.qWait(20)
 assert w.pzt_panel.data.get('refresh_supported'),w.pzt_panel.data
 report['initial_controller']=w.pzt_panel.data['controller']
 for repeat in range(4):
  t=time.monotonic();w.pzt_panel.read_button.click()
  assert w.pzt_panel.pending_request,w.pzt_panel.hint.text()
  request_id=w.pzt_panel.pending_request[0]
  while w.pzt_panel.pending_request and time.monotonic()-t<21:QTest.qWait(20)
  elapsed=time.monotonic()-t;ack=w.pzt_panel.data.get('refresh') or {}
  report['manual'].append(dict(repeat=repeat,elapsed_s=elapsed,ack=ack))
  assert ack.get('request_id')==request_id and ack.get('state')=='COMPLETE',ack
  QTest.qWait(300)
 QTest.qWait(1200)
 last={n:r.get('read_unix_s') for n,r in w.pzt_panel.data['channels'].items()}
 start=time.monotonic()
 while time.monotonic()-start<12:
  QTest.qWait(10)
  now=time.time()
  for name,row in w.pzt_panel.data['channels'].items():
   stamp=row.get('read_unix_s')
   if stamp and stamp!=last.get(name):
    report['samples'].append(dict(channel=name,read_unix_s=stamp,seen_unix_s=now,period_s=stamp-last[name],publish_to_panel_s=now-w.pzt_panel.data['published_unix_s'],read_to_panel_s=now-stamp,state=row['state']))
    last[name]=stamp
 report['final_status']=w.pzt_panel.data
 report['cards']={key:text.text() for key,(_,text) in w.pzt_panel.cards.items()}
 report['mapping']=w.pzt_panel.layout_data
 report['layout_sha256']=hashlib.sha256((ROOT/'configs/pzt_layout.json').read_bytes()).hexdigest()
 report['screenshot']=str(out/'optimized_actual_readback.png');w.grab().save(report['screenshot'])
 samples=report['samples'];report['median_period_s']=statistics.median(s['period_s'] for s in samples)
 report['median_read_to_panel_s']=statistics.median(s['read_to_panel_s'] for s in samples)
 report['p95_read_to_panel_s']=sorted(s['read_to_panel_s'] for s in samples)[int(.95*(len(samples)-1))]
 report['median_warm_button_s']=statistics.median(s['elapsed_s'] for s in report['manual'][1:])
 assert all(s['state']=='READ_OK' for s in samples)
 assert all(r['parameters'] for r in report['final_status']['channels'].values())
 report['passed']=True
finally:
 (out/'ui_latency_probe.json').write_text(json.dumps(report,indent=2))
 w.close();app.processEvents()
print(json.dumps({k:report[k] for k in ['passed','initial_controller','median_period_s','median_read_to_panel_s','p95_read_to_panel_s','median_warm_button_s']}))
