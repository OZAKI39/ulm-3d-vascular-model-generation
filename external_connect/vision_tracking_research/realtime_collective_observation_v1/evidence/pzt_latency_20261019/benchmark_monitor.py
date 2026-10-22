"""Read-only full-vs-display timing and value equivalence. Close all AFG GUIs first."""
import json
import statistics
import time
from afg_control import Rig, load_config, BASE
from afg_telemetry import read_monitor_channel

class ReadOnlyRig(Rig):
    def _write(self, session, command, **kwargs):
        if not all(c.endswith('?') for c in command.split(';')):
            raise RuntimeError('Control command prohibited by read-only benchmark')
        return super()._write(session,command,**kwargs)


def main():
    cfg=load_config();cfg['query_batch_mode']='sequential'
    report=dict(mode='QUERY_ONLY',measurements=[],equivalent=True)
    began=time.perf_counter()
    with ReadOnlyRig(cfg) as rig:
        report['devices']=rig.identities.copy();report['before']=rig.status()
        for repeat in range(5):
            for name in cfg['channels']:
                start=time.perf_counter();full=rig.read_channel(name);full_ms=(time.perf_counter()-start)*1000
                start=time.perf_counter();display=read_monitor_channel(rig,name);display_ms=(time.perf_counter()-start)*1000
                equal=display=={k:full[k] for k in display}
                report['measurements'].append(dict(repeat=repeat,channel=name,full_ms=full_ms,display_ms=display_ms,equivalent=equal))
                if not equal:
                    report['equivalent']=False
                    raise RuntimeError('Live settings changed or summary differs; benchmark stopped')
        report['after']=rig.status();report['settings_unchanged']=report['before']==report['after'];report['rig_touched']=rig.touched
    report['wall_s']=time.perf_counter()-began
    report['full_channel_median_ms']=statistics.median(r['full_ms'] for r in report['measurements'])
    report['display_channel_median_ms']=statistics.median(r['display_ms'] for r in report['measurements'])
    report['four_channel_query_estimate_s']=4*report['display_channel_median_ms']/1000
    (BASE/'monitor_benchmark.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['equivalent','settings_unchanged','rig_touched','full_channel_median_ms','display_channel_median_ms','wall_s']}))

if __name__=='__main__':main()
