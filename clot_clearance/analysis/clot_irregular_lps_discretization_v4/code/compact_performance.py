"""Small performance delivery without duplicating full fracture-event payloads."""
from common import *
def main():
 source=O/'PERFORMANCE_BASELINE_V3.json';full=json.loads(source.read_text());summary_keys=['status','accepted_steps','predictive_rejections','predictive_substeps','wall_s','cpu_s','peak_RSS_bytes','max_relative_energy_residual','time_s','native_topology_matches_current_native','original_onsets_unchanged']
 small={k:v for k,v in full.items() if k!='windows'};small['windows']={}
 for name,window in full['windows'].items():
  out={k:v for k,v in window.items() if k!='summary'};out['summary']={k:window['summary'][k] for k in summary_keys if k in window['summary']};out['summary']['secondary_event_count']=len(window['summary'].get('secondary_events',[]));small['windows'][name]=out
 small['full_payload_archive']=dict(path=str(source),sha256=sha(source),bytes=source.stat().st_size,omitted_from_compact='Duplicated full fracture-event coordinates, bond arrays, and lineage events; timing fields and numerical-equivalence verdicts retained. Full original remains on server.')
 dump(O/'PERFORMANCE_BASELINE_V3_COMPACT.json',small);print('COMPACT_PERFORMANCE', (O/'PERFORMANCE_BASELINE_V3_COMPACT.json').stat().st_size,flush=True)
if __name__=='__main__':main()
