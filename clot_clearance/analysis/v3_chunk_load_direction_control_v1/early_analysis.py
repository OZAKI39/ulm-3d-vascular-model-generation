"""Read only the immutable accepted 0.8 ms checkpoint, not open live streams."""
from common import *
def main():
 folder=O/'early_diagnostics';folder.mkdir(exist_ok=True);state=load_npz(O/'early_checkpoint/accepted.state.npz');assert abs(float(state['time'])-EARLY)<1e-15
 assert json.loads((O/'early_checkpoint/RESTORE_AUDIT.json').read_text())['status']=='PASS'
 ref=load_npz(O/'provenance/reference.npz');audit=module_from_file('early_unchanged_candidate_screen',PREV_A/'topology_audit.py');audit.reference=lambda:ref;audit.keyframes=lambda:iter([(None,state)]);audit.O=folder;audit.main()
 record=json.loads((folder/'V3_CHUNK_CANDIDATES.json').read_text());dump(folder/'ACCEPTED_CHECKPOINT_SOURCE.json',dict(state_sha256=sha(O/'early_checkpoint/accepted.state.npz'),step=int(state['step']),time_s=float(state['time']),classifier_sha256=sha(PREV_A/'topology_audit.py'),candidate_count=record['candidate_count'],candidate_ge14=sum(c['nodes']>=14 for c in record['candidates']),candidate_ge34=record['candidate_ge200'],candidate_ge114=record['candidate_ge300'],no_dynamics_or_parameter_change=True))
if __name__=='__main__':main()
