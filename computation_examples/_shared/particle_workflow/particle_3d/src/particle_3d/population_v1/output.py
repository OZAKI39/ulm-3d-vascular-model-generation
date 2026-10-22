"""Chunked full accepted-time state, source ledger and explicit terminal records."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import gzip,hashlib,json,time
import numpy as np
from .provenance import canonical,write_json,digest


class PopulationOutput:
    def __init__(self,folder,identity,*,chunk_rows=200000):
        self.folder=Path(folder);self.folder.mkdir(parents=True,exist_ok=False)
        (self.folder/'states').mkdir();self.identity=dict(identity);self.chunk_rows=int(chunk_rows)
        self.pending=[];self.pending_rows=0;self.chunks=[];self.index=0;self.futures=[]
        self.executor=ThreadPoolExecutor(2);self.io_seconds=0.;self.total_rows=0
        self.streams={name:gzip.open(self.folder/(name+'.jsonl.gz'),'wb',compresslevel=3)
                      for name in ['source','terminals','interaction_audit']}
        self.hashes={name:hashlib.sha256() for name in self.streams}
        write_json(self.folder/'IDENTITY.json',self.identity)

    def _row(self,name,row):
        data=canonical(row);self.streams[name].write(data);self.hashes[name].update(data)
    def source(self,row):self._row('source',row)
    def terminal(self,row):self._row('terminals',row)

    def snapshot(self,row):
        copied={k:np.asarray(v).copy() if isinstance(v,np.ndarray) else v for k,v in row.items()}
        self.pending.append(copied);self.pending_rows+=len(row['particle_id'])
        if self.pending_rows>=self.chunk_rows:self.flush()

    def audit(self,world):
        a=world.projection
        pair=[r for r in a.get('potential_blocks',[]) if r['kind']=='PAIR' and r['active']]
        contact=a.get('handoff_constraints',[])
        components=a.get('components',[a])
        # Save all interacting IDs and actual numerical constraint diagnostics.
        row=dict(physical_time_s=world.time_s,pair_nearfield=pair,handoff_contacts=contact,
            solver=[dict(particle_ids=r.get('particle_ids',sorted(world.shapes)),
                contact_count=r.get('contact_count',0),contact_ids=r.get('contact_ids',[]),
                multipliers=r.get('multipliers',[]),contact_kkt=r.get('contact_kkt',{}),
                contact_redundancy=r.get('contact_redundancy',{}),
                max_constraint_violation_m_s=r.get('max_constraint_violation_m_s',0.),
                backward_residual=r.get('relative_backward_residual',0.))
                for r in components if r.get('contact_count',0)],
            handoff_event=a.get('handoff_event'),accepted_certificate_count=len(a.get('continuous_certificates',[])))
        self._row('interaction_audit',row)

    def _write_chunk(self,index,rows):
        start=time.perf_counter();arrays={};keys=[k for k in rows[0] if k!='physical_time_s']
        offsets=np.r_[0,np.cumsum([len(r['particle_id']) for r in rows])]
        for k in keys:arrays[k]=np.concatenate([r[k] for r in rows],axis=0)
        arrays['state_time_s']=np.array([r['physical_time_s'] for r in rows],dtype=np.float64)
        arrays['state_offsets']=offsets
        science=hashlib.sha256()
        for k in sorted(arrays):
            a=np.ascontiguousarray(arrays[k]);science.update(canonical(dict(name=k,shape=a.shape,dtype=str(a.dtype))));science.update(a.tobytes())
        path=self.folder/'states'/f'chunk_{index:06d}.npz';np.savez_compressed(path,**arrays)
        return dict(index=index,path=str(path.relative_to(self.folder)),state_count=len(rows),row_count=int(offsets[-1]),
            first_time_s=rows[0]['physical_time_s'],last_time_s=rows[-1]['physical_time_s'],
            array_sha256=science.hexdigest(),file_sha256=digest(path),write_seconds=time.perf_counter()-start)

    def flush(self):
        if not self.pending:return
        self.futures.append(self.executor.submit(self._write_chunk,self.index,self.pending))
        self.index+=1;self.total_rows+=self.pending_rows;self.pending=[];self.pending_rows=0
        if len(self.futures)>=2:self.chunks.append(self.futures.pop(0).result())

    def close(self,summary):
        self.flush();self.chunks.extend(f.result() for f in self.futures);self.executor.shutdown()
        for stream in self.streams.values():stream.close()
        manifest=dict(schema='POPULATION_ACCEPTED_TIME_STATES_V1',identity=self.identity,chunks=self.chunks,
            recorded_particle_states=self.total_rows,io_seconds=sum(c['write_seconds'] for c in self.chunks),
            ledger_scientific_hashes={k:h.hexdigest() for k,h in self.hashes.items()},
            time='ABSOLUTE_PHYSICAL_SECONDS',age='PHYSICAL_TIME_MINUS_BIRTH_TIME',
            interpolation='HELD_VELOCITY_LINEAR_CENTER_PATH_BETWEEN_EVERY_ACCEPTED_ENDPOINT',
            velocity='HELD_ON_PRECEDING_ACCEPTED_INTERVAL; BIRTH_ROW_IS_INITIAL_FREE_VELOCITY',
            contact_ids='interaction_audit.jsonl.gz',**summary)
        write_json(self.folder/'MANIFEST.json',manifest)
        return manifest


def iter_states(folder):
    folder=Path(folder);manifest=json.loads((folder/'MANIFEST.json').read_text())
    for chunk in manifest['chunks']:
        with np.load(folder/chunk['path']) as data:
            for k,t in enumerate(data['state_time_s']):
                a,b=data['state_offsets'][k:k+2]
                yield dict(physical_time_s=float(t),**{name:data[name][a:b] for name in data.files if name not in ['state_time_s','state_offsets']})
