"""Generate isolated V4 surface sources from the frozen V3/V2 implementations."""
from common import *
import shutil
def replace(text,old,new):
 assert text.count(old)==1,(old[:100],text.count(old));return text.replace(old,new)
def main():
 shutil.copy2(V3_A/'shear_field_library.py',A/'shear_field_library.py')
 t=(V2_A/'current_surface.py').read_text()
 t=replace(t,"X=ref['X']; axes=[np.unique(X[:,j]) for j in range(3)]\n        lattice=np.rint((X-X.min(0))/DX).astype(int)","X=ref['X']; grid=ref['regular_X']; axes=[np.unique(grid[:,j]) for j in range(3)]\n        lattice=np.rint((grid-grid.min(0))/DX).astype(int)")
 (A/'supported_surface.py').write_text(t)
 t=(V3_A/'current_surface.py').read_text();t=t.replace("OldSurface=module('v2_surface_frozen',V2_A/'current_surface.py').CurrentSurface","from supported_surface import CurrentSurface as OldSurface\nimport mapped_helpers as mapped")
 t=replace(t,'self.reference_area=np.r_[self.base.reference_area,np.full(self.M,DX**2)]','self.reference_area=np.r_[self.base.reference_area,np.full(self.M,DX**2)]\n        mapped.initialize(self)')
 t=replace(t,"self.lattice=np.rint((ref['X']-ref['X'].min(0))/DX).astype(int)","grid=ref['regular_X'];self.lattice=np.rint((grid-grid.min(0))/DX).astype(int)")
 t=replace(t,"self.lattice[i]+(self.offset[k]/(DX/2)).astype(int)","self.lattice[i]+np.eye(3,dtype=int)[self.axis[k]]*int(self.sign[k])")
 lo=t.index('    def _load_ray_kernel(self):');hi=t.index('    def topology(self,d):')
 t=t[:lo]+'''    def _load_ray_kernel(self):mapped.load_kernel(self)
    def _reference_crossings(self):return mapped.reference_crossings(self)

'''+t[hi:]
 old='''            i=self.owner[face];a=self.axis[face];g=gradients[a];lo,hi=g.indptr[i:i+2]
            js=g.indices[lo:hi];vals=.5*self.sign[face]*g.data[lo:hi]
            js=np.r_[js,i];vals=np.r_[vals,1.]'''
 new='''            i=self.owner[face];js=[i];vals=[1.]
            for a,g in enumerate(gradients):
                lo,hi=g.indptr[i:i+2];js.extend(g.indices[lo:hi]);vals.extend((self.offset[face,a]/DX*g.data[lo:hi]).tolist())
            js=np.asarray(js);vals=np.asarray(vals)'''
 t=replace(t,old,new)
 lo=t.index('        normal=np.zeros((M,3));area=np.zeros(M)');hi=t.index('        geometric=valid[owners]',lo)
 t=t[:lo]+'''        av=np.linalg.det(deformation[owners])[:,None]*np.einsum('nji,nj->ni',inverse[owners],self.reference_area_vectors)
        area=np.linalg.norm(av,axis=1);normal=av/np.maximum(area[:,None],1e-300)
'''+t[hi:]
 t=t.replace('np.sqrt(3)*DX/2+1e-12','self.cell_reference_radius.max()+1e-12')
 lo=t.index('        nb=np.maximum(self.neighbor,0);half_self=');hi=t.index('        ray_ids=np.flatnonzero(opened);',lo)
 t=t[:lo]+'''        nb=np.maximum(self.neighbor,0)
        corner_offset=np.einsum('nij,nkj->nki',deformation,self.cell_corner_offset)
        half_self=np.einsum('nkj,nj->nk',corner_offset[owners],normal).max(1)
        half_other=-np.einsum('nkj,nj->nk',corner_offset[nb],normal).min(1)
        half_other=np.where(valid[nb],half_other,self.cell_reference_radius[nb])
        aperture=np.sum((x[nb]-x[owners])*normal,axis=1)-half_self-half_other
        original=self.neighbor<0;aperture[original]=FULL_OPEN
        opened=eligible&(aperture>OPEN)
        radii=np.linalg.norm(corner_offset,axis=2).max(1);radii[~valid]=self.cell_reference_radius[~valid]
'''+t[hi:]
 lo=t.index('            inputs=[np.ascontiguousarray(center[ray_ids])');hi=t.index('            # Also respect',lo)
 t=t[:lo]+'''            inputs,cl,bl=mapped.mapped_ray(self,ray_ids,center,normal,owners,x,deformation,valid,radii)
'''+t[hi:]
 t=t.replace('area_ratio=area/DX**2','area_ratio=area/self.reference_area[self.L:]')
 t=t.replace("sha(V2_O/'provenance/reference.npz')","sha(O/'provenance/reference.npz')").replace('V3_DUAL_FACE_1','V4_MAPPED_DUAL_FACE_1')
 t=t.replace('use reference dx^3 cells','use actual nonuniform mapped reference cells').replace('These fixed incidence lists are geometry independent.','These lists are frozen only after the V4 reference geometry is fixed.')
 (A/'current_surface.py').write_text(t)
 print('V4_SURFACE_SOURCES_GENERATED',flush=True)
if __name__=='__main__':main()
