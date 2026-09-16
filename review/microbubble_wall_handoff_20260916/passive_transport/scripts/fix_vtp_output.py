from pathlib import Path
import shutil
c=Path(__file__).parent;r=Path('/home/lzy/projects/compre_output/passive_transport_v0/20260915_233516');p=c/'visualize_passive.py';s=p.read_text();marker="def register(name,inputs):\n"
helper='''def write_ascii_vtp(path,coordinates,point_data,polylines=None,polygons=None):
 import xml.etree.ElementTree as ET
 polylines=[] if polylines is None else polylines;polygons=[] if polygons is None else polygons
 root=ET.Element('VTKFile',type='PolyData',version='0.1',byte_order='LittleEndian');pd=ET.SubElement(root,'PolyData');piece=ET.SubElement(pd,'Piece',NumberOfPoints=str(len(coordinates)),NumberOfVerts='0',NumberOfLines=str(len(polylines)),NumberOfStrips='0',NumberOfPolys=str(len(polygons)))
 arrays=ET.SubElement(piece,'PointData')
 for name,values in point_data.items():
  arr=ET.SubElement(arrays,'DataArray',type='Float64',Name=name,format='ascii');arr.text=' '.join(format(float(x),'.17g') for x in values)
 ET.SubElement(piece,'CellData');points=ET.SubElement(piece,'Points');arr=ET.SubElement(points,'DataArray',type='Float64',NumberOfComponents='3',format='ascii');arr.text=' '.join(format(float(x),'.17g') for x in np.asarray(coordinates).ravel())
 for tag,cells in [('Lines',polylines),('Polys',polygons)]:
  group=ET.SubElement(piece,tag);connect=ET.SubElement(group,'DataArray',type='Int64',Name='connectivity',format='ascii');connect.text=' '.join(str(int(i)) for line in cells for i in line);offset=ET.SubElement(group,'DataArray',type='Int64',Name='offsets',format='ascii');offset.text=' '.join(str(int(i)) for i in np.cumsum([len(line) for line in cells]))
 ET.indent(root);ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
 ET.parse(path)
 reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(path));reader.Update();output=reader.GetOutput();assert output.GetNumberOfPoints()==len(coordinates);assert np.array_equal(vtk_to_numpy(output.GetPoints().GetData()),coordinates)
 for key,values in point_data.items():assert np.array_equal(vtk_to_numpy(output.GetPointData().GetArray(key)),values)
 assert output.GetNumberOfLines()==len(polylines) and output.GetNumberOfPolys()==len(polygons)
 return output
'''
s=s.replace(marker,helper+marker)
old="writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(str(out/'FROZEN_LUMEN.vtp'));writer.SetInputData(mesh);assert writer.Write()==1;register('FROZEN_LUMEN.vtp',sources([],['provenance/closed_geometry_m.stl']))"
new="write_ascii_vtp(out/'FROZEN_LUMEN.vtp',vtk_to_numpy(mesh.GetPoints().GetData()).astype(float),{},polygons=tri);register('FROZEN_LUMEN.vtp',sources([],['provenance/closed_geometry_m.stl']))"
assert old in s;s=s.replace(old,new)
a=s.index('def vtp_trajectories(name,rows,inputs):');b=s.index("vtp_trajectories('CASE2_TRAJECTORY.vtp'",a)
s=s[:a]+'''def vtp_trajectories(name,rows,inputs):
 fields=['particle_id','time_s','diameter_um','speed_m_s','fluid_speed_m_s','center_wall_distance_um','surface_wall_gap_um']
 polylines=[np.flatnonzero(rows['particle_id']==pid) for pid in np.unique(rows['particle_id'])]
 write_ascii_vtp(out/name,xyz(rows),{key:rows[key] for key in fields},polylines=polylines)
 register(name,inputs)
'''+s[b:]
p.write_text(s);shutil.copy2(r/'visualization/FROZEN_LUMEN.vtp',r/'provenance/VISUALIZATION_FORMAT_FAILURE/FROZEN_LUMEN.vtp');shutil.copy2(p,r/'scripts/visualize_passive.py')
