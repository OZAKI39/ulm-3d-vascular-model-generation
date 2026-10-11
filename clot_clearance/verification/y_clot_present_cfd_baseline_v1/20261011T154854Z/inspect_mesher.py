import sv,json
o=sv.meshing.TetGenOptions(global_edge_size=.00025,surface_mesh_flag=False,volume_mesh_flag=True)
print('VALUES',o.get_values())
print('DIR',dir(o))
print('LOCAL',sv.meshing.TetGenOptions.local_edge_size.__doc__)
print('SPHERE',sv.meshing.TetGenOptions.sphere_refinement.__doc__)
