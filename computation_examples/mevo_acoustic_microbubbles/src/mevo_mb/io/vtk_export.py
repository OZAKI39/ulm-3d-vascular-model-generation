"""ParaView readable point trajectories, vector forces, acoustic fields and PVD time series."""
from pathlib import Path
from xml.etree.ElementTree import Element, SubElement, ElementTree

import h5py
import numpy as np
import pyvista as pv


def export_run(run: Path, stride: int, field=None):
    root = Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    collection = SubElement(root, "Collection")
    with h5py.File(run/"particles/trajectories.h5") as f:
        times = f["time_s"][:]
        indices = sorted(set(range(0, len(times), stride)) | {len(times)-1})
        for index in indices:
            data = pv.PolyData(f["positions_m"][index])
            data["radius_m"] = f["radii_m"][:]
            data["particle_id"] = f["particle_id"][:]
            for key in ["cluster_labels", "primary_force_N", "secondary_force_N", "flow_velocity_m_s", "total_velocity_m_s"]:
                data[key] = f[key][index]
            name = f"particles_{index:06d}.vtp"
            data.save(run/"vtk"/name)
            SubElement(collection, "DataSet", timestep=str(times[index]), group="", part="0", file=name)
    ElementTree(root).write(run/"vtk/particles.pvd", encoding="utf-8", xml_declaration=True)
    if field is not None:
        grid = pv.RectilinearGrid(*field.axes_m)
        grid["pressure_amplitude_Pa"] = np.abs(field.pressure_Pa).ravel(order="F")
        grid["pressure_phase_rad"] = np.angle(field.pressure_Pa).ravel(order="F")
        grid["pressure_real_Pa"] = field.pressure_Pa.real.ravel(order="F")
        grid["pressure_imag_Pa"] = field.pressure_Pa.imag.ravel(order="F")
        grid.save(run/"vtk/acoustic_field.vtr")
