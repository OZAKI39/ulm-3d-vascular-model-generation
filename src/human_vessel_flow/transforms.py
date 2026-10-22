"""Positions use source millimetres, CFD metres, or explicitly named print poses."""
import numpy as np
from .paths import ROOT, read_json

def geometry_record(root=ROOT):
    return read_json(root/'inputs/geometry/geometry.json')

def source_mm_to_cfd_m(points, root=ROOT):
    origin = np.asarray(geometry_record(root)['origin_source_mm'])
    return (np.asarray(points)-origin)*1e-3

def cfd_m_to_candidate_mm(points, candidate, root=ROOT):
    record = geometry_record(root)
    transform = np.asarray(record['pose_transforms'][str(candidate)]['source_mm_to_candidate_mm'])
    source = np.asarray(points)*1e3 + np.asarray(record['origin_source_mm'])
    return source @ transform[:3,:3].T + transform[:3,3]

def rotate_vectors(vectors, candidate, root=ROOT):
    transform = np.asarray(geometry_record(root)['pose_transforms'][str(candidate)]['source_mm_to_candidate_mm'])
    return np.asarray(vectors) @ transform[:3,:3].T

