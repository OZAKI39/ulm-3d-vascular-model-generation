"""Read native numeric XDMF/HDF5 using the existing HDF5 C library, CPU only.

No package installation, MPI initialization, solver import, or write access.
"""
import ctypes as C
import ctypes.util
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def read_xmf(path):
    path = Path(path)
    name = C.util.find_library('hdf5_serial_hl') or C.util.find_library('hdf5_openmpi_hl')
    if not name:
        raise RuntimeError('EXISTING_HDF5_READER_UNAVAILABLE')
    lib = C.CDLL(name)
    lib.H5Fopen.argtypes = [C.c_char_p, C.c_uint, C.c_longlong]
    lib.H5Fopen.restype = C.c_longlong
    lib.H5Fclose.argtypes = [C.c_longlong]
    lib.H5LTget_dataset_ndims.argtypes = [C.c_longlong, C.c_char_p, C.POINTER(C.c_int)]
    lib.H5LTget_dataset_info.argtypes = [C.c_longlong, C.c_char_p, C.POINTER(C.c_ulonglong), C.POINTER(C.c_int), C.POINTER(C.c_size_t)]
    lib.H5LTread_dataset_double.argtypes = [C.c_longlong, C.c_char_p, C.c_void_p]
    result = {}
    for node in ET.parse(path).iter('DataItem'):
        if node.get('Format') != 'HDF':
            continue
        filename, dataset = node.text.strip().split(':', 1)
        if Path(filename).name != filename or filename in ('.', '..'):
            raise ValueError('UNSAFE_XMF_REFERENCE')
        target = path.parent / filename
        if target.is_symlink():
            raise ValueError('SYMLINK_HDF5')
        fid = lib.H5Fopen(str(target).encode(), 0, 0)
        if fid < 0:
            raise ValueError('HDF5_NOT_CLOSED_OR_READABLE')
        try:
            rank = C.c_int()
            key = dataset.encode()
            if lib.H5LTget_dataset_ndims(fid, key, C.byref(rank)) < 0 or not 1 <= rank.value <= 5:
                raise ValueError('BAD_HDF5_DIMENSIONS')
            dims = (C.c_ulonglong * rank.value)(); kind = C.c_int(); size = C.c_size_t()
            if lib.H5LTget_dataset_info(fid, key, dims, C.byref(kind), C.byref(size)) < 0:
                raise ValueError('BAD_HDF5_DATASET')
            shape = tuple(dims)
            if shape != tuple(map(int, node.get('Dimensions').split())) or np.prod(shape) > 5_000_000:
                raise ValueError('HDF5_SHAPE_OR_SIZE_LIMIT')
            array = np.empty(shape, dtype=np.float64)
            if lib.H5LTread_dataset_double(fid, key, array.ctypes.data) < 0:
                raise ValueError('HDF5_READ_FAILED')
            if dataset.strip('/') == 'id':
                if not np.isfinite(array).all() or np.any(np.abs(array) >= 2**53) or np.any(array != np.rint(array)):
                    raise ValueError('INEXACT_PARTICLE_ID')
                array = array.astype(np.int64).reshape(-1)
            result[dataset.strip('/')] = array
        finally:
            lib.H5Fclose(fid)
    if not result:
        raise ValueError('NO_HDF5_DATA')
    return result
