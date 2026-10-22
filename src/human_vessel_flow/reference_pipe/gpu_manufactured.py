"""Original native backend/option checks, without the historical launcher."""
import re, shlex

def science_options(text):
    tokens = shlex.split(text)
    for flag in ('-mat_type', '-vec_type'):
        index = tokens.index(flag)
        del tokens[index:index + 2]
    return [value for value in tokens if value != '-log_view_gpu_time']

def backend_evidence(text):
    # Only actual PETSc object views qualify; command-line requests do not.
    types = set(re.findall(r'^\s*type:\s*(\w+)\s*$', text, re.M))
    matrix = sorted(types & {'aijcusparse', 'seqaijcusparse', 'mpiaijcusparse'})
    vector = sorted(types & {'cuda', 'seqcuda', 'mpicuda'})
    return dict(verified=bool(matrix and vector), matrix_types=matrix,
                vector_types=vector, device_backend='CUDA' if matrix and vector else None)

def fatal_reason(text):
    if re.search(r'DIVERGED_\w+', text):
        return 'GPU_LINEAR_DIVERGENCE'
    if re.search(r'\b(?:nan|[+-]?inf)\b', text, re.I):
        return 'NONFINITE_NATIVE_OUTPUT'
    if re.search(r'CUDA error|CUDA_ERROR|out of memory|CUDART error|PETSC ERROR|'
                 r'error while loading shared libraries|lib.*not found', text, re.I):
        return 'GPU_BACKEND_OR_DEVICE_ERROR'
    return None
