"""Task-private compiled CUDA kernels; float64 and disabled FMA contraction."""
from functools import lru_cache
from pathlib import Path
import os


@lru_cache(1)
def extension():
    from torch.utils.cpp_extension import load
    root=Path(__file__).parent/'cuda'
    os.environ.setdefault('MAX_JOBS','4')
    os.environ.setdefault('TORCH_CUDA_ARCH_LIST','8.9')
    return load(name='brava_population_v1_cuda',sources=[str(root/'bindings.cpp'),str(root/'kernels.cu')],
        extra_cflags=['-O3','-ffp-contract=off'],
        extra_cuda_cflags=['-O3','--fmad=false','--prec-div=true','--prec-sqrt=true'],verbose=False)
