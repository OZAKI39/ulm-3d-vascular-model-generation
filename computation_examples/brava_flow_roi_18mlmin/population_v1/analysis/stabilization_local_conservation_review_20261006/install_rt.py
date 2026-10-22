from common import *
import subprocess,importlib.metadata
subprocess.run([sys.executable,'-m','pip','install','--disable-pip-version-check','--no-deps','--target',str(HERE/'runtime/rt_dependencies'),'scikit-fem==12.0.2'],check=True)
import skfem
p=Path(skfem.__file__).parent
put('RT_TOOL_IDENTITY.json',dict(library='scikit-fem',version=skfem.__version__,install_scope='project runtime/rt_dependencies only; --no-deps',source='https://github.com/kinnala/scikit-fem',files={str(f.relative_to(p)):sha(f) for f in p.rglob('*.py')},numpy=np.__version__))
print('MATURE_RT_LIBRARY_READY',skfem.__version__)
