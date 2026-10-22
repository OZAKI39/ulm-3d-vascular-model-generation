"""CPU-only presentation environment; never imports the GPU integration runner."""
import sys,time
from .common import BASE,write,identity

def main():
    stage=sys.argv[1];start=time.perf_counter()
    if stage=='figures':
        from .figures import main as figures
        figures()
    elif stage=='media':
        from .figures import main as figures
        from .videos import main as videos
        figures();videos()
    elif stage=='package':
        from .package import main as package
        package()
    elif stage=='report':
        from .report import main as report
        report()
    else:raise ValueError(stage)
    write(BASE/'evidence'/f'{stage}_TIMING.json',dict(wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
