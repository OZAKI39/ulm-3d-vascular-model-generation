"""生成/查看保存的单红细胞证据；不运行求解器。"""
import argparse
import subprocess
from py_scripts.single_rbc_benchmark.physics import load_config
from py_scripts.single_rbc_benchmark.reporting import export

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',default='py_scripts/single_rbc_benchmark.yaml');p.add_argument('--open',action='store_true');a=p.parse_args()
    result=export(load_config(a.config));print(result['html'])
    if a.open:
        target=subprocess.check_output(['wslpath','-w',result['html']],text=True).strip()
        subprocess.Popen(['/mnt/c/Windows/explorer.exe',target],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)

if __name__=='__main__':main()
