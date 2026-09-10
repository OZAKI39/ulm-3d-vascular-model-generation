"""单细胞原生剪切基准。帮助和预检不导入 Mirheo、不运行 GPU。"""
import argparse
import json
from .single_rbc_benchmark.physics import load_config

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',default='py_scripts/single_rbc_benchmark.yaml')
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--preflight-only',action='store_true',help='资源与已有证据核查，绝不启动 GPU')
    g.add_argument('--prepare-cpu',action='store_true',help='仅首次新增案例编译、原生几何/材料力探针和短 CPU 双向耦合预检')
    g.add_argument('--execute',action='store_true',help='必须已有本轮明确授权；冷启动、材料检查、重复和加严计算')
    a=p.parse_args();c=load_config(a.config)
    if a.execute:
        from .single_rbc_benchmark.execution import execute
        result=execute(c)
    else:
        from .single_rbc_benchmark.workflow import preflight
        result=preflight(c,prepare_cpu=a.prepare_cpu)
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
