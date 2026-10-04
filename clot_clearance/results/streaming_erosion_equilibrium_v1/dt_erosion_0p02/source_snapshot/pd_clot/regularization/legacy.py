"""Explicitly labeled compatibility entry; the frozen legacy runner is intact."""
import argparse,json,warnings
from pathlib import Path
from ..fragmentation import run as legacy_run
from ..output import write_json


def run(config,output):
    c=json.loads(json.dumps(config));forced=c['damage'].get('D_break',1.)<1.
    c['forced_failure_verification']=forced
    if forced:
        warnings.warn('D_break < 1: forced_failure_verification = true. Legacy verification only.',RuntimeWarning)
        c['verification_label']='FORCED_FAILURE_VERIFICATION: '+c.get('verification_label','legacy')
    result=legacy_run(c,output)
    result['forced_failure_verification']=forced
    write_json(Path(output)/'SUMMARY.json',result)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(json.loads(a.config.read_text()),a.output)


if __name__=='__main__':main()
