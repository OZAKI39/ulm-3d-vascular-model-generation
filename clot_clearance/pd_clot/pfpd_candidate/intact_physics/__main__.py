import argparse,json,sys
from .api import run
p=argparse.ArgumentParser(description='Exploratory quasistatic intact physics interface; no automatic backend fallback')
p.add_argument('request_json');a=p.parse_args();r=run(json.load(open(a.request_json)))
print(json.dumps(r,ensure_ascii=False,allow_nan=False));raise SystemExit(0 if r['success'] else 2)
