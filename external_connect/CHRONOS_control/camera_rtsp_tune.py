"""Inspect, apply or restore the validated Chronos 1.4 network-sender fix.

apply/restore restart camera video services. They refuse busy cameras and
cameras containing unsaved frames. apply is pinned to the tested binary.
"""
import argparse
from datetime import datetime
import getpass
import json
from pathlib import Path
import shlex
import sys

BASE=Path(__file__).resolve().parent

def run_remote(client, command, timeout=45):
    _,stdout,stderr=client.exec_command(command,timeout=timeout)
    body=stdout.read().decode('utf-8',errors='replace')
    errors=stderr.read().decode('utf-8',errors='replace')
    code=stdout.channel.recv_exit_status()
    return {'exit_code':code,'stdout':body,'stderr':errors}

def put(client, local, remote):
    stdin,stdout,stderr=client.exec_command('cat > '+shlex.quote(remote),timeout=10)
    stdin.write(Path(local).read_bytes());stdin.flush();stdin.channel.shutdown_write()
    if stdout.channel.recv_exit_status():
        raise RuntimeError(stderr.read().decode('utf-8',errors='replace'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['inspect','apply','restore'])
    parser.add_argument('--factory-password',action='store_true',help='Use documented factory root password; otherwise prompt locally')
    args=parser.parse_args()
    import paramiko
    config=json.loads((BASE/'chronos_config.json').read_text(encoding='utf-8-sig'))
    client=paramiko.SSHClient()
    client.load_host_keys(str(BASE/'camera_rtsp'/'known_hosts'))
    client.set_missing_host_key_policy(paramiko.RejectPolicy())
    password='chronos' if args.factory_password else getpass.getpass('Camera root SSH password: ')
    try:
        client.connect(config['camera_ip'],username='root',password=password,
                       look_for_keys=False,allow_agent=False,timeout=5,auth_timeout=5,banner_timeout=5)
        del password
        if args.action=='apply':
            put(client,BASE/'camera_rtsp'/'chronos_rtsp_lowlatency.c','/tmp/chronos_rtsp_lowlatency.c')
            put(client,BASE/'camera_rtsp'/'scope_test.c','/tmp/chronos-lowlatency-test.c')
        name='install' if args.action=='apply' else args.action
        command=(BASE/'camera_rtsp'/f'{name}.sh').read_text(encoding='utf-8-sig')
        result=run_remote(client,command)
        folder=BASE/'data'/'camera_maintenance'
        folder.mkdir(parents=True,exist_ok=True)
        path=folder/f'{datetime.now():%Y%m%d_%H%M%S}_{args.action}.json'
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        print(result['stdout'])
        if result['stderr']:print(result['stderr'],file=sys.stderr)
        print(f'Log: {path}')
        return result['exit_code']
    finally:
        client.close()

if __name__=='__main__':sys.exit(main())
