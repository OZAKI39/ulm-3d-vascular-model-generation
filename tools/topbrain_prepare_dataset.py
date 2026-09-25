#!/usr/bin/env python3
"""Explicit provisioning of the fixed official TopBrain 2025 v2 archive.

Never invoked automatically by the visualization entry point.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE='TopBrain_Data_Release_Batches1n2_081425.zip'
MD5='b703ea31cd1f0e7115a5d3e6e61f59b3'
URL=f'https://zenodo.org/records/16878417/files/{ARCHIVE}?download=1'


def prepare(download=False):
    base=ROOT/'data/external/TopBrain'
    dest=base/'downloads'/ARCHIVE
    free=shutil.disk_usage(ROOT).free
    if free < 12*1024**3:
        raise RuntimeError(f'INSUFFICIENT_DISK_SPACE: free_bytes={free}; require at least 12 GiB')
    base.mkdir(parents=True,exist_ok=True)
    dest.parent.mkdir(exist_ok=True)
    report={'dataset':'TopBrain Challenge Data Release','zenodo_record':16878417,'doi':'10.5281/zenodo.16878417',
            'source_url':URL,'source_archive_expected_md5':MD5,'archive':str(dest),'free_bytes_before':free}
    if download:
        command=['wget','-c','--tries=5','--timeout=45','--progress=dot:giga',URL,'-O',str(dest)]
        result=subprocess.run(command) if shutil.which('wget') else None
        if result is None or result.returncode:
            import requests
            links=[URL]
            try:
                response=requests.get('https://zenodo.org/api/records/16878417',timeout=60)
                report['api_status']=response.status_code
                response.raise_for_status()
                file=next(x for x in response.json()['files'] if x['key']==ARCHIVE)
                links.extend(file['links'].values())
            except Exception as exc:
                report['api_error']=str(exc)
            succeeded=False
            for link in dict.fromkeys(links):
                for attempt in range(5):
                    try:
                        offset=dest.stat().st_size if dest.exists() else 0
                        with requests.get(link,headers={'Range':f'bytes={offset}-'} if offset else {},stream=True,timeout=(30,90)) as response:
                            report['last_http_status']=response.status_code
                            if response.status_code==416 and offset:
                                succeeded=True
                                break
                            response.raise_for_status()
                            if offset and response.status_code!=206:
                                raise RuntimeError('Server ignored resume Range; refusing silent restart')
                            with dest.open('ab' if offset else 'wb') as handle:
                                for chunk in response.iter_content(1024*1024):
                                    handle.write(chunk)
                        succeeded=True
                        break
                    except Exception as exc:
                        report['last_error']=str(exc)
                        print(f'Download retry {attempt+1}: {exc}',flush=True)
                        time.sleep(5)
                if succeeded:
                    break
            if not succeeded:
                report['downloaded_bytes']=dest.stat().st_size if dest.exists() else 0
                (base/'download_failure.json').write_text(json.dumps(report,indent=2))
                raise RuntimeError(f'DOWNLOAD_FAILED: {report}')
    if not dest.is_file():
        raise RuntimeError('ARCHIVE_NOT_FOUND: explicitly pass --download to provision the official dataset')
    with dest.open('rb') as stream:
        actual=hashlib.file_digest(stream,'md5').hexdigest()
    report.update(archive_md5=actual,archive_bytes=dest.stat().st_size)
    if actual!=MD5:
        corrupt=dest.with_name(dest.name+f'.{time.time_ns()}.corrupt')
        dest.rename(corrupt)
        report['corrupt_path']=str(corrupt)
        (base/'download_failure.json').write_text(json.dumps(report,indent=2))
        raise RuntimeError(f'MD5_MISMATCH: expected {MD5}, actual {actual}; not extracted. Saved {corrupt}. Rerun --download.')
    extraction=base/ARCHIVE.removesuffix('.zip')
    marker=extraction/'.verified_archive_md5'
    if not marker.exists():
        extraction.mkdir(exist_ok=True)
        with zipfile.ZipFile(dest) as archive:
            report['uncompressed_bytes']=sum(x.file_size for x in archive.infolist())
            if shutil.disk_usage(base).free < report['uncompressed_bytes']+2*1024**3:
                raise RuntimeError('INSUFFICIENT_DISK_SPACE for extracted contents')
            for info in archive.infolist():
                target=(extraction/info.filename).resolve()
                if not target.is_relative_to(extraction.resolve()) or (info.external_attr>>16)&0o170000==0o120000:
                    raise RuntimeError(f'UNSAFE_ARCHIVE_MEMBER: {info.filename}')
            archive.extractall(extraction)  # zipfile verifies per-member CRC during extraction.
        marker.write_text(MD5+'\n')
    directories={name:list(extraction.rglob(name)) for name in ['imagesTr_topbrain_mr','labelsTr_topbrain_mr','imagesTr_topbrain_ct','labelsTr_topbrain_ct','itksnap_labelmap_txt']}
    if any(len(paths)!=1 for paths in directories.values()):
        raise RuntimeError(f'DATA_STRUCTURE_AMBIGUOUS: {directories}')
    common=Path(os.path.commonpath([str(p[0].parent) for p in directories.values()]))
    mount=ROOT/'data/TopBrain'
    if mount.is_symlink():
        if mount.resolve()!=common.resolve():
            raise RuntimeError(f'Existing mount points elsewhere: {mount}')
    elif mount.exists():
        raise RuntimeError(f'Refusing to replace existing data path: {mount}')
    else:
        mount.symlink_to(common,target_is_directory=True)
    report.update(extraction=str(extraction),topbrain_root=str(common),mount=str(mount),
                  directories={k:str(v[0]) for k,v in directories.items()},
                  counts={k:len(list(v[0].glob('*.nii.gz'))) for k,v in directories.items() if k!='itksnap_labelmap_txt'},
                  licenses=[str(p) for p in extraction.rglob('*') if p.is_file() and 'license' in p.name.lower()])
    sys.path.insert(0,str(ROOT))
    from vascular_processing.topbrain_dataset import discover
    pairs=discover(common,'2025')
    report['paired_cases']=[{'case_id':p.case_id,'modality':p.modality,'image':str(p.image),'label':str(p.label)} for p in pairs]
    report['status']='PASS'
    (base/'dataset_provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    (ROOT/'outputs/topbrain_validation/dataset_provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download',action='store_true')
    prepare(parser.parse_args().download)
