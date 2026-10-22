#!/usr/bin/env python3
"""Install a pinned local Orca runtime without root or system-library changes."""
import hashlib,json,os,platform,shutil,subprocess,urllib.request
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
VERSION='2.4.2'
NAME=f'OrcaSlicer_Linux_AppImage_Ubuntu2404_V{VERSION}.AppImage'
URL=f'https://github.com/OrcaSlicer/OrcaSlicer/releases/download/v{VERSION}/{NAME}'
PACKAGES=['libwebkit2gtk-4.1-0','libjavascriptcoregtk-4.1-0','libgstreamer-plugins-base1.0-0','libgstreamer-gl1.0-0','libsoup-3.0-0','libharfbuzz-icu0','libmanette-0.2-0','libenchant-2-2','libhyphen0','libsecret-1-0','liborc-0.4-0t64']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    if platform.machine()!='x86_64' or not shutil.which('apt-get'):raise RuntimeError('This bootstrap targets Ubuntu 24.04 x86_64; provide an installed Orca CLI on another platform.')
    out=ROOT/'tools/external/OrcaSlicer';out.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(out).free<2_000_000_000:raise RuntimeError('INSUFFICIENT_RUNTIME_DISK_SPACE')
    archive=out/f'OrcaSlicer-{VERSION}.AppImage'
    request=urllib.request.Request(f'https://api.github.com/repos/OrcaSlicer/OrcaSlicer/releases/tags/v{VERSION}',headers={'User-Agent':'vascular-support-planner'})
    release=json.load(urllib.request.urlopen(request,timeout=45));asset=next(a for a in release['assets'] if a['name']==NAME)
    expected=asset.get('digest','').removeprefix('sha256:')
    if not archive.exists() or archive.stat().st_size!=asset['size']:
        subprocess.run(['curl','-L','--fail','--retry','3','-C','-',URL,'-o',str(archive)],check=True)
    actual=sha(archive)
    if expected and actual!=expected:raise RuntimeError('ORCA_RELEASE_SHA256_MISMATCH')
    archive.chmod(0o755)
    if not (out/'squashfs-root/bin/orca-slicer').exists():subprocess.run([str(archive),'--appimage-extract'],cwd=out,check=True,stdout=subprocess.DEVNULL)
    libs=out/'hostlibs';libs.mkdir(exist_ok=True)
    subprocess.run(['apt-get','download',*PACKAGES],cwd=libs,check=True)
    for package in libs.glob('*.deb'):subprocess.run(['dpkg-deb','-x',str(package),str(libs/'root')],check=True)
    launcher=out/'run-orca.sh';launcher.write_text('''#!/usr/bin/env bash
set -euo pipefail
ORCA_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export LD_LIBRARY_PATH="$ORCA_ROOT/hostlibs/root/usr/lib/x86_64-linux-gnu:$ORCA_ROOT/squashfs-root/lib/orca-runtime:$ORCA_ROOT/squashfs-root/bin${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export LC_ALL=C
exec "$ORCA_ROOT/squashfs-root/bin/orca-slicer" "$@"
''');launcher.chmod(0o755)
    probe=subprocess.run([str(launcher),'--help'],capture_output=True,text=True,check=True);(out/'cli_help.txt').write_text(probe.stdout+'\n'+probe.stderr)
    manifest=dict(version=VERSION,release_url=release['html_url'],archive_sha256=actual,official_digest=asset.get('digest'),binary_sha256=sha(out/'squashfs-root/bin/orca-slicer'),
        packages={p.name:sha(p) for p in libs.glob('*.deb')},local_only=True)
    (out/'runtime_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(launcher)

if __name__=='__main__':main()
