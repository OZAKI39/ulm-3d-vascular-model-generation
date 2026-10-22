"""Decode exported video bytes; retain actual first/middle/last frames."""
from pathlib import Path
import hashlib
import imageio_ffmpeg
from PIL import Image, ImageDraw, ImageFont
from .particle82_provenance import atomic_json,sha256


def decode_video(path, expected_frames=None, inspection=None):
    path=Path(path);reader=imageio_ffmpeg.read_frames(str(path),pix_fmt='rgb24');metadata=next(reader)
    selected={0,expected_frames//2,expected_frames-1} if expected_frames else {0}
    distinct=set();saved=[];count=0
    if inspection is not None:Path(inspection).mkdir(parents=True,exist_ok=True)
    for index,pixels in enumerate(reader):
        count=index+1;distinct.add(hashlib.sha256(pixels).hexdigest())
        if index in selected and inspection is not None:
            image=Image.frombytes('RGB',metadata['size'],pixels)
            target=Path(inspection)/f'{path.stem}_frame_{index:04d}.png';image.save(target);saved.append(str(target))
    if expected_frames is not None and count!=expected_frames:raise ValueError('Decoded frame count mismatch')
    return dict(file=path.name,sha256=sha256(path),decoded_frames=count,distinct_frames=len(distinct),
        dimensions=list(metadata['size']),fps=metadata['fps'],codec=metadata['codec'],inspection_frames=saved,
        all_frames_decoded=True)


def audit_directory(output):
    import json
    output=Path(output);records=[]
    for path in sorted((output/'animations').glob('*.mp4')):
        manifest=json.loads((output/'frames'/(path.stem+'.json')).read_text())
        record=decode_video(path,manifest['frame_count'],output/'inspection');records.append(record)
        if record['distinct_frames']<manifest['frame_count']//2:raise ValueError('Animation is unexpectedly static')
    if len(records)!=3:raise ValueError('Expected exactly three scientific animations')
    sheet=Image.new('RGB',(1800,1050),'#090f19');d=ImageDraw.Draw(sheet)
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',17)
    for row,record in enumerate(records):
        d.text((12,350*row+8),record['file'],font=font,fill='#e5edf7')
        for col,p in enumerate(record['inspection_frames']):
            im=Image.open(p).convert('RGB');im.thumbnail((596,310));sheet.paste(im,(600*col,350*row+35))
    sheet.save(output/'inspection'/'ANIMATION_STORYBOARD.png')
    atomic_json(output/'MEDIA_VALIDATION.json',dict(all_pass=True,records=records,
        audit='EVERY_EXPORTED_MP4_FRAME_DECODED; STORYBOARD_FROM_ACTUAL_VIDEO_BYTES'))
    return records
