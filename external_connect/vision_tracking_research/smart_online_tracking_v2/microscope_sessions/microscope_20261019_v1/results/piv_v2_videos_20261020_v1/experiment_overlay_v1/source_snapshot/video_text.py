"""Actual Arial glyphs, blended only over their small image region."""
from functools import lru_cache
from pathlib import Path
import hashlib
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont


@lru_cache(maxsize=1)
def arial_path():
    candidates = [os.environ.get('F10Y_ARIAL_FONT'), '/mnt/c/Windows/Fonts/arial.ttf',
                  '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf',
                  '/usr/share/fonts/truetype/msttcorefonts/arial.ttf',
                  str(Path.home()/'.local/share/fonts/arial.ttf')]
    for name in candidates:
        if name and Path(name).is_file() and ImageFont.truetype(name,16).getname()[0]=='Arial':
            return str(Path(name).resolve())
    raise RuntimeError('Arial font not found. Set F10Y_ARIAL_FONT to the installed arial.ttf.')


@lru_cache(maxsize=32)
def font(size):
    return ImageFont.truetype(arial_path(),max(6,round(size)))


@lru_cache(maxsize=256)
def glyph(text,size):
    face=font(size)
    x0,y0,x1,y1=face.getbbox(text,anchor='ls')
    mask=Image.new('L',(max(1,x1-x0),max(1,y1-y0)))
    ImageDraw.Draw(mask).text((-x0,-y0),text,font=face,fill=255,anchor='ls')
    return np.asarray(mask,dtype=np.uint16),x0,y0


def text_width(text,size):
    return font(size).getlength(text)


def draw_text(image,text,x,y,size=16,color=(232,226,216)):
    mask,dx,dy=glyph(str(text),round(size))
    x,y=round(x)+dx,round(y)+dy
    h,w=mask.shape
    left,top=max(0,x),max(0,y);right,bottom=min(image.shape[1],x+w),min(image.shape[0],y+h)
    if right<=left or bottom<=top:return
    alpha=mask[top-y:bottom-y,left-x:right-x,None]
    roi=image[top:bottom,left:right]
    roi[:]=((roi.astype(np.uint16)*(255-alpha)+np.asarray(color,dtype=np.uint16)*alpha+127)//255).astype(np.uint8)


@lru_cache(maxsize=1)
def font_record():
    path=arial_path()
    return dict(family='Arial',path=path,sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                captions='Essential measurements and status only; interpretation stored in metadata')


def register_qt_arial():
    from PySide6.QtGui import QFontDatabase
    if 'Arial' not in QFontDatabase.families():
        if QFontDatabase.addApplicationFont(arial_path())<0:
            raise RuntimeError('Could not load Arial for the observation interface.')
