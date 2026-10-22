"""Standalone read-only media checks and contact sheets for human visual review."""
from pathlib import Path
import sys,json,hashlib
import fitz
import numpy as np
from PIL import Image,ImageDraw

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 b=Path(sys.argv[1]);out=b/'evidence/media_review';out.mkdir(parents=True,exist_ok=True);rows=[];thumbs=[]
 for p in sorted((b/'figures').glob('*.pdf')):
  with fitz.open(p) as d:
   assert len(d)==1;page=d[0];pix=page.get_pixmap(matrix=fitz.Matrix(1,1),alpha=False)
   im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples);assert np.asarray(im).std()>5
   blocks=page.get_text('blocks');outside=[list(block[:4]) for block in blocks if block[0]<-.5 or block[1]<-.5 or block[2]>page.rect.width+.5 or block[3]>page.rect.height+.5]
   assert not outside,(p.name,outside)
   png=p.with_suffix('.png');assert png.exists();rows.append(dict(pdf=p.name,pdf_sha256=sha(p),png_sha256=sha(png),page_count=1,text_blocks_within_page=True))
   canvas=Image.new('RGB',(800,620),'white');im.thumbnail((780,570));canvas.paste(im,((800-im.width)//2,35));ImageDraw.Draw(canvas).text((10,10),p.stem,fill='black');thumbs.append(canvas)
 for first in range(0,len(thumbs),4):
  sheet=Image.new('RGB',(1600,1240),'#dddddd')
  for k,im in enumerate(thumbs[first:first+4]):sheet.paste(im,(800*(k%2),620*(k//2)))
  sheet.save(out/f'figures_{first//4+1}.png')
 videos=json.loads((b/'videos/MANIFEST.json').read_text());assert len(videos)==3
 sheet=Image.new('RGB',(1800,900),'#080d16')
 for k,v in enumerate(videos):
  assert v['fully_decoded'] and v['frames']>1;assert sha(b/'videos'/(v['name']+'.mp4'))==v['sha256']
  frames=sorted((b/'videos').glob(v['name']+'_*.png'));assert len(frames)==3
  for j,p in enumerate(frames):
   im=Image.open(p).convert('RGB');im.thumbnail((600,300));sheet.paste(im,(600*j,300*k))
 sheet.save(out/'videos.png')
 (out/'AUTOMATED.json').write_text(json.dumps(dict(figures=rows,videos=videos,script_sha256=sha(Path(__file__)),manual_review_pending=True),indent=2)+'\n')
 print('Checked',len(rows),'PDF/PNG pairs and',len(videos),'decoded video manifests')
if __name__=='__main__':main()
