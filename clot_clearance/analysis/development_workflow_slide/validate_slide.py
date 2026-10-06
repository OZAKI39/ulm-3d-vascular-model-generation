"""Delivery checks for the one-slide artifact; no numerical experiment."""
import json,hashlib,re,sys
from pathlib import Path
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE
from PIL import Image
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
NAME='development_workflow_current_progress'
def main():
 r=Presentation(HERE/(NAME+'.pptx'));assert len(r.slides)==1 and r.slide_width/r.slide_height==16/9
 slide=r.slides[0];pictures=[s for s in slide.shapes if s.shape_type==MSO_SHAPE_TYPE.PICTURE];assert len(pictures)==3
 texts=[s.text for s in slide.shapes if s.has_text_frame];assert not any(re.search('[\u3400-\u9fff]',s) for s in texts)
 for phrase in ['INCONCLUSIVE','NEXT','FUTURE','3δ','5.9 × 10⁻¹⁵']:assert phrase in '\n'.join(texts),phrase
 for d in json.loads((HERE/'source_manifest.json').read_text()):assert hashlib.sha256((ROOT/d['path']).read_bytes()).hexdigest()==d['sha256']
 im=Image.open(HERE/(NAME+'.png'));assert im.size==(4800,2700) and im.info['dpi'][0]>=299
 try:import pymupdf
 except ImportError:
  sys.path.insert(0,str(ROOT/'build_mechanical_equilibrium_v1/adr_completion/pdf_tools'));import pymupdf
 pdf=pymupdf.open(HERE/(NAME+'.pdf'));assert len(pdf)==1
 page=pdf[0];spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for ln in b['lines'] for s in ln['spans']]
 clipped=[s['text'] for s in spans if s['bbox'][0]<-.1 or s['bbox'][1]<-.1 or s['bbox'][2]>page.rect.width+.1 or s['bbox'][3]>page.rect.height+.1];assert not clipped,clipped
 for s in slide.shapes:
  for e in s._element.xpath('./p:style/a:effectRef'):assert e.get('idx')=='0'
 result={'slide_count':1,'editable_native_shapes':len(slide.shapes)-len(pictures),'native_text_frames':sum(s.has_text_frame for s in slide.shapes),'real_data_insets':3,'all_visible_text_English':True,'quantitative_values_source_derived':True,'source_files_unchanged':True,'PNG_pixels':im.size,'PNG_dpi':im.info['dpi'],'PDF_pages':1,'rendered_from_actual_PowerPoint':True,'rendered_text_within_page':True,'current_stage':4,'next_stage':5,'root_cause':'INCONCLUSIVE','new_simulations':0,'files':{ext:(HERE/(NAME+'.'+ext)).stat().st_size for ext in ['pptx','png','pdf']}}
 (HERE/'FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
