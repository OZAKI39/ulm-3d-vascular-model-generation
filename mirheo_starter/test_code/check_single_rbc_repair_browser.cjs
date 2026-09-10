/* Offline Windows Chrome CDP check. Uses existing Node/Chrome, no server/npm. */
'use strict';
const {spawn}=require('node:child_process');
const fs=require('node:fs'),path=require('node:path'),os=require('node:os');
const {pathToFileURL}=require('node:url');
const {createHash}=require('node:crypto');
const {StringDecoder}=require('node:string_decoder');
const [htmlPath,resultPath,out]=process.argv.slice(2);
if(!htmlPath||!resultPath||!out)throw Error('Usage: node check_single_rbc_benchmark_browser.cjs HTML RESULTS_JSON NEW_OUTPUT_DIRECTORY');
fs.mkdirSync(out,{recursive:false});
// Chromium file locking on WSL UNC shares is unreliable. Use a new Windows-local
// temporary profile; the reviewed HTML and evidence still remain in the workspace.
const profileDir=fs.mkdtempSync(path.join(os.tmpdir(),'rbc-review-'));
const child=spawn('C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',[
  '--headless=new','--disable-gpu','--use-angle=swiftshader','--enable-unsafe-swiftshader',
  '--no-first-run','--no-default-browser-check','--disable-background-networking',
  '--remote-debugging-pipe','--user-data-dir='+profileDir,'about:blank'],
  {stdio:['ignore','ignore','pipe','pipe','pipe']});
let id=0,buffer='',sessionId;const requests=new Map(),exceptions=[],remote=[],decoder=new StringDecoder('utf8');
child.stderr.on('data',x=>fs.appendFileSync(path.join(out,'browser_stderr.log'),x));
child.stdio[4].on('data',x=>{buffer+=decoder.write(x);let end;while((end=buffer.indexOf('\0'))>=0){const m=JSON.parse(buffer.slice(0,end));buffer=buffer.slice(end+1);const q=requests.get(m.id);if(q){clearTimeout(q.timer);requests.delete(m.id);m.error?q.reject(Error(JSON.stringify(m.error))):q.resolve(m.result);}if(m.method==='Runtime.exceptionThrown')exceptions.push(m.params);if(m.method==='Network.requestWillBeSent'&&/^https?:/i.test(m.params.request.url))remote.push(m.params.request.url);}});
const call=(method,params={},session=sessionId)=>new Promise((resolve,reject)=>{const seq=++id;const timer=setTimeout(()=>{requests.delete(seq);reject(Error('CDP timeout '+method));},20000);requests.set(seq,{resolve,reject,timer});child.stdio[3].write(JSON.stringify({id:seq,method,params,...(session?{sessionId:session}:{})})+'\0');});
const wait=ms=>new Promise(r=>setTimeout(r,ms));
async function ev(expression){const r=await call('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
const report={recorded_at:new Date().toISOString(),html_path:htmlPath,html_sha256:createHash('sha256').update(fs.readFileSync(htmlPath)).digest('hex'),checker_sha256:createHash('sha256').update(fs.readFileSync(__filename)).digest('hex'),human_review:'PENDING',checks:{}};
function check(ok,message){if(!ok)throw Error(message);}
(async()=>{const timeout=setTimeout(()=>child.kill(),60000);try{
  report.browser=await call('Browser.getVersion',{},null);
  const {targetId}=await call('Target.createTarget',{url:'about:blank'},null);
  sessionId=(await call('Target.attachToTarget',{targetId,flatten:true},null)).sessionId;
  await call('Page.enable');await call('Runtime.enable');await call('Network.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1440,height:1080,deviceScaleFactor:1,mobile:false});
  await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
  for(let i=0;i<100;i++){if(await ev('window.rbcReady===true'))break;await wait(150);}
  report.checks.repair_pending_is_explicit=await ev("JSON.parse(document.getElementById('repair-results').textContent).qualified_speedup===null&&document.getElementById('repair-current').textContent.includes('尚未完成同等质量比较')");
  report.checks.legacy_is_labelled=await ev("document.getElementById('repair-legacy').textContent.includes('旧 campaign 的原始离线结果')");
  const repairShot=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'repair-overview.png'),Buffer.from(repairShot.data,'base64'));
  const embedded=await ev("JSON.parse(document.getElementById('rbc-audit-data').textContent)");
  const external=JSON.parse(fs.readFileSync(resultPath,'utf8'));
  report.checks.embedded_equals_results=JSON.stringify(embedded)===JSON.stringify(external);
  report.checks.fifteen_plots=await ev("window.rbcReady===true&&document.querySelectorAll('.js-plotly-plot').length===15");
  report.checks.three_conclusions=await ev("document.querySelectorAll('.lead p').length===3");
  report.checks.every_cost_row=await ev("document.querySelectorAll('#runs-table tbody tr').length")===embedded.runs.length;
  report.checks.null_speedup_when_unqualified=embedded.qualified_speedup===null||(embedded.workflow==='PASS'&&embedded.model_comparability==='PASS'&&embedded.benchmark_screen==='PASS');
  const hasMir=embedded.runs.some(r=>r.task==='main_mirheo_1'&&r.frames.length);
  report.checks.no_fake_mirheo_mesh=await ev("document.getElementById('mirheo3d').data.some(t=>t.name==='实际膜网格')")===hasMir;
  report.checks.chinese_plot_titles=await ev("document.getElementById('deformation').layout.title.text==='统一形变 D'");
  report.checks.human_pending=embedded.human_review==='PENDING';
  report.checks.zero_external_requests=remote.length===0;
  report.checks.no_js_exceptions=exceptions.length===0;
  await ev("document.getElementById('strain-slider').value=document.getElementById('strain-slider').max;document.getElementById('strain-slider').dispatchEvent(new Event('input'))");await wait(400);
  report.checks.slider=await ev("Number(document.getElementById('strain-label').textContent)>0");
  report.checks.common_axis_scales=await ev("['hemo3d','mirheo3d'].every(id=>['xaxis','yaxis','zaxis'].every(a=>JSON.stringify(document.getElementById(id).layout.scene[a].range)==='[0,24]'))");
  report.checks.missing_endpoint_stays_blank=await ev("(()=>{const r=window.rbcData.runs.find(r=>r.task==='main_mirheo_1');const has=r?.frames.some(f=>f.phase==='shear'&&Math.abs(f.strain-4)<1e-7)||false;return document.getElementById('mirheo3d').data.some(t=>t.name==='实际膜网格')===has;})()");
  report.checks.actual_feedback_panels=await ev("['feedback-hemo','feedback-mirheo'].every(id=>document.getElementById(id).classList.contains('js-plotly-plot'))");
  await ev("document.getElementById('play-button').click()");await wait(1200);
  report.checks.play_advances_real_strain=await ev("document.getElementById('play-button').textContent==='暂停'&&Number(document.getElementById('strain-slider').value)>0&&Number(document.getElementById('strain-slider').value)<Number(document.getElementById('strain-slider').max)");
  await ev("document.getElementById('play-button').click()");
  report.checks.pause_control=await ev("document.getElementById('play-button').textContent==='播放'");
  await ev("document.getElementById('repeat-select').value='2';document.getElementById('repeat-select').dispatchEvent(new Event('change'))");await wait(400);
  report.checks.second_repeat_uses_saved_vertices=await ev("(()=>{const r=window.rbcData.runs.find(r=>r.task==='main_mirheo_2');const f=r.frames.find(f=>f.phase==='shear'&&Math.abs(f.strain)<1e-7);const t=document.getElementById('mirheo3d').data.find(t=>t.name==='实际膜网格');return JSON.stringify(t.x)===JSON.stringify(f.vertices.map(v=>v[0]));})()");
  await ev("document.getElementById('repeat-select').value='diagnostic';document.getElementById('repeat-select').dispatchEvent(new Event('change'))");await wait(400);
  await ev("document.getElementById('strain-slider').value='34';document.getElementById('strain-slider').dispatchEvent(new Event('input'))");await wait(400);
  report.checks.diagnostic_saved_frame=await ev("(()=>{const r=window.rbcData.runs.find(r=>r.task==='diagnostic_mirheo_full');const f=r.frames.find(f=>f.phase==='shear'&&Math.abs(f.strain-3.4)<1e-7);const t=document.getElementById('mirheo3d').data.find(t=>t.name==='实际膜网格');return Number(document.getElementById('strain-label').textContent)===3.4&&JSON.stringify(t.x)===JSON.stringify(f.vertices.map(v=>v[0]))&&document.getElementById('animation-note').textContent.includes('不同执行配置');})()");
  report.checks.legends_do_not_overlap_axis_titles=await ev("[...document.querySelectorAll('.js-plotly-plot')].every(p=>{const l=p.querySelector('g.legend'),t=p.querySelector('g.g-xtitle');if(!l||!t)return true;const a=l.getBoundingClientRect(),b=t.getBoundingClientRect();return !a.height||!b.height||a.top>=b.bottom+2||a.bottom<=b.top-2;})");
  await ev("document.getElementById('hemo3d').scrollIntoView({block:'start'})");await wait(200);
  const shot=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'cells.png'),Buffer.from(shot.data,'base64'));
  await ev("document.getElementById('deformation').scrollIntoView({block:'start'})");await wait(100);
  const dynamics=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'dynamics.png'),Buffer.from(dynamics.data,'base64'));
  await ev("document.getElementById('cost').scrollIntoView({block:'start'})");await wait(100);
  const costs=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'costs.png'),Buffer.from(costs.data,'base64'));
  await ev('window.scrollTo(0,0)');await wait(100);
  const top=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'overview.png'),Buffer.from(top.data,'base64'));
  report.checks.zero_external_requests=remote.length===0;
  report.checks.no_js_exceptions=exceptions.length===0;
  for(const [k,v]of Object.entries(report.checks))check(v,k);
  report.status='PASS';
}catch(e){report.status='FAIL';report.error=String(e);process.exitCode=1;}finally{
  report.exceptions=exceptions;report.remote_requests=remote;
  fs.writeFileSync(path.join(out,'browser_check.json'),JSON.stringify(report,null,2));
  try{await call('Browser.close',{},null);}catch(_){}child.kill();clearTimeout(timeout);
  await wait(250);try{fs.rmSync(profileDir,{recursive:true,force:true});}catch(_){}
  process.stdout.write(JSON.stringify({status:report.status,checks:report.checks,error:report.error})+'\n');
}})();
