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
(async()=>{const timeout=setTimeout(()=>child.kill(),45000);try{
  report.browser=await call('Browser.getVersion',{},null);
  const {targetId}=await call('Target.createTarget',{url:'about:blank'},null);
  sessionId=(await call('Target.attachToTarget',{targetId,flatten:true},null)).sessionId;
  await call('Page.enable');await call('Runtime.enable');await call('Network.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1440,height:1080,deviceScaleFactor:1,mobile:false});
  await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
  for(let i=0;i<100;i++){if(await ev("document.readyState==='complete' && Boolean(document.querySelector('.js-plotly-plot')?._fullData?.length)"))break;await wait(150);}
  const frames=JSON.parse(fs.readFileSync(resultPath)).frames.filter(f=>f.vertices);
  const first=frames[0].vertices.map(v=>v[0]),last=frames.at(-1).vertices.map(v=>v[0]);
  const numerical=JSON.parse(fs.readFileSync(path.join(path.dirname(resultPath),'numerical_checks.json')));
  report.checks.no_endpoint_fabrication=await ev("document.querySelector('body > p').innerText.includes("+JSON.stringify(numerical.CLOUD_RBC_RUN_COMPLETE)+")");
  report.checks.probe_scope_explicit=await ev("document.body.innerText.includes('NOT_VERIFIED')");
  report.checks.human_pending=await ev("document.body.innerText.includes('PENDING')");
  report.checks.actual_first_coordinates=await ev("JSON.stringify(Array.from(document.querySelector('.js-plotly-plot')._fullData[0].x))==="+JSON.stringify(JSON.stringify(first)));
  report.checks.real_frames_only=await ev("document.querySelector('.js-plotly-plot')._transitionData._frames.length==="+frames.length);
  const before=await ev("JSON.stringify(document.querySelector('.js-plotly-plot').layout.scene)");
  await ev("Plotly.animate(document.querySelector('.js-plotly-plot'),["+JSON.stringify(String(frames.length-1))+"],{mode:'immediate',frame:{duration:0,redraw:true},transition:{duration:0}})");
  await wait(300);
  report.checks.actual_last_coordinates=await ev("JSON.stringify(Array.from(document.querySelector('.js-plotly-plot')._fullData[0].x))==="+JSON.stringify(JSON.stringify(last)));
  report.checks.common_axes=await ev("JSON.stringify(document.querySelector('.js-plotly-plot').layout.scene)")===before;
  report.checks.no_external_requests=remote.length===0;report.checks.no_js_exceptions=exceptions.length===0;
  const shot=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});fs.writeFileSync(path.join(out,'review.png'),Buffer.from(shot.data,'base64'));
  for(const [name,value]of Object.entries(report.checks))check(value,name);
  report.status='PASS';
}catch(e){report.status='FAIL';report.error=String(e);process.exitCode=1;}finally{
  report.fixture_only=htmlPath.includes('tests');report.exceptions=exceptions;report.remote_requests=remote;
  fs.writeFileSync(path.join(out,'browser_check.json'),JSON.stringify(report,null,2));
  try{await call('Browser.close',{},null);}catch(_){}child.kill();clearTimeout(timeout);
  await wait(250);try{fs.rmSync(profileDir,{recursive:true,force:true});}catch(_){}
  process.stdout.write(JSON.stringify({status:report.status,checks:report.checks,error:report.error})+'\n');
}})();
