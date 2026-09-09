/* Offline Windows Chrome CDP check. Uses existing Node/Chrome, no server/npm. */
'use strict';
const {spawn}=require('node:child_process');
const fs=require('node:fs'),path=require('node:path'),os=require('node:os');
const {pathToFileURL}=require('node:url');
const {createHash}=require('node:crypto');
const {StringDecoder}=require('node:string_decoder');
const [htmlPath,resultPath,out]=process.argv.slice(2);
if(!htmlPath||!resultPath||!out)throw Error('Usage: node check_solver_benchmark_browser.cjs HTML RESULTS_JSON NEW_OUTPUT_DIRECTORY');
fs.mkdirSync(out,{recursive:false});
// Chromium file locking on WSL UNC shares is unreliable. Use a new Windows-local
// temporary profile; the reviewed HTML and evidence still remain in the workspace.
const profileDir=fs.mkdtempSync(path.join(os.tmpdir(),'hemo-review-'));
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
const report={html_path:htmlPath,html_sha256:createHash('sha256').update(fs.readFileSync(htmlPath)).digest('hex'),human_review:'PENDING',checks:{}};
function check(ok,message){if(!ok)throw Error(message);}
(async()=>{const timeout=setTimeout(()=>child.kill(),60000);try{
  report.browser=await call('Browser.getVersion',{},null);
  const {targetId}=await call('Target.createTarget',{url:'about:blank'},null);
  sessionId=(await call('Target.attachToTarget',{targetId,flatten:true},null)).sessionId;
  await call('Page.enable');await call('Runtime.enable');await call('Network.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1440,height:1080,deviceScaleFactor:1,mobile:false});
  await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
  for(let i=0;i<100;i++){if(await ev('window.benchmarkReady===true'))break;await wait(150);}
  const embedded=await ev("JSON.parse(document.getElementById('benchmark-audit-data').textContent)");
  const external=JSON.parse(fs.readFileSync(resultPath,'utf8'));
  report.checks.embedded_equals_results=JSON.stringify(embedded)===JSON.stringify(external);
  report.checks.six_plots=await ev("window.benchmarkReady===true&&document.querySelectorAll('.js-plotly-plot').length===6");
  report.checks.three_conclusions=await ev("document.querySelectorAll('.lead p').length===3");
  report.checks.every_cost_row=await ev("document.querySelectorAll('#costs tbody tr').length")===embedded.comparison.results.length;
  report.checks.null_speedup_when_unqualified=embedded.comparison.qualified_speedups.every(x=>x.qualified_speedup===null)||embedded.comparison.status==='COMPLETE_LOCAL_SCREEN';
  const hasMir=embedded.comparison.results.some(r=>r.backend==='Mirheo'&&r.accuracy?.measured_profile_si);
  report.checks.no_fake_thermal_trace=await ev("document.getElementById('temperature').data.length")===(hasMir?1:0);
  report.checks.human_pending=embedded.human_review==='PENDING';
  report.checks.zero_external_requests=remote.length===0;
  report.checks.no_js_exceptions=exceptions.length===0;
  await ev("document.getElementById('profile').scrollIntoView({block:'start'})");await wait(200);
  await ev("(async()=>{await Plotly.relayout('profile',{'xaxis.range':[1,2]});return true})()");
  await ev("document.getElementById('reset-plots').click()");await wait(200);
  report.checks.reset=await ev("document.getElementById('profile').layout.xaxis.autorange===true");
  const shot=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'profile.png'),Buffer.from(shot.data,'base64'));
  await ev('window.scrollTo(0,0)');await wait(100);
  const top=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  fs.writeFileSync(path.join(out,'overview.png'),Buffer.from(top.data,'base64'));
  for(const [k,v]of Object.entries(report.checks))check(v,k);
  report.status='PASS';
}catch(e){report.status='FAIL';report.error=String(e);process.exitCode=1;}finally{
  report.exceptions=exceptions;report.remote_requests=remote;
  fs.writeFileSync(path.join(out,'browser_check.json'),JSON.stringify(report,null,2));
  try{await call('Browser.close',{},null);}catch(_){}child.kill();clearTimeout(timeout);
  await wait(250);try{fs.rmSync(profileDir,{recursive:true,force:true});}catch(_){}
  process.stdout.write(JSON.stringify({status:report.status,checks:report.checks,error:report.error})+'\n');
}})();
