/* Offline native Windows Chrome CDP review, adapted from the existing geometry
 * checker. No npm dependencies, web server or CUDA/GPU simulation invocation. */
'use strict';
const {spawn}=require('node:child_process');
const fs=require('node:fs'),path=require('node:path');
const {pathToFileURL}=require('node:url');
const {createHash}=require('node:crypto');
const {StringDecoder}=require('node:string_decoder');
const [htmlPath,out]=process.argv.slice(2);
if(!htmlPath||!out)throw Error('Usage: node check_sdpd_diagnostics_browser.cjs <HTML> <new output directory>');
fs.mkdirSync(out,{recursive:false});
const profile=path.join(out,'temporary_browser_profile');
const child=spawn('C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',[
  '--headless=new','--disable-gpu','--use-angle=swiftshader','--enable-unsafe-swiftshader',
  '--no-first-run','--no-default-browser-check','--disable-background-networking',
  '--remote-debugging-pipe','--user-data-dir='+profile,'about:blank'],{stdio:['ignore','ignore','pipe','pipe','pipe']});
let id=0,buffer='',sessionId;const requests=new Map(),exceptions=[],remote=[],decoder=new StringDecoder('utf8');
child.stderr.on('data',x=>fs.appendFileSync(path.join(out,'browser_stderr.log'),x));
child.stdio[4].on('data',x=>{buffer+=decoder.write(x);let end;while((end=buffer.indexOf('\0'))>=0){const m=JSON.parse(buffer.slice(0,end));buffer=buffer.slice(end+1);const q=requests.get(m.id);if(q){clearTimeout(q.timer);requests.delete(m.id);m.error?q.reject(Error(JSON.stringify(m.error))):q.resolve(m.result);}if(m.method==='Runtime.exceptionThrown')exceptions.push(m.params);if(m.method==='Network.requestWillBeSent'&&/^https?:/i.test(m.params.request.url))remote.push(m.params.request.url);}});
const call=(method,params={},session=sessionId)=>new Promise((resolve,reject)=>{const seq=++id;const timer=setTimeout(()=>{requests.delete(seq);reject(Error('CDP timeout '+method));},30000);requests.set(seq,{resolve,reject,timer});child.stdio[3].write(JSON.stringify({id:seq,method,params,...(session?{sessionId:session}:{})})+'\0');});
const wait=ms=>new Promise(r=>setTimeout(r,ms));
async function ev(expression){const r=await call('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
const report={html_path:htmlPath,html_sha256:createHash('sha256').update(fs.readFileSync(htmlPath)).digest('hex'),renderer:'Windows headless Chrome / SwiftShader software renderer',human_review:'PENDING',checks:{}};
const check=(condition,message)=>{if(!condition)throw Error(message);};
(async()=>{
  const timeout=setTimeout(()=>child.kill(),70000);
  try{
    report.browser=await call('Browser.getVersion',{},null);
    const {targetId}=await call('Target.createTarget',{url:'about:blank'},null);
    sessionId=(await call('Target.attachToTarget',{targetId,flatten:true},null)).sessionId;
    await call('Page.enable');await call('Runtime.enable');await call('Network.enable');
    await call('Emulation.setDeviceMetricsOverride',{width:1440,height:1050,deviceScaleFactor:1,mobile:false});
    await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
    for(let i=0;i<100;i++){if(await ev('window.diagnosticsReady===true'))break;await wait(200);}
    report.checks.rendered=await ev("window.diagnosticsReady===true && document.querySelectorAll('.js-plotly-plot').length===8");
    report.audit=await ev("JSON.parse(document.getElementById('diagnostics-audit-data').textContent)");
    report.checks.evidence_and_budget=report.audit.evidence_count>=12&&report.audit.shared_budget.total_charged_or_reserved_s<=3600;
    report.checks.tables=await ev("document.querySelectorAll('#sampling-table tbody tr').length>10 && document.querySelectorAll('#repair-table tbody tr').length===18");
    report.checks.human_review_pending=report.audit.human_review==='PENDING';
    report.checks.profile_dropdown=await ev(`(async()=>{const a=JSON.parse(document.getElementById('diagnostics-audit-data').textContent);const e=document.getElementById('profile-select');e.value=String(a.profile_options.length-1);e.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,200));const selected=a.profile_options.at(-1).indices;return document.getElementById('profile').data.every((t,i)=>t.visible===selected.includes(i));})()`);
    await ev("document.getElementById('temperature').scrollIntoView({block:'start'})");
    await wait(200);
    const before=await ev("JSON.stringify(document.getElementById('temperature').layout.xaxis.range)");
    const rect=await ev("(()=>{let r=document.getElementById('temperature').getBoundingClientRect();return {x:r.x+550,y:r.y+180}})()");
    await call('Input.dispatchMouseEvent',{type:'mouseMoved',x:rect.x,y:rect.y});
    await call('Input.dispatchMouseEvent',{type:'mouseWheel',x:rect.x,y:rect.y,deltaY:-160,deltaX:0});await wait(500);
    report.checks.mouse_wheel_zoom=before!==await ev("JSON.stringify(document.getElementById('temperature').layout.xaxis.range)");
    report.checks.reset=await ev("(async()=>{document.getElementById('reset-plots').click();await new Promise(r=>setTimeout(r,200));return document.getElementById('temperature').layout.xaxis.autorange===true})()");
    const previousLegend=await ev("document.getElementById('temperature').data[0].visible");
    const legendBox=await ev("(()=>{const b=document.querySelector('#temperature .legendtoggle').getBoundingClientRect();return {x:b.x+b.width/2,y:b.y+b.height/2};})()");
    await call('Input.dispatchMouseEvent',{type:'mousePressed',x:legendBox.x,y:legendBox.y,button:'left',clickCount:1});
    await call('Input.dispatchMouseEvent',{type:'mouseReleased',x:legendBox.x,y:legendBox.y,button:'left',clickCount:1});await wait(500);
    report.checks.legend=previousLegend!==await ev("document.getElementById('temperature').data[0].visible");
    // Capture the ordinary initial view after testing, without leftover zoom,
    // historical traces or the transient Plotly legend notification.
    await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
    for(let i=0;i<100;i++){if(await ev('window.diagnosticsReady===true'))break;await wait(200);}
    report.checks.default_view_restored=await ev("window.diagnosticsReady===true && document.getElementById('profile-caption').textContent.includes('sdpd_flow') && document.getElementById('temperature').data.every(t=>t.visible===(t.name.startsWith('sdpd_')?true:'legendonly'))");
    for(const [name,id] of [['summary',null],['profile','profile'],['temperature','temperature'],['eos','eos'],['probe','probe']]){
      await ev(id?`document.getElementById('${id}').scrollIntoView({block:'start'})`:'window.scrollTo(0,0)');await wait(300);
      const shot=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
      fs.writeFileSync(path.join(out,name+'.png'),Buffer.from(shot.data,'base64'),{flag:'wx'});
    }
    report.checks.offline_no_remote_requests=remote.length===0;report.checks.no_javascript_exceptions=exceptions.length===0;
    check(Object.values(report.checks).every(Boolean),'A rendering, interaction or offline check failed');report.status='PASS';
  }catch(error){report.status='FAIL';report.error=String(error);process.exitCode=1;}
  finally{
    report.exceptions=exceptions;report.remote_requests=remote;
    try{await call('Browser.close',{},null);}catch(_){}
    if(child.exitCode===null)await Promise.race([new Promise(r=>child.once('exit',r)),wait(3000)]);
    if(child.exitCode===null)child.kill();clearTimeout(timeout);
    try{fs.rmSync(profile,{recursive:true,force:true,maxRetries:5,retryDelay:100});report.profile_cleanup='PASS';}catch(error){report.profile_cleanup=String(error);}
    fs.writeFileSync(path.join(out,'browser_checks.json'),JSON.stringify(report,null,2),{flag:'wx'});
    console.log(JSON.stringify({status:report.status,checks:report.checks,output:out,error:report.error}));
  }
})();
