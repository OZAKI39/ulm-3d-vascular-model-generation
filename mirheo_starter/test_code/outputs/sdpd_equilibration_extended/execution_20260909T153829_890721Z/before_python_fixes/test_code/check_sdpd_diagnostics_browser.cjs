/* Offline native Windows Chrome CDP review, adapted from the existing geometry
 * checker. No npm dependencies, web server or CUDA/GPU simulation invocation. */
'use strict';
const {spawn}=require('node:child_process');
const fs=require('node:fs'),path=require('node:path');
const {pathToFileURL}=require('node:url');
const {createHash}=require('node:crypto');
const {StringDecoder}=require('node:string_decoder');
const [htmlPath,out]=process.argv.slice(2);
const extended=process.env.MIRHEO_REVIEW_KIND==='extended';
const equilibration=extended||process.env.MIRHEO_REVIEW_KIND==='equilibration';
const readyExpression=equilibration?'window.equilibrationReady===true':'window.diagnosticsReady===true';
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
    for(let i=0;i<100;i++){if(await ev(readyExpression))break;await wait(200);}
    if(extended){
      report.checks.rendered=await ev("window.equilibrationReady===true && document.querySelectorAll('.js-plotly-plot').length===4");
      report.audit=await ev("JSON.parse(document.getElementById('equilibration-audit-data').textContent)");
      report.checks.fixed_extended_plan=report.audit.planned_steps===800000&&JSON.stringify(report.audit.formal_window_star)==='[0.6,0.8]';
      report.checks.seven_separate_observables=await ev("document.querySelectorAll('#outcomes tbody tr').length===7");
      report.checks.no_new_measurements_fabricated=report.audit.actual_steps===null&&report.audit.raw_csv_rows===0&&report.audit.new_temperature_mean===null;
      report.checks.restart_not_falsely_validated=report.audit.restart_validity==='RESTART_NOT_VALIDATED';
      report.checks.history_complete_and_labeled=await ev("(()=>{const t=document.getElementById('temperature')._fullData[0];return t.name.includes('历史参考')&&t.y.length===2001&&Math.max(...Array.from(t.y))>8})()");
      report.checks.late_fixed_500_samples=await ev("document.getElementById('late')._fullData[0].y.length===500");
      report.checks.pressure_500_actual_samples=await ev("document.getElementById('blocks')._fullData[0].y.length===500");
      report.checks.complete_budget_request=report.audit.requested_additional_s===1379&&report.audit.extra_authorized_gpu_seconds===0;
      report.checks.human_review_pending=report.audit.human_review==='PENDING'&&report.audit.selection===null;
    }else if(equilibration){
      report.checks.rendered=await ev("window.equilibrationReady===true && document.querySelectorAll('.js-plotly-plot').length===4");
      report.audit=await ev("JSON.parse(document.getElementById('equilibration-audit-data').textContent)");
      report.checks.fixed_plan=report.audit.planned_steps===400000&&JSON.stringify(report.audit.formal_window_star)==='[0.3,0.4]';
      report.checks.separate_outcomes=await ev("document.querySelectorAll('#outcomes tbody tr').length===4 && document.querySelectorAll('#budget-table tbody tr').length===6");
      report.checks.pending_not_fabricated=report.audit.actual_steps!==null||(report.audit.new_temperature_mean===null&&report.audit.stationary===null&&report.audit.temperature_match===null&&report.audit.raw_csv_rows===0);
      report.checks.history_labeled=await ev("document.getElementById('temperature').data.some(t=>t.name.includes('历史参考'))");
      report.checks.full_startup_peak_kept=await ev("document.getElementById('temperature')._fullData.some(t=>t.name.includes('历史参考') && Math.max(...Array.from(t.y))>8)");
      report.checks.actual_curve_complete=report.audit.actual_steps===null||await ev("(()=>{const a=JSON.parse(document.getElementById('equilibration-audit-data').textContent);const f=document.getElementById('temperature');const t=f._fullData.find(t=>t.name==='本次连续演化');return !!t && t.visible!==false && t.visible!=='legendonly' && t.y.length===a.raw_csv_rows && f.layout.xaxis.range[0]===0;})()");
      report.checks.formal_window_visible=report.audit.actual_steps!==400000||await ev("(()=>{const f=document.getElementById('late');const y=Array.from(f._fullData[0].y);const lo=Math.min(.98,...y),hi=Math.max(1.02,...y),pad=.25*(hi-lo),r=f._fullLayout.yaxis.range;return y.length===500 && r[0]>=lo-pad && r[1]<=hi+pad;})()");
      report.checks.human_review_pending=report.audit.human_review==='PENDING'&&report.audit.selection===null;
    }else{
    report.checks.rendered=await ev("window.diagnosticsReady===true && document.querySelectorAll('.js-plotly-plot').length===8");
    report.audit=await ev("JSON.parse(document.getElementById('diagnostics-audit-data').textContent)");
    report.checks.evidence_and_budget=report.audit.evidence_count>=12&&report.audit.shared_budget.total_charged_or_reserved_s<=3600;
    report.checks.tables=await ev("document.querySelectorAll('#sampling-table tbody tr').length>10 && document.querySelectorAll('#repair-table tbody tr').length===18");
    report.checks.human_review_pending=report.audit.human_review==='PENDING';
    report.checks.profile_dropdown=await ev(`(async()=>{const a=JSON.parse(document.getElementById('diagnostics-audit-data').textContent);const e=document.getElementById('profile-select');e.value=String(a.profile_options.length-1);e.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,200));const selected=a.profile_options.at(-1).indices;return document.getElementById('profile').data.every((t,i)=>t.visible===selected.includes(i));})()`);
    }
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
    for(let i=0;i<100;i++){if(await ev(readyExpression))break;await wait(200);}
    report.checks.default_view_restored=await ev(equilibration?readyExpression:"window.diagnosticsReady===true && document.getElementById('profile-caption').textContent.includes('sdpd_flow') && document.getElementById('temperature').data.every(t=>t.visible===(t.name.startsWith('sdpd_')?true:'legendonly'))");
    for(const [name,id] of (equilibration?[['summary',null],['temperature','temperature'],['late','late'],['blocks','blocks'],['evolution','evolution'],['budget','budget-table']]:[['summary',null],['profile','profile'],['temperature','temperature'],['eos','eos'],['probe','probe']])){
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
