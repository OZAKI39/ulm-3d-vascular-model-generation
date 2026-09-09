/* Offline native Windows Chrome CDP review, adapted from the existing geometry
 * checker. No npm dependencies, web server or CUDA/GPU simulation invocation. */
'use strict';
const {spawn}=require('node:child_process');
const fs=require('node:fs'),path=require('node:path');
const {pathToFileURL}=require('node:url');
const {createHash}=require('node:crypto');
const {StringDecoder}=require('node:string_decoder');
const [htmlPath,out]=process.argv.slice(2);
if(!htmlPath||!out)throw Error('Usage: node check_fluid_physics_browser.cjs <HTML> <new output directory>');
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
const timeout=setTimeout(()=>{child.kill();process.exitCode=2;},150000);
(async()=>{
  try{
    report.browser=await call('Browser.getVersion',{},null);const {targetId}=await call('Target.createTarget',{url:'about:blank'},null);({sessionId}=await call('Target.attachToTarget',{targetId,flatten:true},null));
    await call('Runtime.enable');await call('Page.enable');await call('Network.enable');await call('Network.setBlockedURLs',{urls:['http://*','https://*']});
    await call('Emulation.setDeviceMetricsOverride',{width:1600,height:1150,deviceScaleFactor:1,mobile:false});await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
    let loaded=false;for(let i=0;i<100;i++){loaded=await ev('window.physicsReady===true');if(loaded)break;await wait(350);}check(loaded,'Plotly did not finish');
    const audit=await ev(`(()=>{const a=JSON.parse(document.getElementById('physics-audit-data').textContent),v=document.getElementById('vessel-view');return {expected:a.geometry_triangle_count,faces:v.data.filter(t=>t.type==='mesh3d').reduce((s,t)=>s+t.i.length,0),canvases:v.querySelectorAll('canvas').length,aspect:v.layout.scene.aspectmode,task_ids:a.actual_task_ids,charts:Object.keys(a.charts).map(id=>({id,traces:document.getElementById(id).data.length})),state:a.summary_status};})()`);
    report.audit=audit;check(audit.faces===audit.expected&&audit.canvases>0&&audit.aspect==='data','Geometry rendering mismatch');check(audit.charts.every(x=>x.traces>0),'Missing measured chart');report.checks.complete_geometry_and_real_charts=true;
    async function screenshot(name,id){if(id)await ev(`document.getElementById(${JSON.stringify(id)}).scrollIntoView({block:'start'})`);await wait(500);const r=await call('Page.captureScreenshot',{format:'png'});fs.writeFileSync(path.join(out,name+'.png'),Buffer.from(r.data,'base64'),{flag:'wx'});}
    await screenshot('physics_overview');await screenshot('measured_profile','profile');await screenshot('pressure_response','eos');await screenshot('vessel_measurement_planes','vessel-view');
    report.checks.hide_wall=await ev(`(async()=>{const a=JSON.parse(document.getElementById('physics-audit-data').textContent),e=document.getElementById('wall-visible');e.checked=false;e.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,150));return a.controls.walls.every(i=>document.getElementById('vessel-view').data[i].visible===false);})()`);
    report.checks.port_selection=await ev(`(async()=>{const a=JSON.parse(document.getElementById('physics-audit-data').textContent),e=document.getElementById('port-select'),id=Object.keys(a.controls.ports)[0];e.value=id;e.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,150));return Object.entries(a.controls.ports).every(([p,is])=>is.every(i=>document.getElementById('vessel-view').data[i].opacity===(p===id?1:.12)));})()`);
    report.checks.wall_opacity=await ev(`(async()=>{const a=JSON.parse(document.getElementById('physics-audit-data').textContent),e=document.getElementById('wall-opacity');e.value='.7';e.dispatchEvent(new Event('input'));await new Promise(r=>setTimeout(r,150));return a.controls.walls.every(i=>document.getElementById('vessel-view').data[i].opacity===.7);})()`);
    await ev("document.getElementById('vessel-view').scrollIntoView({block:'start'})");
    const box=await ev("(()=>{const b=document.getElementById('vessel-view').getBoundingClientRect();return {x:b.x+500,y:b.y+300};})()");
    const before=await ev("JSON.stringify(document.getElementById('vessel-view').layout.scene.camera)");
    await call('Input.dispatchMouseEvent',{type:'mousePressed',x:box.x,y:box.y,button:'left',clickCount:1});for(let i=1;i<=8;i++)await call('Input.dispatchMouseEvent',{type:'mouseMoved',x:box.x+i*12,y:box.y+i*4,button:'left',buttons:1});await call('Input.dispatchMouseEvent',{type:'mouseReleased',x:box.x+96,y:box.y+32,button:'left',clickCount:1});await wait(250);
    const rotated=await ev("JSON.stringify(document.getElementById('vessel-view').layout.scene.camera)");report.checks.mouse_rotation=before!==rotated;
    await call('Input.dispatchMouseEvent',{type:'mouseMoved',x:box.x,y:box.y});await wait(300);
    report.wheel_debug_before=await ev(`({camera:document.getElementById('vessel-view').layout.scene.camera,scrollY:window.scrollY,target:document.elementFromPoint(${box.x},${box.y})?.tagName,events:window.vesselWheelEvents})`);
    await call('Input.dispatchMouseEvent',{type:'mouseWheel',x:box.x,y:box.y,deltaY:-240,deltaX:0});await wait(800);report.checks.wheel_zoom=rotated!==await ev("JSON.stringify(document.getElementById('vessel-view').layout.scene.camera)");
    report.wheel_debug_after=await ev("({camera:document.getElementById('vessel-view').layout.scene.camera,scrollY:window.scrollY,events:window.vesselWheelEvents})");
    report.checks.reset_view=await ev(`(async()=>{document.getElementById('reset-view').click();await new Promise(r=>setTimeout(r,150));const e=document.getElementById('vessel-view').layout.scene.camera.eye;return e.x===1.4&&e.y===1.5&&e.z===1.1;})()`);
    report.checks.offline_no_remote_requests=remote.length===0;report.checks.no_javascript_exceptions=exceptions.length===0;
    check(Object.values(report.checks).every(Boolean),'An interaction/offline check failed');report.status='PASS';
  }catch(error){report.status='FAIL';report.error=String(error);process.exitCode=1;}
  finally{
    report.exceptions=exceptions;report.remote_requests=remote;
    fs.writeFileSync(path.join(out,'browser_checks.json'),JSON.stringify(report,null,2),{flag:'wx'});
    try{await call('Browser.close',{},null);}catch(_){}
    await Promise.race([new Promise(r=>child.on('exit',r)),wait(3000)]);if(child.exitCode===null)child.kill();clearTimeout(timeout);
    console.log(JSON.stringify({status:report.status,checks:report.checks,output:out,error:report.error}));
  }
})();
