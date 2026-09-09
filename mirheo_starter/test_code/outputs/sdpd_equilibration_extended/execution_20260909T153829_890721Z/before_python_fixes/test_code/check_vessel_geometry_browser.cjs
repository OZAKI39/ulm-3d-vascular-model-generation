/* Optional offline Windows Chrome check through a temporary CDP pipe.
 * No server, npm packages, desktop session, CUDA or hardware GPU required.
 */
'use strict';
const {spawn} = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const {createHash} = require('node:crypto');
const {StringDecoder} = require('node:string_decoder');
const [htmlPath, outputDir] = process.argv.slice(2);
if (!htmlPath || !outputDir) throw Error('Usage: node check_vessel_geometry_browser.cjs <HTML> <new output directory>');
fs.mkdirSync(outputDir, {recursive:false});
const profile = path.join(outputDir, 'temporary_browser_profile');
const browser = 'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe';
const args = ['--headless=new', '--disable-gpu', '--use-angle=swiftshader',
  '--enable-unsafe-swiftshader', '--no-first-run', '--no-default-browser-check',
  '--disable-background-networking', '--remote-debugging-pipe', '--user-data-dir='+profile,
  'about:blank'];
const child = spawn(browser, args, {stdio:['ignore','ignore','pipe','pipe','pipe']});
let sequence=0, buffer='', sessionId;
const requests=new Map(), exceptions=[], remoteRequests=[], decoder=new StringDecoder('utf8');
child.stderr.on('data', chunk=>fs.appendFileSync(path.join(outputDir,'browser_stderr.log'),chunk));
child.on('error', error=>{for (const req of requests.values()) req.reject(error);});
child.stdio[4].on('data', chunk=>{
  buffer += decoder.write(chunk);
  let end;
  while ((end=buffer.indexOf('\0'))>=0) {
    const message=JSON.parse(buffer.slice(0,end)); buffer=buffer.slice(end+1);
    const req=requests.get(message.id);
    if (req) {clearTimeout(req.timer);requests.delete(message.id);message.error?req.reject(Error(JSON.stringify(message.error))):req.resolve(message.result);}
    if (message.method==='Runtime.exceptionThrown') exceptions.push(message.params);
    if (message.method==='Network.requestWillBeSent' && /^https?:/i.test(message.params.request.url)) remoteRequests.push(message.params.request.url);
  }
});
function call(method,params={},session=sessionId) {
  return new Promise((resolve,reject)=>{
    const id=++sequence;
    const timer=setTimeout(()=>{requests.delete(id);reject(Error('CDP timeout: '+method));},30000);
    requests.set(id,{resolve,reject,timer});
    child.stdio[3].write(JSON.stringify({id,method,params,...(session?{sessionId:session}:{})})+'\0');
  });
}
async function evaluate(expression) {
  const result=await call('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});
  if (result.exceptionDetails) throw Error(JSON.stringify(result.exceptionDetails));
  return result.result.value;
}
const wait = ms=>new Promise(resolve=>setTimeout(resolve,ms));
function check(value,message) {if(!value)throw Error(message);}
const report={status:'FAIL',html_path:htmlPath,html_sha256:createHash('sha256').update(fs.readFileSync(htmlPath)).digest('hex'),
  human_review:'PENDING',renderer:'Windows headless Chrome / SwiftShader software rendering',checks:{}};
const deadline=setTimeout(()=>{child.kill();process.exitCode=2;},150000);
(async()=>{
  try {
    report.browser=await call('Browser.getVersion',{},null);
    const {targetId}=await call('Target.createTarget',{url:'about:blank'},null);
    ({sessionId}=await call('Target.attachToTarget',{targetId,flatten:true},null));
    await call('Runtime.enable'); await call('Page.enable'); await call('Network.enable');
    await call('Network.setBlockedURLs',{urls:['http://*','https://*']});
    await call('Emulation.setDeviceMetricsOverride',{width:1600,height:1200,deviceScaleFactor:1,mobile:false});
    await call('Page.navigate',{url:pathToFileURL(htmlPath).href});
    let loaded=false;
    for(let i=0;i<100;i++) {
      loaded=await evaluate("Boolean(document.getElementById('render-status')?.textContent.includes('已加载'))");
      if(loaded)break;
      await wait(350);
    }
    check(loaded,'Plotly scene did not finish loading');
    report.mesh=await evaluate(`(()=>{const a=JSON.parse(document.getElementById('geometry-audit-data').textContent),v=document.getElementById('vessel-view');return {faces:v.data.filter(t=>t.type==='mesh3d').reduce((n,t)=>n+t.i.length,0),expected:a.summary.triangle_count,canvases:v.querySelectorAll('canvas').length,aspect:v.layout.scene.aspectmode,ports:Object.keys(a.controls.ports).length};})()`);
    check(report.mesh.faces===report.mesh.expected && report.mesh.canvases>0 && report.mesh.aspect==='data','Mesh, canvas or aspect mismatch');
    report.checks.complete_mesh_rendered=true;
    await wait(600);
    let shot=await call('Page.captureScreenshot',{format:'png'});
    fs.writeFileSync(path.join(outputDir,'geometry_overview.png'),Buffer.from(shot.data,'base64'),{flag:'wx'});
    report.checks.hide_wall=await evaluate(`(async()=>{const a=JSON.parse(document.getElementById('geometry-audit-data').textContent),e=document.getElementById('wall-visible');e.checked=false;e.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,200));return a.controls.walls.every(i=>document.getElementById('vessel-view').data[i].visible===false);})()`);
    report.checks.wall_opacity=await evaluate(`(async()=>{const a=JSON.parse(document.getElementById('geometry-audit-data').textContent),visible=document.getElementById('wall-visible'),e=document.getElementById('wall-opacity');visible.checked=true;visible.dispatchEvent(new Event('change'));e.value='0.7';e.dispatchEvent(new Event('input'));await new Promise(r=>setTimeout(r,200));return a.controls.walls.every(i=>document.getElementById('vessel-view').data[i].opacity===0.7);})()`);
    report.checks.port_selection=await evaluate(`(async()=>{const a=JSON.parse(document.getElementById('geometry-audit-data').textContent),e=document.getElementById('port-select'),id=Object.keys(a.controls.ports)[0];e.value=id;e.dispatchEvent(new Event('change'));await new Promise(r=>setTimeout(r,300));const v=document.getElementById('vessel-view');return Object.entries(a.controls.ports).every(([p,ids])=>ids.every(i=>v.data[i].opacity===(p===id?1:0.12))) && document.querySelector('tr[data-entity="'+id+'"]').classList.contains('selected');})()`);
    shot=await call('Page.captureScreenshot',{format:'png'});
    fs.writeFileSync(path.join(outputDir,'selected_port.png'),Buffer.from(shot.data,'base64'),{flag:'wx'});
    // Exercise real mouse dragging and wheel events, not just layout assignment.
    const box=await evaluate("(()=>{const b=document.getElementById('vessel-view').getBoundingClientRect();return {x:b.x+500,y:b.y+300};})()");
    const before=await evaluate("JSON.stringify(document.getElementById('vessel-view').layout.scene.camera)");
    await call('Input.dispatchMouseEvent',{type:'mousePressed',x:box.x,y:box.y,button:'left',clickCount:1});
    for(let i=1;i<=10;i++)await call('Input.dispatchMouseEvent',{type:'mouseMoved',x:box.x+i*12,y:box.y+i*4,button:'left',buttons:1});
    await call('Input.dispatchMouseEvent',{type:'mouseReleased',x:box.x+120,y:box.y+40,button:'left',clickCount:1});
    await wait(500);
    const dragged=await evaluate("JSON.stringify(document.getElementById('vessel-view').layout.scene.camera)");
    report.checks.mouse_rotation=before!==dragged;
    await call('Input.dispatchMouseEvent',{type:'mouseWheel',x:box.x,y:box.y,deltaY:-240,deltaX:0});
    await wait(500);
    const zoomed=await evaluate("JSON.stringify(document.getElementById('vessel-view').layout.scene.camera)");
    report.checks.mouse_wheel_zoom=zoomed!==dragged;
    report.checks.reset_view=await evaluate(`(async()=>{document.getElementById('reset-view').click();await new Promise(r=>setTimeout(r,200));const e=document.getElementById('vessel-view').layout.scene.camera.eye;return e.x===1.4&&e.y===1.5&&e.z===1.1;})()`);
    report.checks.offline_no_remote_requests=remoteRequests.length===0;
    report.checks.no_uncaught_javascript_errors=exceptions.length===0;
    report.checks.human_review_remains_pending=await evaluate("JSON.parse(document.getElementById('geometry-audit-data').textContent).states.human_review==='PENDING'");
    check(Object.values(report.checks).every(Boolean),'One or more browser interactions failed');
    report.status='PASS';
  } catch(error) {report.error=String(error);process.exitCode=2;}
  finally {
    report.exceptions=exceptions;report.remote_requests=remoteRequests;
    report.completed_utc=new Date().toISOString();
    fs.writeFileSync(path.join(outputDir,'browser_checks.json'),JSON.stringify(report,null,2)+'\n',{flag:'wx'});
    try {await call('Browser.close',{},null);} catch {}
    await Promise.race([new Promise(resolve=>child.once('exit',resolve)),wait(3000)]);
    if(child.exitCode===null)child.kill();
    for(const req of requests.values())clearTimeout(req.timer);
    clearTimeout(deadline);
    try {fs.rmSync(profile,{recursive:true,force:true,maxRetries:3,retryDelay:300});} catch {}
    console.log(JSON.stringify({status:report.status,checks:report.checks,error:report.error,report:path.join(outputDir,'browser_checks.json')}));
  }
})();
