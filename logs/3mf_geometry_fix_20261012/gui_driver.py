"""Isolated Bambu GUI inspection; no printer connection or slicing actions."""
from pathlib import Path
import json, os, subprocess, sys, time
from PIL import ImageGrab
from Xlib import X, XK, display
from Xlib.ext import xtest
out=Path(__file__).resolve().parent/'gui';root=out.parents[2]
cmd=sys.argv[1]
if cmd=='start':
    xdg=out/'xdg';locale=out/'locale';locale.mkdir(exist_ok=True)
    for name in ['en_US','en_GB']:
        dest=locale/(name+'.UTF-8')
        if not dest.exists():subprocess.run(['localedef','-i',name,'-f','UTF-8',str(dest)],check=True)
    stop=out/'STOP'
    if stop.exists():stop.unlink()
    disp=next(':'+str(n) for n in range(94,120) if not Path('/tmp/.X11-unix/X'+str(n)).exists())
    x=subprocess.Popen(['Xvfb',disp,'-screen','0','1440x1000x24','-nolisten','tcp'],stdout=(out/'xvfb.log').open('w'),stderr=subprocess.STDOUT);time.sleep(1)
    env=dict(os.environ,DISPLAY=disp,XDG_CONFIG_HOME=str(xdg),LOCPATH=str(locale),LC_ALL='en_US.UTF-8',LIBGL_ALWAYS_SOFTWARE='1',WEBKIT_DISABLE_SANDBOX_THIS_IS_DANGEROUS='1',SSL_CERT_FILE='/etc/ssl/certs/ca-certificates.crt')
    args=[str(root/'tools/bambu-studio-2.7.1.57/squashfs-root/AppRun'),*sys.argv[2:]]
    g=subprocess.Popen(args,env=env,stdout=(out/'gui_stdout.log').open('w'),stderr=(out/'gui_stderr.log').open('w'))
    (out/'session.json').write_text(json.dumps(dict(display=disp,pid=g.pid,xvfb_pid=x.pid,command=args)))
    print('GUI_STARTED',g.pid,disp,flush=True)
    try:
        for _ in range(1800):
            if stop.exists() or g.poll() is not None:break
            time.sleep(1)
    finally:
        print('GUI_EXIT',g.poll(),flush=True)
        if g.poll() is None:g.terminate()
        x.terminate()
    sys.exit()
s=json.loads((out/'session.json').read_text());d=display.Display(s['display'])
def key(name,press=True):xtest.fake_input(d,X.KeyPress if press else X.KeyRelease,d.keysym_to_keycode(XK.string_to_keysym(name)));d.sync()
def chord(mod,name):key(mod);key(name);key(name,False);key(mod,False)
def text(t):
    from Xlib import Xatom
    from Xlib.protocol import event
    window=d.screen().root.create_window(0,0,1,1,0,d.screen().root_depth)
    selection=d.intern_atom('CLIPBOARD');utf8=d.intern_atom('UTF8_STRING');targets=d.intern_atom('TARGETS')
    window.set_selection_owner(selection,X.CurrentTime);d.sync();chord('Control_L','v')
    end=time.monotonic()+1.0
    while time.monotonic()<end:
        while d.pending_events():
            e=d.next_event()
            if e.type==X.SelectionRequest:
                prop=e.property or e.target
                if e.target==targets:e.requestor.change_property(prop,Xatom.ATOM,32,[utf8,Xatom.STRING])
                else:e.requestor.change_property(prop,e.target,8,t.encode())
                e.requestor.send_event(event.SelectionNotify(time=e.time,requestor=e.requestor,selection=e.selection,target=e.target,property=prop),propagate=False);d.flush()
        time.sleep(.01)
    window.destroy();d.sync()
if cmd=='click':
    x,y=map(int,sys.argv[2:4]);xtest.fake_input(d,X.MotionNotify,x=x,y=y);xtest.fake_input(d,X.ButtonPress,1);xtest.fake_input(d,X.ButtonRelease,1);d.sync();time.sleep(2)
elif cmd=='key':
    key(sys.argv[2]);key(sys.argv[2],False);time.sleep(2)
elif cmd=='open':
    d.set_screen_saver(0,0,X.DontPreferBlanking,X.DontAllowExposures);d.force_screen_saver(X.ScreenSaverReset)
    for x,y in [(150,115),(207,167)]:
        xtest.fake_input(d,X.MotionNotify,x=x,y=y);xtest.fake_input(d,X.ButtonPress,1);xtest.fake_input(d,X.ButtonRelease,1);d.sync();time.sleep(.5)
    time.sleep(1);chord('Control_L','l');chord('Control_L','a');text(sys.argv[2]);key('Return');key('Return',False);time.sleep(5)
elif cmd=='location':
    chord('Control_L','l');chord('Control_L','a');text(sys.argv[2]);key('Return');key('Return',False);time.sleep(4)
elif cmd=='text':text(sys.argv[2])
print('WINDOWS',[(w.id,w.get_wm_name(),w.get_geometry().width,w.get_geometry().height) for w in d.screen().root.query_tree().children if w.get_wm_name()])
ImageGrab.grab(xdisplay=s['display']).save(out/'current.png')
