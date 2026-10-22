set -eu
python3 - <<'PY'
import dbus
b=dbus.SystemBus()
c=dbus.Interface(b.get_object('ca.krontech.chronos.control','/ca/krontech/chronos/control'),'ca.krontech.chronos.control')
v=dbus.Interface(b.get_object('ca.krontech.chronos.video','/ca/krontech/chronos/video'),'ca.krontech.chronos.video')
assert c.get(['state'])['state']=='idle', 'Camera is busy'
s=v.status()
assert s['totalFrames']==0 and not s['filesave'] and not s['liverecord'], 'Video data must be saved first'
PY
dropin=/etc/systemd/system/chronos-video.service.d/60-chronos-lowlatency.conf
if test -f "$dropin"; then
    grep -q '^# Managed by CHRONOS_control camera_rtsp_tune.py$' "$dropin" || { echo 'Unrecognized drop-in'; exit 3; }
    rm -f "$dropin"
fi
systemctl daemon-reload
systemctl restart chronos-video
for service in chronos-gui chronos-http; do
    if test -f "/opt/chronos-rtsp-lowlatency/$service.was-active"; then systemctl start "$service"; fi
done
echo 'Original video startup configuration restored; optimization files retained for audit'
sha256sum /usr/bin/cam-pipeline
