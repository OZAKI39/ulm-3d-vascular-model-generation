set -eu
camera_dir=/opt/chronos-rtsp-lowlatency
dropin=/etc/systemd/system/chronos-video.service.d/60-chronos-lowlatency.conf
expected=b04f19cd1c49dc2bb730bd2cbb5becae3c94c4cf7bc70108968ad4e0099b3d44
actual=$(sha256sum /usr/bin/cam-pipeline | cut -d ' ' -f 1)
test "$actual" = "$expected" || { echo 'Video binary differs from the tested version; stop for review'; exit 2; }
python3 - <<'PY'
import dbus
b=dbus.SystemBus()
c=dbus.Interface(b.get_object('ca.krontech.chronos.control','/ca/krontech/chronos/control'),'ca.krontech.chronos.control')
v=dbus.Interface(b.get_object('ca.krontech.chronos.video','/ca/krontech/chronos/video'),'ca.krontech.chronos.video')
p=c.get(['cameraSerial','state'])
s=v.status()
assert str(p['cameraSerial'])=='01909', 'Not the validated camera'
assert str(p['state'])=='idle', 'Camera is busy'
assert s['totalFrames']==0 and not s['filesave'] and not s['liverecord'], 'Video data must be saved first'
PY
if test -e "$dropin"; then
    grep -q '^# Managed by CHRONOS_control camera_rtsp_tune.py$' "$dropin" || { echo 'Unrecognized existing drop-in'; exit 3; }
fi
mkdir -p "$camera_dir" /etc/systemd/system/chronos-video.service.d
cp /tmp/chronos_rtsp_lowlatency.c "$camera_dir/chronos_rtsp_lowlatency.c"
gcc -Wall -Wextra -Werror -fPIC -shared -O2 "$camera_dir/chronos_rtsp_lowlatency.c" -o "$camera_dir/lowlatency.so.new" $(pkg-config --cflags --libs gstreamer-0.10 gstreamer-base-0.10) -ldl
gcc /tmp/chronos-lowlatency-test.c -o "$camera_dir/scope-test" $(pkg-config --cflags --libs gstreamer-0.10 gstreamer-base-0.10)
LD_PRELOAD="$camera_dir/lowlatency.so.new" "$camera_dir/scope-test"
mv "$camera_dir/lowlatency.so.new" "$camera_dir/lowlatency.so"
chmod 755 "$camera_dir/lowlatency.so"
for service in chronos-gui chronos-http; do
    if systemctl is-active "$service" >/dev/null; then touch "$camera_dir/$service.was-active"; fi
done
printf '# Managed by CHRONOS_control camera_rtsp_tune.py\n[Service]\nEnvironment="LD_PRELOAD=/opt/chronos-rtsp-lowlatency/lowlatency.so"\n' > "$dropin.new"
mv "$dropin.new" "$dropin"
# Remove only the exact diagnostic files created during this task.
rm -f /run/chronos-rtsp-latency/deadline
rm -f /run/systemd/system/chronos-video.service.d/60-chronos-latency.conf
systemctl daemon-reload
if systemctl restart chronos-video; then
    for service in chronos-gui chronos-http; do
        if test -f "$camera_dir/$service.was-active"; then systemctl start "$service"; fi
    done
    if systemctl is-active chronos-video >/dev/null; then
        echo 'Applied: RTSP sender sync=false; original video executable unchanged'
        sha256sum /usr/bin/cam-pipeline "$camera_dir/lowlatency.so"
        exit 0
    fi
fi
rm -f "$dropin"
systemctl daemon-reload
systemctl restart chronos-video
for service in chronos-gui chronos-http; do
    if test -f "$camera_dir/$service.was-active"; then systemctl start "$service"; fi
done
echo 'Apply failed; original service configuration restored'
exit 1
