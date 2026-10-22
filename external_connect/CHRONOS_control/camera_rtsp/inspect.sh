printf '\n=== Versions ===\n'
dpkg-query -W chronos-video chronos-control chronos-gui
printf '\n=== Original program ===\n'
sha256sum /usr/bin/cam-pipeline
printf '\n=== Optimization configuration ===\n'
if test -f /etc/systemd/system/chronos-video.service.d/60-chronos-lowlatency.conf; then
    cat /etc/systemd/system/chronos-video.service.d/60-chronos-lowlatency.conf
else
    echo 'Persistent optimization not installed'
fi
printf '\n=== Services ===\n'
systemctl is-active chronos-video chronos-gui chronos-http
printf '\n=== Video state ===\n'
gdbus call --system --dest ca.krontech.chronos.video --object-path /ca/krontech/chronos/video --method ca.krontech.chronos.video.status
printf '\n=== Recent optimization messages ===\n'
journalctl -u chronos-video -n 100 --no-pager | grep CHRONOS_RTSP_LOWLATENCY || true
