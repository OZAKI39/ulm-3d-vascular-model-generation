/* Chronos RTSP latency fix, 2026-10-15.
 * Locally developed adapter using the public GStreamer 0.10 API.
 * Applies only to the existing H.264 network sender. The camera executable,
 * encoder, recording settings, LCD and HDMI sinks are not modified.
 * Revert by removing the service's LD_PRELOAD drop-in and restarting it.
 */
#define _GNU_SOURCE
#include <gst/gst.h>
#include <gst/base/gstbasesink.h>
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

GstElement *gst_element_factory_make(const gchar *factory, const gchar *name)
{
    GstElement *(*original)(const gchar *, const gchar *);
    GstElement *element;
    original = dlsym(RTLD_NEXT, "gst_element_factory_make");
    if (!original) {
        fprintf(stderr, "CHRONOS_RTSP_LOWLATENCY: cannot resolve original factory\n");
        return NULL;
    }
    element = original(factory, name);
    if (element && factory && name && !strcmp(factory, "multiudpsink") &&
        !strcmp(name, "net-sink") && GST_IS_BASE_SINK(element)) {
        gst_base_sink_set_sync(GST_BASE_SINK(element), FALSE);
        fprintf(stderr, "CHRONOS_RTSP_LOWLATENCY: net-sink sync=%d (send without clock wait)\n",
                gst_base_sink_get_sync(GST_BASE_SINK(element)));
    }
    return element;
}
