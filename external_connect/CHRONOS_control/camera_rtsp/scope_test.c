#include <gst/gst.h>
#include <gst/base/gstbasesink.h>

int main(int argc, char **argv) {
    GstElement *network, *other, *display;
    gst_init(&argc,&argv);
    network=gst_element_factory_make("multiudpsink","net-sink");
    other=gst_element_factory_make("multiudpsink","another-sink");
    display=gst_element_factory_make("fakesink","net-sink");
    if(!network||!other||!display) return 2;
    if(gst_base_sink_get_sync(GST_BASE_SINK(network))) return 3;
    if(!gst_base_sink_get_sync(GST_BASE_SINK(other))) return 4;
    /* fakesink has sync=false by default; ensure it still constructs normally. */
    g_print("PASS: only multiudpsink named net-sink disables sync; unrelated sink retains sync\n");
    gst_object_unref(network);gst_object_unref(other);gst_object_unref(display);
    return 0;
}
