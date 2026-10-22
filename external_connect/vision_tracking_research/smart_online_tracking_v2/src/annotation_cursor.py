"""Keep the workbench brush circular when napari falls back to a cross."""
from qtpy.QtCore import Qt


class CircularBrushCursor:
    def __init__(self,workbench):
        self.w=workbench
        viewer=workbench.viewer
        # Run after napari's native cursor update, without patching its package.
        for emitter in (viewer.cursor.events.style,viewer.cursor.events.size,
                        viewer.scene.camera.events.zoom,viewer.layers.selection.events.active):
            emitter.connect(self.update,position='last')
        self.update()

    def update(self,event=None):
        w=self.w;viewer=w.viewer
        if w.closing:return
        layer=viewer.layers.selection.active
        if layer not in w.brush_layers() or layer.mode not in ('paint','erase') or viewer.cursor.style!='circle':return
        actual=viewer.cursor.size*(viewer.scene.camera.zoom if viewer.cursor.scaled else 1)
        # A minimum display ring makes a subpixel brush locatable, but never
        # changes the data-space brush radius. Disclose this in the brush help.
        overlay=viewer.canvas.overlays._brush_circle
        overlay.size=max(8,round(actual));overlay.position_is_frozen=False
        overlay.visible=True
        viewer.window._qt_viewer.canvas.native.setCursor(Qt.CursorShape.BlankCursor)
        w.update_brush_help()
