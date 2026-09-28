"""Display-only typography and units for the compact BraVa viewer.

World coordinates, vessel scalar arrays and their colour lookup table stay in um.
Only coordinate text and a separate scalar-bar lookup table are expressed in mm.
"""
from functools import lru_cache
from pathlib import Path
import subprocess

import numpy as np

from utils.rodent_vasculature import interactive as ui

ORIENTATION_LEGEND_LABEL = 'Centerline orientation'
CANDIDATE_LEGEND_LABEL = 'MeVO candidate ROI'
LEGEND_CENTRE_Y = .875
COMPACT_LEGEND_FONT_SIZE = 14
COMPACT_AXIS_TITLE_FONT_SIZE = 16
SCALE_TEXT_GAP_PX = 5.
SCALE_TITLE_FONT_SIZE = 16
SCALE_LABEL_FONT_SIZE = 14


@lru_cache(maxsize=1)
def display_font():
    """Prefer an installed Helvetica file; VTK's built-in Arial is the fallback."""
    try:
        listing = subprocess.run(['fc-list', '--format', '%{family}\t%{file}\n'],
            capture_output=True, text=True, check=True, timeout=5).stdout
    except (OSError, subprocess.SubprocessError):
        listing = ''
    fonts = []
    for line in listing.splitlines():
        if '\t' not in line:
            continue
        families, path = line.split('\t', 1)
        if Path(path).is_file():
            fonts.append(({name.strip().lower() for name in families.split(',')}, path))
    for family in ('Helvetica', 'Arial'):
        paths = [path for names, path in fonts if family.lower() in names]
        if paths:
            paths.sort(key=lambda path: (any(word in path.lower() for word in ('bold', 'italic', 'oblique')), path))
            return family, paths[0]
    return 'Arial', None


def apply_font(prop):
    _, path = display_font()
    if path:
        prop.SetFontFamilyAsString('File')
        prop.SetFontFile(path)
    else:
        prop.SetFontFamilyToArial()


def _scene_arrow_legend(actor, x, y, length_px, viewport_size):
    """Project the same vtkArrowSource used by PolyData.glyph's default arrow."""
    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkUnsignedCharArray
    from vtkmodules.vtkCommonDataModel import vtkCellArray
    from vtkmodules.vtkFiltersSources import vtkArrowSource

    source = vtkArrowSource()
    source.Update()
    mesh = pv.wrap(source.GetOutput()).copy()
    mesh.rotate_x(20., inplace=True)
    mesh = mesh.compute_normals(cell_normals=True, point_normals=False)
    # Painter's ordering and face lighting retain the cone/cylinder appearance
    # in this 2-D swatch, with the same orange base colour as the scene arrows.
    order = np.argsort(mesh.cell_centers().points[:, 2])
    cells = vtkCellArray()
    colours = vtkUnsignedCharArray()
    colours.SetNumberOfComponents(3)
    base_colour = np.array(pv.Color('#FF9F1C').int_rgb)
    light = np.array([-.25, .45, 1.]); light /= np.linalg.norm(light)
    normals = mesh.cell_data['Normals']
    for index in order:
        cells.InsertNextCell(mesh.GetCell(int(index)).GetPointIds())
        shade = .45 + .55 * max(0., float(np.dot(normals[index], light)))
        colours.InsertNextTuple3(*(base_colour * shade).astype(int))
    mesh.SetPolys(cells)
    mesh.GetCellData().SetScalars(colours)
    width, height = viewport_size
    mesh.points[:, 0] = x + mesh.points[:, 0] * length_px / width
    mesh.points[:, 1] = y + mesh.points[:, 1] * length_px / height
    mesh.points[:, 2] = 0.
    mapper = actor.GetMapper()
    mapper.SetInputData(mesh)
    mapper.SetScalarModeToUseCellData()
    mapper.SetColorModeToDirectScalars()
    mapper.ScalarVisibilityOn()


def top_row_legend(renderer, prefix):
    """Centre the complete legend row using actual font metrics, below the top."""
    from vtkmodules.vtkRenderingCore import vtkTextRenderer

    entries = sorted((int(name.removeprefix(prefix + '_text_')), actor)
        for name, actor in renderer.actors.items() if name.startswith(prefix + '_text_'))
    if not entries:
        return
    width, height = renderer.GetSize()
    if width <= 0 or height <= 0:
        return
    text_widths = []
    for _, actor in entries:
        if actor.GetInput() == 'parent -> current':
            actor.SetInput(ORIENTATION_LEGEND_LABEL)
        elif actor.GetInput() == 'candidate ROI':
            actor.SetInput(CANDIDATE_LEGEND_LABEL)
        apply_font(actor.GetTextProperty())
        actor.GetTextProperty().SetFontSize(COMPACT_LEGEND_FONT_SIZE)
        bbox = [0, 0, 0, 0]
        vtkTextRenderer.GetInstance().GetBoundingBox(actor.GetTextProperty(), actor.GetInput(),
            bbox, renderer.GetRenderWindow().GetDPI())
        text_widths.append(bbox[1] - bbox[0] + 1)
    swatch_width = ui.LEGEND_SWATCH_WIDTH_FRACTION * width
    text_gap = ui.LEGEND_TEXT_GAP_FRACTION * width
    gap = max(6., min(28., (width * .97 - sum(text_widths)
        - len(entries) * (swatch_width + text_gap)) / max(1, len(entries) - 1)))
    total_width = sum(text_widths) + len(entries) * (swatch_width + text_gap) + (len(entries) - 1) * gap
    cursor = (width - total_width) * .5
    for (index, actor), text_width in zip(entries, text_widths):
        position = ((cursor + swatch_width + text_gap) / width, LEGEND_CENTRE_Y)
        symbol_actor = renderer.actors[f'{prefix}_symbol_{index}']
        if actor.GetInput() == ORIENTATION_LEGEND_LABEL:
            _scene_arrow_legend(symbol_actor, cursor / width, LEGEND_CENTRE_Y, swatch_width, (width, height))
        else:
            symbol = symbol_actor.GetMapper().GetInput()
            delta = np.asarray(position) - actor.GetPosition()
            points = symbol.GetPoints()
            for point_index in range(points.GetNumberOfPoints()):
                x, y, z = points.GetPoint(point_index)
                points.SetPoint(point_index, x + delta[0], y + delta[1], z)
            points.Modified()
            symbol.Modified()
        actor.SetPosition(*position)
        actor.SetTextScaleModeToNone()
        cursor += swatch_width + text_gap + text_width + gap
    background = renderer.actors.get(f'{prefix}_background')
    if background is not None:
        background.SetVisibility(False)


def coordinate_labels_mm(controller):
    """Format physical coordinates without moving or rescaling the vessel data."""
    # Static exports and the interactive compact viewer share this controller.
    # Keep one set of X/Y/Z ticks visible throughout camera movement.
    controller.single_edge_per_axis = True
    for axis_index, annotations in enumerate((controller.x_annotations, controller.y_annotations, controller.z_annotations)):
        for annotation in annotations:
            values = np.linspace(annotation.point1_xyz[axis_index], annotation.point2_xyz[axis_index],
                                 len(annotation.tick_actors)) / 1000.
            for value, actor in zip(values, annotation.tick_actors):
                actor.SetInput(f'{value:.2f}')
                actor.SetTextScaleModeToNone()
                actor.GetTextProperty().SetFontSize(14)
                apply_font(actor.GetTextProperty())
            annotation.title_actor.SetInput(f'{annotation.axis_name} (mm)')
            annotation.title_actor.GetTextProperty().SetFontSize(COMPACT_AXIS_TITLE_FONT_SIZE)
            apply_font(annotation.title_actor.GetTextProperty())
            annotation.tick_label_offset_px = 20.
            annotation.title_offset_px = 32.
            annotation.vertical_axis_name = 'Y' if controller.y_vertical else 'Z'
            annotation.endpoint_label_inset_px = 28.
    controller.update()


def bottom_diameter_scale_mm(renderer):
    """Give the bar its own mm LUT; never change the vessel mapper's um LUT."""
    from vtkmodules.vtkRenderingCore import vtkTextActor

    for actor in list(renderer.actors.values()):
        if not actor.IsA('vtkScalarBarActor'):
            continue
        # select() builds a fresh um bar; this guard also makes restyling safe.
        if actor.GetTitle() == ui.DIAMETER_COLORBAR_TITLE:
            original = actor.GetLookupTable()
            lookup = original.NewInstance()
            lookup.DeepCopy(original)
            lookup.SetRange(*(value / 1000. for value in original.GetRange()))
            actor.SetLookupTable(lookup)
        # Explicit overlay text avoids VTK squeezing the title and tick labels
        # into the colour strip. Its centre is the orientation widget's origin
        # at 10% of subplot height; normalized positions also survive resizing.
        actor.SetTitle('')
        actor.DrawTickLabelsOff()
        actor.SetOrientationToHorizontal()
        actor.SetBarRatio(1.)
        actor.SetTitleRatio(0.)
        actor.SetPosition(.245, .0875)
        actor.SetWidth(.65); actor.SetHeight(.025)
        actor.SetLabelFormat('%.2f')
        for prop in (actor.GetTitleTextProperty(), actor.GetLabelTextProperty(), actor.GetAnnotationTextProperty()):
            apply_font(prop)
        labels = [('title', 'Diameter (mm)', .57, .141, SCALE_TITLE_FONT_SIZE)]
        for index, value in enumerate(np.linspace(*actor.GetLookupTable().GetRange(), 5)):
            labels.append((f'tick_{index}', f'{value:.2f}', .245 + index * .65 / 4,
                           .064, SCALE_LABEL_FONT_SIZE))
        for name, text, x, y, size in labels:
            label = vtkTextActor()
            label.SetInput(text)
            label.SetTextScaleModeToNone()
            label.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
            label.SetPosition(x, y)
            prop = label.GetTextProperty()
            apply_font(prop)
            prop.SetFontSize(size)
            prop.SetColor(*(244 / 255.,) * 3)
            prop.SetJustificationToCentered()
            prop.SetVerticalJustificationToCentered()
            renderer.add_actor(label, name=f'diameter_mm_{name}', reset_camera=False,
                               pickable=False, render=False)


class StableOrientationLabels:
    """Fixed-size, centre-anchored letters moving continuously with the triad.

    The original widget and its coloured 3-D axes are retained. vtkCaptionActor2D
    auto-fits and reanchors text as the axes rotate; replacing only those captions
    with projected text removes the resulting size/alignment jumps.
    """
    def __init__(self, widget):
        from vtkmodules.vtkRenderingCore import vtkTextActor

        self.widget = widget
        self.renderer = widget.GetRenderer()
        self.axes = widget.GetOrientationMarker()
        self.camera = self.renderer.GetActiveCamera()
        self.axes.AxisLabelsOff()
        self.actors = []
        for letter in 'XYZ':
            actor = vtkTextActor()
            actor.SetInput(f'{letter} (mm)')
            actor.SetTextScaleModeToNone()
            actor.GetPositionCoordinate().SetCoordinateSystemToDisplay()
            prop = actor.GetTextProperty()
            apply_font(prop)
            prop.SetFontSize(ui.ORIENTATION_AXIS_FONT_SIZE)
            prop.SetBold(True)
            prop.SetColor(1., 1., 1.)
            prop.SetJustificationToCentered()
            prop.SetVerticalJustificationToCentered()
            self.renderer.AddActor2D(actor)
            self.actors.append(actor)
        self.observers = [(self.camera, self.camera.AddObserver('ModifiedEvent', self.update)),
                          (self.renderer, self.renderer.AddObserver('StartEvent', self.update))]
        self.update()

    def update(self, *_):
        lengths = self.axes.GetTotalLength()
        for index, actor in enumerate(self.actors):
            endpoint = [0., 0., 0.]
            endpoint[index] = lengths[index] * 1.20
            self.renderer.SetWorldPoint(*endpoint, 1.)
            self.renderer.WorldToDisplay()
            x, y, _ = self.renderer.GetDisplayPoint()
            if np.isfinite(x) and np.isfinite(y):
                actor.SetPosition(float(x), float(y))

    def dispose(self):
        for source, observer in self.observers:
            source.RemoveObserver(observer)
        self.observers.clear()
        for actor in self.actors:
            self.renderer.RemoveActor2D(actor)


class DualPanelAnnotations:
    """Shared annotation layout; hosts supply plotter, controllers and label state."""
    def style_annotations(self,panel):
        renderer=self.plotter.renderers[panel]
        coordinate_labels_mm(self.controllers[panel])
        top_row_legend(renderer,('full_scene_legend','sampling_roi_legend')[panel])
        self.legend_layout_sizes[panel]=renderer.GetSize()
        if panel==1:bottom_diameter_scale_mm(renderer)
        if self.orientation_labels[panel] is None:
            self.orientation_labels[panel]=StableOrientationLabels(renderer.axes_widget)

    def align_bottom_scale(self,*_):
        """Maintain centred legends and tightly spaced, aligned scale annotations."""
        for panel,renderer in enumerate(self.plotter.renderers):
            size=renderer.GetSize()
            if size!=self.legend_layout_sizes[panel]:
                top_row_legend(renderer,('full_scene_legend','sampling_roi_legend')[panel])
                self.legend_layout_sizes[panel]=size
        labels=self.orientation_labels[1]
        if labels is None:return
        renderer=self.plotter.renderers[1]
        height=renderer.GetSize()[1]
        if height<=0:return
        labels.renderer.SetWorldPoint(0.,0.,0.,1.)
        labels.renderer.WorldToDisplay()
        target=labels.renderer.GetDisplayPoint()[1]-renderer.GetOrigin()[1]
        for actor in renderer.actors.values():
            if actor.IsA('vtkScalarBarActor'):
                rect=[0,0,0,0];actor.GetScalarBarRect(rect,renderer)
                if rect[3]<=0:continue
                delta=target-(rect[1]+rect[3]/2)
                if abs(delta)>.75:
                    x,y=actor.GetPosition();actor.SetPosition(x,y+delta/height)
                # Anchor text to the colour strip itself, with a five-pixel
                # visible gap instead of widely separated fixed viewport rows.
                half_height=rect[3]/2
                title=renderer.actors.get('diameter_mm_title')
                if title is not None:
                    title.GetTextProperty().SetVerticalJustificationToBottom()
                    title.SetPosition(title.GetPosition()[0],(target+half_height+SCALE_TEXT_GAP_PX)/height)
                for index in range(5):
                    label=renderer.actors.get(f'diameter_mm_tick_{index}')
                    if label is not None:
                        label.GetTextProperty().SetVerticalJustificationToTop()
                        label.SetPosition(label.GetPosition()[0],(target-half_height-SCALE_TEXT_GAP_PX)/height)

