"""Two-point scale calibration on an unresized Chronos stage-micrometer PNG."""
import argparse
import math
from pathlib import Path

from chronos_common import BASE, load_image, save_image, utc_now, write_json


def scale_from_points(points, distance_um):
    if not math.isfinite(distance_um) or distance_um <= 0:
        raise ValueError("Known distance must be positive and finite")
    x1, y1, x2, y2 = points
    if not all(math.isfinite(v) for v in points):
        raise ValueError("Coordinates must be finite")
    pixels = math.hypot(x2 - x1, y2 - y1)
    if pixels < 2:
        raise ValueError("Choose two distinct marks at least 2 pixels apart")
    return pixels, distance_um / pixels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    parser.add_argument("--distance-um", required=True, type=float)
    parser.add_argument("--points", nargs=4, type=float, metavar=("X1", "Y1", "X2", "Y2"))
    parser.add_argument("--note", default="", help="Objective, relay lens, zoom and camera settings")
    parser.add_argument("--output", default=str(BASE / "calibration.json"))
    args = parser.parse_args()
    import cv2
    image = load_image(args.image)
    h, w = image.shape[:2]
    points = args.points
    if points is None:
        clicked = []
        # Fit to a fixed display size, then map clicks to ORIGINAL image pixels.
        factor = min(1.0, 1200 / w, 850 / h)
        preview = cv2.resize(image, (round(w * factor), round(h * factor)))
        sx, sy = w / preview.shape[1], h / preview.shape[0]
        name = "Click two ruler marks | Enter: accept | R: reset | Esc: cancel"
        cv2.namedWindow(name, cv2.WINDOW_AUTOSIZE)

        def click(event, x, y, flags, param):
            if event == cv2.EVENT_LBUTTONDOWN and len(clicked) < 2:
                clicked.append((x * sx, y * sy))

        cv2.setMouseCallback(name, click)
        while True:
            display = preview.copy()
            for x, y in clicked:
                cv2.circle(display, (round(x / sx), round(y / sy)), 5, (0, 0, 255), -1)
            if len(clicked) == 2:
                a, b = clicked
                cv2.line(display, (round(a[0] / sx), round(a[1] / sy)),
                         (round(b[0] / sx), round(b[1] / sy)), (0, 255, 0), 1)
            cv2.imshow(name, display)
            key = cv2.waitKey(20) & 0xFF
            if key == 27 or cv2.getWindowProperty(name, cv2.WND_PROP_VISIBLE) < 1:
                cv2.destroyAllWindows()
                return 1
            if key in (ord("r"), ord("R")):
                clicked.clear()
            if key in (10, 13) and len(clicked) == 2:
                points = [v for p in clicked for v in p]
                break
        cv2.destroyAllWindows()
    if not (0 <= points[0] < w and 0 <= points[2] < w and 0 <= points[1] < h and 0 <= points[3] < h):
        parser.error("Calibration points must be inside the original image")
    pixels, scale = scale_from_points(points, args.distance_um)
    path = Path(args.output)
    result = {"created_utc": utc_now(), "source_image": str(Path(args.image).resolve()),
              "image_width": w, "image_height": h, "points_xy": points,
              "known_distance_um": args.distance_um, "distance_px": pixels,
              "um_per_pixel": scale, "note": args.note,
              "assumption": "Isotropic scale, unchanged optics, ROI, resolution and stream scaling."}
    write_json(path, result)
    a = (round(points[0]), round(points[1]))
    b = (round(points[2]), round(points[3]))
    cv2.line(image, a, b, (0, 0, 255), 2)
    save_image(path.with_suffix(".png"), image)
    print(f"Scale: {scale:.8g} um/pixel | {pixels:.3f} pixels = {args.distance_um:g} um\nSaved: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
