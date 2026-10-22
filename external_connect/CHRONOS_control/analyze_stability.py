"""Baseline contour measurements on saved Chronos images; no AI model."""
import argparse
import csv
import json
from pathlib import Path

from chronos_common import load_image, save_image, utc_now, write_json


def measure_frame(frame, roi=None, polarity="dark", threshold=None, min_area=20):
    import cv2
    import numpy as np
    h, w = frame.shape[:2]
    x, y, rw, rh = roi or (0, 0, w, h)
    if min(x, y) < 0 or min(rw, rh) <= 0 or x + rw > w or y + rh > h:
        raise ValueError("ROI must lie inside the image")
    gray = cv2.cvtColor(frame[y:y+rh, x:x+rw], cv2.COLOR_BGR2GRAY)
    kind = cv2.THRESH_BINARY_INV if polarity == "dark" else cv2.THRESH_BINARY
    if threshold is None:
        kind |= cv2.THRESH_OTSU
    value, mask = cv2.threshold(gray, threshold or 0, 255, kind)
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    valid = []
    for contour in contours:
        bx, by, bw, bh = cv2.boundingRect(contour)
        # Reject clipped objects and foreground reaching the ROI boundary.
        if bx == 0 or by == 0 or bx + bw >= rw or by + bh >= rh:
            continue
        if cv2.contourArea(contour) >= min_area:
            valid.append(contour)
    result = {"detected": False, "threshold_used": value, "candidate_count": len(valid)}
    overlay = frame.copy()
    cv2.rectangle(overlay, (x, y), (x+rw-1, y+rh-1), (255, 128, 0), 1)
    if not valid:
        return result, overlay
    contour = max(valid, key=cv2.contourArea)
    moments = cv2.moments(contour)
    if moments["m00"] <= 0:
        return result, overlay
    cx, cy = x + moments["m10"] / moments["m00"], y + moments["m01"] / moments["m00"]
    area = cv2.contourArea(contour)
    perimeter = cv2.arcLength(contour, True)
    sides = cv2.minAreaRect(contour)[1]
    result.update(detected=True, centroid_x_px=cx, centroid_y_px=cy, area_px2=area,
                  aspect_ratio=max(sides) / min(sides) if min(sides) > 0 else None,
                  circularity=4 * np.pi * area / perimeter ** 2 if perimeter else None)
    shifted = contour + np.array([[[x, y]]], dtype=contour.dtype)
    cv2.drawContours(overlay, [shifted], -1, (0, 255, 0), 1)
    cv2.drawMarker(overlay, (round(cx), round(cy)), (0, 0, 255), cv2.MARKER_CROSS, 12, 1)
    return result, overlay


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--images", required=True, help="Directory of raw snapshots")
    parser.add_argument("--roi", nargs=4, type=int, metavar=("X", "Y", "WIDTH", "HEIGHT"))
    parser.add_argument("--polarity", choices=["dark", "bright"], default="dark")
    parser.add_argument("--threshold", type=int, help="Fixed 0..255 threshold; omit to use Otsu per frame")
    parser.add_argument("--min-area", type=float, default=20)
    parser.add_argument("--calibration", help="Use this calibration JSON; otherwise report pixels only")
    parser.add_argument("--output", help="New output directory; default: beside input images")
    args = parser.parse_args()
    if args.threshold is not None and not 0 <= args.threshold <= 255:
        parser.error("threshold must be 0..255")
    if args.min_area <= 0:
        parser.error("min-area must be positive")
    import numpy as np
    images = sorted(p for p in Path(args.images).iterdir() if p.suffix.lower() in (".png", ".jpg", ".jpeg", ".tif", ".tiff"))
    if not images:
        parser.error("No images found")
    calibration = json.loads(Path(args.calibration).read_text(encoding="utf-8-sig")) if args.calibration else None
    scale = float(calibration["um_per_pixel"]) if calibration else None
    if scale is not None and (not np.isfinite(scale) or scale <= 0):
        parser.error("Invalid scale in calibration file")
    from datetime import datetime
    output = Path(args.output) if args.output else Path(args.images).parent / ("analysis_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    expected_shape = None
    for index, path in enumerate(images):
        frame = load_image(path)
        shape = [frame.shape[1], frame.shape[0]]
        if expected_shape is not None and shape != expected_shape:
            raise ValueError("Image resolution changed within this sequence")
        expected_shape = shape
        if calibration and shape != [calibration["image_width"], calibration["image_height"]]:
            raise ValueError("Image resolution differs from calibration; recalibrate first")
        values, overlay = measure_frame(frame, args.roi, args.polarity, args.threshold, args.min_area)
        row = {"filename": path.name, **values}
        metadata_path = path.with_suffix(".json")
        row["received_utc"] = json.loads(metadata_path.read_text(encoding="utf-8"))["received_utc"] if metadata_path.exists() else ""
        if values["detected"] and scale is not None:
            row.update(centroid_x_um=values["centroid_x_px"]*scale,
                       centroid_y_um=values["centroid_y_px"]*scale, area_um2=values["area_px2"]*scale**2)
        rows.append(row)
        if index < 10:
            save_image(output / "overlays" / (path.stem + ".png"), overlay)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (output / "measurements.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    summary = {"created_utc": utc_now(), "images": len(rows), "detected": sum(r["detected"] for r in rows),
               "parameters": vars(args), "statistics": {},
               "limitations": "Largest fully contained outer contour per frame; no identity tracking. Inspect overlays. Holes are included in outer area. No ground-truth accuracy claim."}
    for column in ("centroid_x_px", "centroid_y_px", "area_px2", "aspect_ratio", "circularity",
                   "centroid_x_um", "centroid_y_um", "area_um2"):
        data = [r[column] for r in rows if r.get(column) is not None]
        if data:
            summary["statistics"][column] = {"n": len(data), "mean": float(np.mean(data)),
                                            "sample_std": float(np.std(data, ddof=1)) if len(data)>1 else None}
    write_json(output / "summary.json", summary)
    print(f"Detected {summary['detected']}/{len(rows)} images. Review overlays before using measurements.\nOutput: {output}")
    return 0 if summary["detected"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
