"""Display an optical clock, or summarize paired reference/stream clock readings."""
import argparse
import statistics
import time

from chronos_common import BASE, utc_now, write_json


def latency_summary(pairs):
    import math
    delays = [float(reference) - float(stream) for reference, stream in pairs]
    if not delays or any(not math.isfinite(d) or d < 0 for d in delays):
        raise ValueError("Use finite readings from the same clock session, reference >= streamed reading")
    return {"n": len(delays), "samples_ms": delays, "mean_ms": statistics.mean(delays),
            "median_ms": statistics.median(delays), "min_ms": min(delays), "max_ms": max(delays),
            "sample_std_ms": statistics.stdev(delays) if len(delays) > 1 else None}


def show_clock():
    import tkinter as tk
    root = tk.Tk()
    root.title("Chronos optical latency reference | elapsed milliseconds")
    root.geometry("850x280")
    label = tk.Label(root, font=("Consolas", 64, "bold"), bg="white", fg="black")
    label.pack(fill="both", expand=True)
    started = time.perf_counter()

    def update():
        label.config(text=f"{int((time.perf_counter()-started)*1000):09d} ms")
        root.after(8, update)

    update()
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clock", action="store_true")
    parser.add_argument("--pair", nargs=2, type=float, action="append", metavar=("REFERENCE_MS", "STREAM_MS"))
    parser.add_argument("--output", default=str(BASE / "data" / "latency.json"))
    args = parser.parse_args()
    if args.clock:
        show_clock()
        return 0
    if not args.pair:
        parser.error("Use --clock, or one/more --pair REFERENCE_MS STREAM_MS")
    result = latency_summary(args.pair)
    result.update(created_utc=utc_now(), paired_readings_ms=args.pair,
                  method="Clock directly visible and its camera preview captured together in the same photograph/screenshot.",
                  limitations="Approximate optical scene-to-display delay, affected by screen refresh, shutter and reading precision; not a synchronized hardware measurement.")
    write_json(args.output, result)
    print(f"Approximate optical delay: {result['mean_ms']:.1f} ms (n={result['n']})\nSaved: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
