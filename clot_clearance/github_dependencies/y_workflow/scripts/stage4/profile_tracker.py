"""Short instrumented real-CFD case, separate from formal branch statistics."""
import cProfile
import io
import pstats
import argparse
from scripts.stage4.common import new_run, finish
from scripts.stage4.experiments import setup,run_case


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case",default="B0")
    parser.add_argument("--duration",type=float,default=.02)
    parser.add_argument("--secondary",action="store_true",help="Profile the existing secondary-ON model without changing its guards")
    args = parser.parse_args()
    provided = setup()
    kind = "performance/profile" if args.case == "B0" else "performance/profile_"+args.case
    if args.secondary:
        kind += "_secondary"
    out = new_run(kind,dict(duration_s=args.duration,case=args.case,secondary=args.secondary,
                  purpose="Measured runtime diagnosis"),["phase_screening"])
    profiler = cProfile.Profile()
    profiler.enable()
    row = run_case(out,args.case,provided[1]["particles"]["seed"],secondary=args.secondary,
                   duration=args.duration,provided=provided)
    profiler.disable()
    profiler.dump_stats(str(out / "tracker.prof"))
    text = io.StringIO()
    pstats.Stats(profiler,stream=text).sort_stats("cumulative").print_stats(35)
    (out / "PROFILE.txt").write_text(text.getvalue())
    print(text.getvalue(),flush=True)
    finish(out,dict(status="PASS",row=row,scope="Instrumented short runtime sample, not branch-control evidence"))


if __name__ == "__main__":
    main()
