"""Distinguish the weak central Y extremum from the outer pressure maxima."""
import time
import numpy as np
from .common import BASE, context, read, write, identity


def main():
    begin=time.perf_counter();c=context();out=BASE/"diagnostics"
    original=read(out/"AXIS_PAIR_PRESSURE_MAXIMA.json");review={}
    for axis,d in original.items():
        ai="XYZ".index(axis);ids=[2,3] if axis=="X" else [0,1]
        rows={v["delta_plus_minus_deg"]:dict(v) for v in d["rows"]}
        for degree,row in rows.items():
            x=np.zeros((len(row["all_local_maxima_m"]),3));x[:,ai]=row["all_local_maxima_m"]
            p,_=c["discs"].evaluate(x)
            row["all_local_maxima_pressure_per_Pa"]=abs(np.array([1,np.exp(1j*np.deg2rad(degree))])@p[ids]).tolist()
        zero=rows[0];values=np.array(zero["all_local_maxima_pressure_per_Pa"])
        eligible=np.flatnonzero(values>=values.max()*(1-1e-10))
        initial=max([zero["all_local_maxima_m"][i] for i in eligible])
        previous=initial
        for degree in [0,*range(10,61,10),*range(-10,-61,-10)]:
            if degree==-10:previous=initial
            row=rows[degree];candidates=row["all_local_maxima_m"]
            selected=min(candidates,key=lambda q:abs(q-previous))
            assert degree==0 or abs(selected-previous)<2e-5, "Branch lost; never join unrelated maxima"
            row["selected_peak_m"]=selected;previous=selected
        series=[dict(delta_plus_minus_deg=k,**{n:v for n,v in row.items() if n!='delta_plus_minus_deg'}) for k,row in sorted(rows.items())]
        slope=np.polyfit([r["delta_plus_minus_deg"] for r in series],[r["selected_peak_m"]*1e6 for r in series],1)[0]
        review[axis]=dict(rows=series,branch_slope_um_per_degree=float(slope),initial_peak_m=initial,
                         selection="Strongest axial pressure maximum at delta=0, positive-position tie rule, then nearest continuous branch with <=20um displacement per 10 degrees",
                         scope="Field-only local pressure maximum; not calibrated trapping or a proven stable 3D force equilibrium")
    write(out/"AXIS_PAIR_PRESSURE_MAXIMA_REVIEWED.json",review)
    write(out/"PEAK_BRANCH_REVIEW.json",dict(identity=identity(),wall_s=time.perf_counter()-begin,
          superseded="AXIS_PAIR_PRESSURE_MAXIMA.json Y linear_fit_um_per_degree is INVALID: nearest-to-center selection switches between a tiny central extremum and outer maxima",
          correction="All original extrema retained; explicitly track the strong positive outer branch at phase zero. No pressure data altered.",
          results={a:v["branch_slope_um_per_degree"] for a,v in review.items()},Y_scope="HOMOGENEOUS_WATER_PROXY_SPECIFIC"))
    print({a:v["branch_slope_um_per_degree"] for a,v in review.items()})


if __name__=="__main__":main()
