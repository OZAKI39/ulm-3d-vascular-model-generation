"""Prepare a new case; never modify the accepted frozen_reference."""
from pathlib import Path
import hashlib, json, shutil, xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / 'flow_cases/mean-2p0-mmps'
Q = 1.551359160885543e-14

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    baseline = ROOT / 'frozen_reference'
    assert not (CASE/'run').exists(), 'Refuse to overwrite an existing run'
    (CASE/'run').mkdir(parents=True)
    (CASE/'reports').mkdir()
    shutil.copytree(baseline/'SV_MESH', CASE/'SV_MESH')
    tree = ET.parse(baseline/'run/solver.xml')
    tree.find(".//Add_BC[@name='INLET']/Value").text = '-1.551359160885543e-14'
    tree.write(CASE/'run/solver.xml', encoding='utf-8', xml_declaration=True)
    shutil.copyfile(baseline/'run/PETSC_OPTIONS.txt', CASE/'run/PETSC_OPTIONS.txt')
    old = json.loads((ROOT/'configs/sv1_3q/policy.json').read_text())['production_policy']
    policy = {k: old[k] for k in ('dt_s','h10_m','Dh_m','nu_m2_s','save_interval_steps',
        'steady_last_intervals','velocity_change_limit','flow_change_limit','mass_limit','maximum_total_steps')}
    policy.update(Q_target_m3_s=Q, Umean_m_s=.002, A_in_m2=7.756795804427715e-12,
                  advective_limit_s=.5*old['h10_m']/.002,
                  viscous_limit_s=.05*old['Dh_m']**2/old['nu_m2_s'],
                  Re=.002*old['Dh_m']/old['nu_m2_s'], initial_state='zero; no restart or scaled field',
                  dt_rule='min(0.5*h10/Umean,0.05*Dh^2/nu)', maximum_wall_time_s=14400)
    (CASE/'policy.json').write_text(json.dumps(policy, indent=2)+'\n')
    manifest = {str(p.relative_to(baseline)):digest(p) for p in sorted(baseline.rglob('*')) if p.is_file()}
    (CASE/'reports/old_baseline_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
    inputs = {str(p.relative_to(CASE)):digest(p) for d in ('SV_MESH','run') for p in sorted((CASE/d).rglob('*')) if p.is_file()}
    (CASE/'input_hashes.json').write_text(json.dumps(inputs,indent=2)+'\n')
    print(json.dumps(policy,indent=2))

if __name__ == '__main__': main()
