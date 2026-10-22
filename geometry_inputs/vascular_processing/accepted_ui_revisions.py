"""Explicit September UI revision compatibility; never authorizes data changes."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
REVISIONS={
 's1-2_swc_roi_generate_human.py':('3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729','fa07775019b6cbbbaa5c89e775997f462f0dfa308c1aec7646d9002c00dc559e'),
 's1-3_swc_roi_generate_MeVO.py':('348f4e9ae709f6d9aa934bb00d56a5881269dca3648bb13508f6f997347fccfa','b72db7638ced07d153b0612f657d072d8df5d9018406cb0174761e1eeb90c002'),
 'utils/rodent_vasculature/interactive.py':('42c648e192240548c7ee7c1934f8a026ecd52bf1b23b1b6c2909e075bd4a67b2','f23e935bacb41f7f7ed75669b7793bd714b9cc4eca36665114c475bff58d80a7')}
def approved_revision(path,before,after):
    try:name=str(Path(path).resolve().relative_to(ROOT))
    except ValueError:return False
    if name == 's1-3_swc_roi_generate_MeVO.py' and before == '921fe37d18d104481fcc37c38735eb59bfbe5f49d548d35ebf002a0516b0e791':
        return after == REVISIONS[name][1]
    return REVISIONS.get(name)==(before,after)
