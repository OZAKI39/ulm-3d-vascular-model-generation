"""Generate only the new stage's C++ constants from the frozen source artifact."""
import hashlib
import json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=root/'reference/CF2003_FREE_SHEAR_REFERENCE_V0.json'
d=json.loads(p.read_text());assert d['REFERENCE_RESOLUTION_PASS']
s='#pragma once\n#include <array>\nnamespace reducedwall::reference {\n'
s+='// Frozen JSON SHA256: '+hashlib.sha256(p.read_bytes()).hexdigest()+'\n'
for name,key in [('u','u_decimal_literals'),('omega','omega_decimal_literals')]:
    s+='inline constexpr std::array<double,20> '+name+' = {'+', '.join(d[key])+'};\n'
s+='inline constexpr double epsilon_min=0.001, epsilon_max=0.2;\n}\n'
(root/'src/cf2003_coefficients_generated.hpp').write_text(s)
print('Generated C++ constants from validated frozen JSON')
