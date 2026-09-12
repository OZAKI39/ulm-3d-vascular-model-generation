# BC API feasibility — PASS (public API; Stage A runtime verification required)

HemoCell current bundled Palabos provides arbitrary-triangle Guo off-lattice boundary reconstruction. The official aneurysm example uses velocity inlets and constant-density pressure outlets. No new population reconstruction is introduced.

Cartesian addPressureBoundary0P/0N/1P/1N/2P/2N specialize missing populations by coordinate direction and require a planar Cartesian Box3D. Small boxes do not remove that assumption. They are not used here.

Only frozen Step2 cap triangle IDs receive port tags; all other surface triangles remain stationary no-slip. No hole discovery, synthetic circles, new cap, or box-face port is used.

Step2 closed lumen defines physical volume. Step2 opened mask and all 497 port labels are preserved as verified port diagnostics. Native Guo completes outerBorder ghost nodes (official aneurysm arrangement); they are numerical support, not added physical lumen. Runtime checks must prove unchanged closed flags and cap IDs. A width-2 envelope changes storage only.

## inlet
Normal: [-7.007266199383469e-07, -6.027718016794421e-07, 0.9999999999995729]; frozen label 1, 56 cap triangles. Target: target_volumetric_flow. Public API: **YES**, Q-normalized native VelocityPlugProfile3D. Uses actual surface normals, without dominant-axis substitution. Finite-grid cap/wall junction errors and short transients remain to be checked.

## outlet_01
Normal: [0.018388404769920114, 0.8046101061009614, 0.5935186970350784]; frozen label 2, 44 cap triangles. Target: pressure. Public API: **YES**, native DensityNeumannBoundaryProfile3D: fixed rho, normal momentum extrapolation. Uses actual surface normals, without dominant-axis substitution. Finite-grid cap/wall junction errors and short transients remain to be checked.

## outlet_02
Normal: [-0.9010100696103183, -0.2746666844492976, -0.3357663874101153]; frozen label 3, 42 cap triangles. Target: pressure. Public API: **YES**, native DensityNeumannBoundaryProfile3D: fixed rho, normal momentum extrapolation. Uses actual surface normals, without dominant-axis substitution. Finite-grid cap/wall junction errors and short transients remain to be checked.

## outlet_03
Normal: [0.9913380590746697, 0.12599709605873874, -0.03706189977364442]; frozen label 4, 49 cap triangles. Target: pressure. Public API: **YES**, native DensityNeumannBoundaryProfile3D: fixed rho, normal momentum extrapolation. Uses actual surface normals, without dominant-axis substitution. Finite-grid cap/wall junction errors and short transients remain to be checked.

## Local source evidence
- `/home/lzy/projects/hemocell_starter/examples/pipeflow/pipeflow.cpp`: periodic/body-force pipe; unsuitable physical forcing here (SHA256 e071b25952730ce8c5dc9ebaafc07239527941d8c9fd8ba7e73d5f748c3b6a59).
- `/home/lzy/projects/hemocell_starter/examples/pipeflow_with_preinlet/pipeflow_with_preinlet.cpp`: preinlet and Cartesian 0N pressure (SHA256 65f5f7bc157596842b3020d91a7b7252fe52abcf63f6a6d1c847bba2be8d377f).
- `/home/lzy/projects/hemocell_starter/examples/curvedflow_with_preinlet/curvedflow_with_preinlet.cpp`: hardcoded Cartesian pressure box; not arbitrary cap (SHA256 dce7e4fe0bc89a31bfdaf216828dc21b27c411a0105ca86091252856c41f2233).
- `/home/lzy/projects/hemocell_starter/cases/stl_preinlet/stl_preinlet.cpp`: STL with Cartesian preinlet pressure (SHA256 dc986f70fce0e793036ad22eb355aa96dcee0de2d5ef0bc2ed61710c1113776f).
- `/home/lzy/projects/hemocell_starter/cases/AR2/AR2.cpp`: Cartesian 0N pressure reference (SHA256 053e8edd52459e40ba4b4a59aa88144c88c23f563bc2f03349b7badf7fac8ded).
- `/home/lzy/projects/hemocell_starter/cases/preinlet_shear/preinlet_shear.cpp`: Cartesian 0P pressure / planar velocity (SHA256 4756c9b00d2119919849178b6f06b976be4c0855a26efe6c4b803669c096050f).
- `/home/lzy/projects/hemocell_starter/palabos/src/boundaryCondition/boundaryInstantiator3D.h`: 1037-1040 planar domain precondition; 1272-1312 directions 0/1/2 orientations +/- (SHA256 6ecd29e8617767ee2d03a01e5d212a91f31ef7a7fe8025a3c601e7397e20c356).
- `/home/lzy/projects/hemocell_starter/palabos/src/boundaryCondition/zouHeBoundary3D.hh`: WrappedZouHeBoundaryManager pressure is direction/orientation-specialized (SHA256 29518d68768a1462bd1e9177e462a17018ce8609e4be208f15c64e01b58c7624).
- `/home/lzy/projects/hemocell_starter/palabos/examples/showCases/aneurysm/aneurysm.cpp`: 95-137,285-365: native Guo + plug/Poiseuille + constant-density Neumann-velocity multi-opening flow (SHA256 0f0d2185c0b4bca4aaed0c70cebac6b675ee7f3a5456422a1f231c64c6a1e670).
- `/home/lzy/projects/hemocell_starter/palabos/src/offLattice/triangleBoundary3D.h`: public tagDomain, getTag, BoundaryProfiles3D.defineProfile (SHA256 ad7bc10a7a7e95c7948502c51f9b4681fcae296863e5b61e771b1f67d6babd05).
- `/home/lzy/projects/hemocell_starter/palabos/src/offLattice/triangleBoundary3D.hh`: 724+: tagDomain only triangles whose three vertices pass; TriangleFlowShape exact triangle intersections and profiles (SHA256 9a3c13ad993c4b0b1d9f5c6bbfd5237acd5223c601dbeb15c78cb0af7c3e6e38).
- `/home/lzy/projects/hemocell_starter/palabos/src/offLattice/offLatticeBoundaryProfiles3D.hh`: VelocityPlugProfile; 527-556 DensityNeumann fixed density, NOT zero-gradient pressure (SHA256 ebd012d0ed8992ada95f47024fdb4f9ea5b379a06554ae6883a485aabfcfd3f8).
- `/home/lzy/projects/hemocell_starter/palabos/src/offLattice/guoOffLatticeModel3D.hh`: 114-190 dry-node stencil; 485+ density override; 513+ normal momentum; existing population reconstruction (SHA256 2af6eb135ac05e1089e72f140c55260a52b4bf9e05b028d889a2db207890a94f).
- `/home/lzy/projects/hemocell_starter/palabos/src/offLattice/offLatticeBoundaryCondition3D.hh`: insert() no packed-field arguments supported; computePressure subtracts MEAN density, unsuitable for our fixed gauge reference (SHA256 00caee998b86faf32e1403685b7f93177126ea8136b8b9bebe3911a750153804).
- `/home/lzy/projects/hemocell_starter/mechanics/constantConversion.cpp`: 36-54: nuLU=(tau-.5)/3, dt=nuLU*dx²/nuPhys, dm=rhoPhys*dx³, df=dm*dx/dt² (SHA256 28cf72c29d22dd6dc7357d1237051fdc07131bf9c59beb1b891cbafa1c5fe941).
- `/home/lzy/projects/hemocell_starter/palabos/src/latticeBoltzmann/nearestNeighborLattices3D.hh`: D3Q19 weights, velocities, cs²=1/3 (SHA256 b3463c440b222aa8a4fe1e658e66d3b5a1ac29f4ff8aaebc0905cef660f681d3).
- `/home/lzy/projects/hemocell_starter/palabos/src/latticeBoltzmann/dynamicsTemplates.h`: BGK second-order equilibrium; second moment = cs²*rho*I + rho*u*u (SHA256 df0263344a10114e26311c32542278e190b90b8891af66df8858b9ae1db5ad18).
- `/home/lzy/projects/hemocell_starter/palabos/src/basicDynamics/isoThermalDynamics.hh`: 360-365 BGKdynamics::computeEquilibrium delegates to bgk_ma2_equilibrium (SHA256 898118394f0598dde2cc003c30443b833e7fdfa72c02598aa8a975d40d4b8fc6).
