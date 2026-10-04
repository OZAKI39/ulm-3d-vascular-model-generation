# Invariant-Based Spectral Split for Phase-Field Fracture

This repository provides a standalone C++ implementation of an invariant-based spectral decomposition for phase-field fracture modeling. The formulation follows the spectral framework introduced in:

Nardi, D.C., Ferreira, A.R., Leonel, E.D.,  
*"Revisiting Miehe’s Spectral Split: A basis-independent energy decomposition model for phase-field fracture"*, 2026.

The implementation is designed as a modular constitutive routine for computing the stress tensor and the corresponding consistent tangent operator. It can be readily coupled to standard finite element frameworks.

The repository also provides the Gmsh mesh files used in the numerical examples discussed in the reference study.

---

## 📂 Repository Structure

- `src/` — C++ source file containing the constitutive routine  
- `mesh/` — Gmsh (`.msh`) files used in the numerical examples  
- `README.md` — Project documentation  
- `LICENSE` — License information  

---

## 📊 Numerical Examples

The `mesh/` directory contains the Gmsh files corresponding to the following benchmark problems:

1. **Mode I Fracture of a Single-Edge Notched Specimen Under Tension–Compression Cycles**

2. **Mode II Fracture of a Single-Edge Notched Specimen Under Tension–Compression Cycles**

3. **L-shaped Specimen Under Alternating Tension–Compression Loading**

4. **Transversely Isotropic Single-Edge Notched Specimen: Monotonic and Cyclic Loading**

5. **Three-Point Bending of a Layered Beam with Transverse Isotropy**

These meshes reproduce the geometries and boundary-value problems analyzed in the reference study and may be directly used within independent finite element implementations.

---

- `modeI_SEN.msh`
- `modeII_SEN.msh`
- `L_shape.msh`
- `SEN_TI.msh`
- `layered_beam_TPB.msh`
