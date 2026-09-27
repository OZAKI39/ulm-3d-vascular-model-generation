# Retired TopBrain tree pipeline

Skan was evaluated to retain Lee medial samples while using an established branch representation. VMTK was evaluated as an independent surface-based centerline method. Neither route established a defensible production tree for the strict TopBrain M2/M3 anatomy.

The strict mask itself contains digital topology loops. Root localization and cycle removal are separate questions: an anatomical root can lie inside an existing path, but subdividing that path does not remove loops. No lossless automatic tree repair was justified. The canonical VMTK check still failed anatomical coverage and lumen containment.

TopBrain-to-tree conversion, its comparison and forensic tools, tests and isolated environment have therefore been retired. The new route transfers TopBrain M1/M2/M3 semantics with Open3D similarity alignment and POT Partial Fused Gromov-Wasserstein onto the original BraVa tree. BraVa remains the source of coordinates, radii, connectivity and modeling geometry.

Native dataset discovery, ITK-SNAP label parsing, affine coordinates, source provenance, mask display and the shared renderer remain supported. The independent CFD surface-extension workflow still uses VMTK and is not part of this retirement. Original datasets and established UI outputs are retained. Deletion inventory and preserved source provenance are recorded under outputs/topbrain_brava_transfer/.
