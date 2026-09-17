#!/usr/bin/env python
"""
The main entry point for the DCCO vascular generator.
"""

import argparse
from pathlib import Path
import sys
if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ulm_vascular_model_generator.utils.config.yaml_config import DEFAULT_CONFIG_PATH, load_generator_config
from ulm_vascular_model_generator.utils.core.tree import DCCOTree
from ulm_vascular_model_generator.utils.generation.stages import generate_staged_tree
from ulm_vascular_model_generator.utils.io.generation_validation import write_generation_validation_outputs
from ulm_vascular_model_generator.utils.io.run_output import create_run_directory, sync_swc_to_output_root
from ulm_vascular_model_generator.utils.io.swc import summarize, vessels_to_swc, write_swc
from ulm_vascular_model_generator.utils.io.vessel_transport_export import vessel_transport_path_from_swc, write_vessel_transport_npz


def main():
    # ============================
    # Parse command line arguments
    # ============================
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH, 
                        help=f"YAML configuration file. Default: {DEFAULT_CONFIG_PATH}")
    args = parser.parse_args()

    # ============================
    # Load configuration and generate vessels
    # ============================
    output, cfg, stages = load_generator_config(args.config)
    if stages:
        tree = generate_staged_tree(stages)
    else:
        tree = DCCOTree(cfg)
        tree.generate()
    vessels = tree.vessels

    # ===========================
    # Write output files
    # ===========================
    synced_swc_output = output
    run_directory = create_run_directory(
        synced_swc_output.parent,
        tree.cfg.seed,
        tree.cfg.geometry_mode,
    )
    output = run_directory / synced_swc_output.name
    nodes = vessels_to_swc(vessels)
    write_swc(nodes, output)
    vessel_data_output = vessel_transport_path_from_swc(output)
    transport_metadata = {
        "config_path": str(args.config),
        "run_directory": str(run_directory),
        "swc_path": str(output),
        "synced_swc_path": str(synced_swc_output),
        "geometry_mode": tree.cfg.geometry_mode,
        "pressure_and_shear": "not_computed_by_vascular_generator",
    }
    transport_metadata.update(tree.haemodynamics.transport_metadata())
    write_vessel_transport_npz(
        vessels, vessel_data_output,
        metadata=transport_metadata,
    )
    validation_json, validation_csv = write_generation_validation_outputs(
        vessels,
        tree.cfg,
        run_directory,
    )
    sync_swc_to_output_root(output, synced_swc_output)
    print(f"Read config {args.config}")
    print(f"Created result directory {run_directory}")
    print(f"Wrote {output}")
    print(f"Updated {synced_swc_output}")
    print(f"Wrote {vessel_data_output}")
    print(f"Wrote {validation_json}")
    print(f"Wrote {validation_csv}")
    print(summarize(nodes))
    print(
        f"  vessels={len(vessels)}, "
        f"terminals={sum(not vessel.children for vessel in vessels)}"
    )

if __name__ == "__main__":
    main()
