#!/usr/bin/env python
from pathlib import Path
import argparse

from antigen.datasets import find_datasets
from antigen.manifest import save_manifest
from antigen.utils import setup_logging
from antigen import config
from antigen.recipe import Recipe
from antigen.cli import add_calibration_args, add_common_args

DESCRIPTION = r"""
Purpose:
    This script does a basic reduction for GCMS/VIRUS2/VIRUSW
.
"""


def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_args(parser)
    add_calibration_args(parser)

    parser.add_argument(
        '-m', '--instrument',
        type=str,
        default='GCMS',
        help='Name of the instrument used for the observation '
             '(default: %(default)s). Example: GCMS.'
    )
    parser.add_argument(
        '-j', '--good_arc_residual_limit',
        type=float,
        default=0.2,
        help='Residual limit for an arc line to be considered well-fit'
             '(default: %(default)s). Example: 0.2.'
    )
    parser.add_argument(
        '-g', '--binned',
        action='store_true',
        default=False,
        help='Data is binned in the x-direction?'
    )

    return parser.parse_args()


def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose, debug=args.debug)
    logger.info('Starting application...')

    save_path = Path(args.out_folder or args.reduced_dir).expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    # Process manifests
    dataset_manifests = find_datasets(args.in_folder, args.obs_date, args.obs_name, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label, instrument=args.instrument)
    for manifest in dataset_manifests:
        # Save manifest
        manifest_filename = f'manifest_base_reduction.yml'
        save_filepath = save_path / manifest_filename
        save_manifest(manifest, str(save_filepath))

        logger.info(f'Processing base reduction for {manifest_filename}')

        # Get instrument configuration
        config_dict = config.build_config_for_element(
            manifest['unit_instrument'].lower(),
            manifest['unit_id'].upper()
        )

        # Load recipe
        base_path = config.get_base_config_path()
        recipe = Recipe.load("base_reduction", base_path)

        # Collect inputs
        inputs = recipe.collect_inputs(args, manifest, config_dict)

        # Validate inputs
        if errors := recipe.validate_inputs(inputs):
            raise ValueError("Input validation failed:\n" + "\n".join(errors))

        # Generate and save markdown description
        md = recipe.describe_markdown(inputs, viz_type='mermaid', output_folder=save_path)

        (save_path / "base_reduction.md").write_text(md, encoding="utf-8")

        # Run recipe
        outputs = recipe.run(inputs, save_path)

    return None


if __name__ == "__main__":
    main()