#!/usr/bin/env python
from pathlib import Path
import argparse

from antigen.datasets import find_datasets
from antigen.manifest import save_manifest, read_manifest
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
        '-M', '--manifest',
        type=str,
        default=None,
        help='Full path to a YAML manifest file to use instead of auto-discovery (find_datasets). If provided, in_folder/obs_date/etc. are ignored.'
    )

    parser.add_argument(
        '-m', '--instrument',
        type=str,
        default='GCMS',
        help='Name of the instrument used for the observation '
             '(default: %(default)s). Example: GCMS.'
    )
    parser.add_argument(
        '-u', '--unit-id',
        type=str,
        default=None,
        help='If provided, only process datasets whose manifest unit_id matches this value (case-insensitive).'
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

    base_save_path = Path(args.out_folder or args.reduced_dir).expanduser().resolve()
    base_save_path.mkdir(parents=True, exist_ok=True)

    # Process manifests
    dataset_manifests = find_datasets(args.in_folder, args.obs_date, args.obs_name, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label, instrument=args.instrument)
    for manifest in dataset_manifests:
        # If a specific unit-id was requested, skip non-matching manifests
        if getattr(args, 'unit_id', None):
            try:
                m_unit = str(manifest.get('unit_id', '')).upper()
            except Exception:
                m_unit = ''
            if m_unit != str(args.unit_id).upper():
                logger.info(f"Skipping manifest for unit_id={m_unit or 'UNKNOWN'} (requested {args.unit_id})")
                continue
        # Create a per-spectrograph output subfolder to avoid mixing outputs
        try:
            sub_name = f"{manifest['unit_instrument'].upper()}_{manifest['unit_id'].upper()}"
        except Exception:
            # Fallback to unit_id only if instrument missing
            sub_name = str(manifest.get('unit_id', 'UNKNOWN')).upper()
        out_dir = base_save_path / sub_name
        out_dir.mkdir(parents=True, exist_ok=True)

        # Save manifest in the subfolder
        manifest_filename = f'manifest_base_reduction.yml'
        save_filepath = out_dir / manifest_filename
        save_manifest(manifest, str(save_filepath))

        logger.info(f'Processing base reduction for {manifest_filename} in {out_dir}')

        # Get instrument configuration
        try:
            config_dict = config.build_config_for_element(
                manifest['unit_instrument'].lower(),
                manifest['unit_id'].upper()
            )
        except Exception as e:
            logger.error(f"Error building config for {manifest['unit_id']}: {e}")
            continue

        # Binning?
        if args.binned == True:
            config_dict['detector_dimensions']['X'] = int(config_dict['detector_dimensions']['X'] / 2)

        # Load recipe
        base_path = config.get_base_config_path()
        recipe = Recipe.load("base_reduction", base_path)

        # Ensure recipe steps use the per-spectrograph folder by overriding CLI out_folder
        old_out_folder = getattr(args, 'out_folder', None)
        args.out_folder = str(out_dir)
        try:
            # Collect inputs
            inputs = recipe.collect_inputs(args, manifest, config_dict)

            # Validate inputs
            if errors := recipe.validate_inputs(inputs):
                raise ValueError("Input validation failed:\n" + "\n".join(errors))

            # Generate and save markdown description
            md = recipe.describe_markdown(inputs, viz_type='mermaid', output_folder=out_dir)
            (out_dir / "base_reduction.md").write_text(md, encoding="utf-8")

            # Run recipe
            outputs = recipe.run(inputs, out_dir)
        finally:
            # Restore original out_folder for safety before next loop iteration
            args.out_folder = old_out_folder

    return None


if __name__ == "__main__":
    main()