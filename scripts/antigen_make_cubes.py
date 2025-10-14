#!/usr/bin/env python
from pathlib import Path
import argparse

from antigen.datasets import build_dataset_from_reduced_files
from antigen.manifest import save_manifest
from antigen.utils import setup_logging
from antigen import config
from antigen.recipe import Recipe
from antigen.cli import add_common_args

DESCRIPTION = r"""
Purpose:
    This script creates 3D datacubes from reduced fiber spectra with DAR correction.
    It processes fiber spectra and creates a spatial datacube using interpolation methods.
"""


def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_args(parser)

    parser.add_argument(
        '-m', '--instrument',
        type=str,
        default='GCMS',
        help='Name of the instrument used for the observation '
             '(default: %(default)s). Example: GCMS.'
    )
    parser.add_argument(
        '-p', '--pixel_size',
        type=float,
        default=1.0,
        help='Pixel size in arcseconds for output cube '
             '(default: %(default)s). Example: 1.0.'
    )
    parser.add_argument(
        '-i', '--interpolation_method',
        type=str,
        default='linear',
        choices=['linear', 'nearest', 'cubic', 'rbf', 'gdw'],
        help='Interpolation method for cube construction '
             '(default: %(default)s). Options: linear, nearest, cubic, rbf, gdw.'
    )
    parser.add_argument(
        '-k', '--k_neighbors',
        type=int,
        default=5,
        help='Number of neighbors for GDW interpolation '
             '(default: %(default)s). Example: 5.'
    )
    parser.add_argument(
        '-s', '--sigma_gdw',
        type=float,
        default=1.0,
        help='Standard deviation for GDW interpolation '
             '(default: %(default)s). Example: 1.0.'
    )
    parser.add_argument(
        '-r', '--rbf_function',
        type=str,
        default='multiquadric',
        help='Radial basis function type for RBF interpolation '
             '(default: %(default)s). Example: multiquadric.'
    )
    parser.add_argument(
        '--reduced_files',
        nargs='+',
        required=True,
        help='List of reduced spectrum files to process.'
    )

    return parser.parse_args()


def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose, debug=args.debug)
    logger.info('Starting cube creation...')

    save_path = Path(args.out_folder or './cubes').expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    # Process manifests
    cube_manifests = build_dataset_from_reduced_files(args.reduced_dir)

    for manifest in cube_manifests:
        
        # Save manifest
        manifest_filename = f'manifest_make_cubes.yml'
        save_filepath = save_path / manifest_filename
        save_manifest(manifest, str(save_filepath))

        logger.info(f'Processing cube creation for {manifest_filename}')

        # Get instrument configuration
        config_dict = config.build_config_for_element(
            manifest['unit_instrument'].lower(),
            manifest['unit_id'].upper()
        )

        # Load recipe
        base_path = config.get_base_config_path()
        recipe = Recipe.load("make_cubes", base_path)

        # Collect inputs
        inputs = recipe.collect_inputs(args, manifest, config_dict)

        # Validate inputs
        if errors := recipe.validate_inputs(inputs):
            raise ValueError("Input validation failed:\n" + "\n".join(errors))

        # Generate and save markdown description
        md = recipe.describe_markdown(inputs, viz_type='mermaid', output_folder=save_path)
        (save_path / "make_cubes.md").write_text(md, encoding="utf-8")

        # Run recipe
        outputs = recipe.run(inputs, save_path)
        
        logger.info(f"Cube creation completed. Output saved to {save_path}")

    return None


if __name__ == "__main__":
    main()
