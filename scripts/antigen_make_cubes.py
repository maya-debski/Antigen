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

    parser.add_argument('-r', '--reduced_dir', required=True,
                       help='Path to directory with reduced standard star frames.')
    parser.add_argument('-o', '--output_folder', default=None,
                       help='Path to output folder for response function.')
    parser.add_argument('-s', '--object_name', required=False, default=None,
                       help='Name for standard star files in reduced-dir (ex: Feige).')
    parser.add_argument('-e', '--extraction_radius', type=float, default=10.0,
                       help='Extraction radius in arcsec.')
    parser.add_argument('-p', '--pixel_size', type=float, default=1.0,
                       help='Cube pixel size in arcsec.')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='if True, print more process details and logger.info to console')
    parser.add_argument('-d', '--debug', action='store_true',
                       help='if True, print more process details and logger.debug to console')
    parser.add_argument('-i', '--interpolation_method', type=str, default='gdw',
                        choices=['linear', 'nearest', 'cubic', 'rbf', 'gdw'],
                        help='Interpolation method for cube construction (default: %(default)s)')
    parser.add_argument('-k', '--k_neighbors', type=int, default=11,
                        help='Number of neighbors for GDW interpolation (default: %(default)s)')
    parser.add_argument('--sigma_gdw', type=float, default=2.5,
                        help='Standard deviation for GDW interpolation in arcseconds (default: %(default)s)')
    parser.add_argument('--rbf_function', type=str, default='multiquadric',
                        choices=['multiquadric', 'inverse', 'gaussian', 'linear', 'cubic', 'quintic', 'thin_plate'],
                        help='Radial basis function type for RBF interpolation (default: %(default)s)')
    parser.add_argument('-g', '--binned', action='store_true', default=False,
                        help='Data is binned in the x-direction?')
    # VIRUS2: control dropping dedicated sky fibers (head_id like S1, _S17_)
    parser.add_argument('--no-drop-sky-fibers-in-cubes', dest='drop_sky_fibers_in_cubes', action='store_false',
                        help='Do not drop VIRUS2 dedicated sky fibers when building cubes (default is to drop them).')
    parser.set_defaults(drop_sky_fibers_in_cubes=True)
    return parser.parse_args()


def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose, debug=args.debug)
    logger.info('Starting cube creation...')

    save_path = Path(args.output_folder or './cubes').expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    # Process manifests
    cube_manifests = build_dataset_from_reduced_files(args.reduced_dir)

    # If an object name is provided, only process manifests matching that target
    if args.object_name:
        filtered = [m for m in cube_manifests if args.object_name.lower() in m.get('target', '').lower()]
        if not filtered:
            logger.warning(f"No manifests matched object_name='{args.object_name}'. Nothing to do.")
            return None
        logger.info(f"Filtering manifests by object_name='{args.object_name}'. {len(filtered)} match(es) will be processed.")
        cube_manifests = filtered

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

        # Binning?
        if args.binned == True:
            config_dict['detector_dimensions']['X'] = int(config_dict['detector_dimensions']['X'] / 2)

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
