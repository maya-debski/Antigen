#!/usr/bin/env python
from pathlib import Path
import argparse

from antigen.datasets import build_dataset_from_reduced_files
from antigen.manifest import save_manifest
from antigen.utils import setup_logging
from antigen import config
from antigen.recipe import Recipe

DESCRIPTION = r"""
Purpose:
    This script builds an instrument response function from reduced
    standard star observations for the VIRUS2/GCMS instrument.

What it does:
    - Scans a directory of reduced datasets for standard star frames.
    - Generates a dataset manifest for each candidate standard star.
    - Selects the datasets matching the requested standard star name.
    - For each match:
        * Builds an instrument response function by comparing the observed
          spectrum with a reference CALSPEC spectrum.
        * Saves the response manifest and response products to the output folder.
"""

def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, 
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--reduced_dir', required=True,
                       help='Path to directory with reduced standard star frames.')
    parser.add_argument('-o', '--output_folder', default=None,
                       help='Path to output folder for response function.')
    parser.add_argument('-s', '--standard_name', required=True,
                       help='Name for standard star files in reduced-dir (ex: Feige).')
    parser.add_argument('-e', '--extraction_radius', type=float, default=10.0,
                       help='Extraction radius in arcsec.')
    parser.add_argument('-p', '--pixel_size', type=float, default=1.0,
                       help='Cube pixel size in arcsec.')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='if True, print more process details and logger.info to console')
    return parser.parse_args()

def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose)
    logger.info('Starting application...')

    save_path = Path(args.output_folder or args.reduced_dir).expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    # Process manifests
    response_manifest = build_dataset_from_reduced_files(args.reduced_dir)

    for manifest in response_manifest:
        if args.standard_name.lower() in manifest['target'].lower():
            # Save manifest
            manifest_filename = f'manifest_response_{manifest["target"]}.yml'
            save_filepath = save_path / manifest_filename
            save_manifest(manifest, str(save_filepath))
            
            logger.info(f'Processing response for {manifest["target"]}')
            
            # Get instrument configuration
            config_dict = config.build_config_for_element(
                manifest['unit_instrument'].lower(),
                manifest['unit_id'].upper()
            )

            # Load recipe
            base_path = config.get_base_config_path()
            recipe = Recipe.load("measure_response", base_path)
            
            # Collect inputs
            inputs = recipe.collect_inputs(args, manifest, config_dict)
            
            # Validate inputs
            if errors := recipe.validate_inputs(inputs):
                raise ValueError("Input validation failed:\n" + "\n".join(errors))
            
            # Generate and save markdown description
            md = recipe.describe_markdown(inputs, viz_type='mermaid', output_folder=save_path)

            (save_path / "measure_response.md").write_text(md, encoding="utf-8")
            
            # Run recipe
            outputs = recipe.run(inputs, save_path)

    return None
    
if __name__ == "__main__":
    main()