#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from antigen.datasets import build_dataset_from_reduced_files
from antigen.manifest import save_manifest
from antigen.utils import setup_logging


# ---- A minimal SPEC describing the context this script needs ----
SPEC_CALIBRATE = [
    {"path":"cli.standard_name","source":"args","key":"standard_name","required":True,"units":"name","desc":"Spectrophotometric standard star identifier."},
    {"path":"cli.output_folder","source":"args","key":"output_folder","required":True,"units":"path","desc":"Directory for response curve and QA outputs."},
    {"path":"cli.extraction_radius","source":"args","key":"extraction_radius","required":True,"units":"pix","desc":"Aperture radius around DAR track for extraction.","validate":lambda v: float(v)>0 if v is not None else True},

    {"path":"dataset.unit_instrument","source":"dataset","key":"unit_instrument","required":True,"desc":"Instrument name (e.g., VIRUS2)."},
    {"path":"dataset.unit_id","source":"dataset","key":"unit_id","required":True,"desc":"Instrument unit identifier (e.g., VP1B)."},
    {"path":"dataset.ndithers","source":"dataset","key":"ndithers","required":True,"units":"count","desc":"Total number of dithers in the observation.","validate":lambda v: int(v)>=1 and float(v)==int(v)},
    {"path":"dataset.dither_number","source":"dataset","key":"dither_number","required":True,"units":"index","desc":"Index of this exposure in the dither pattern (1-based).","validate":lambda v: isinstance(v,(list,tuple)) and len(v)>0},
    {"path":"dataset.in_folder","source":"dataset","key":"in_folder","required":True,"default":"./","units":"path","desc":"Root folder containing reduced frames."},
    {"path":"dataset.reduced_files","source":"dataset","key":"reduced_files","required":True,"units":"list[str]","desc":"List of reduced spectral frames to process.","validate":lambda v: isinstance(v,(list,tuple)) and len(v)>0},

    {"path":"config.fiber_radius","source":"config","key":"fiber_radius","required":True,"units":"pix","desc":"Effective fiber footprint radius.","validate":lambda v: float(v)>0},
]



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

Inputs:
    - Reduced directory: Path containing reduced standard star datasets
      (already processed with bias, flats, arcs, etc.).
    - Output folder: Destination for response files (defaults to reduced directory).
    - Standard star name: Substring to identify the target (e.g., "Feige").
    - PSF and cube-building settings such as extraction radius and pixel size.

Outputs:
    - A YAML manifest file for each standard star response dataset.
    - Instrument response function FITS and diagnostic plots, written to the output folder.

Example:
    Build a response function from reduced Feige standard star frames:

    $ antigen_build_response.py -r /data/VIRUS2/reduced/ -o /data/VIRUS2/response/ -s Feige

Notes:
    - This script assumes that the input data have already been reduced
      (use the main reduction pipeline first).
    - Logs will report any datasets that do not match the requested standard.
    - Response functions are required for flux calibration of science targets.
"""




def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--reduced_dir', required=True,
                        help='Path to directory with reduced standard star frames.')
    parser.add_argument('-o', '--output_folder', default=None,
                        help='Path to output folder for response function.')
    parser.add_argument('-s', '--standard_name', required=True,
                        help='Name for standard star files in reduced-dir (ex: Feige).')

    # PSF/extraction settings
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
    logger.info(f'Starting application...')

    if args.output_folder is None:
        args.output_folder = args.reduced_dir
    save_path = Path(args.output_folder).expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    # Build response manifests
    response_manifest = build_dataset_from_reduced_files(args.reduced_dir)
    from antigen.recipe import Recipe, BuildGenericContext, PrepareData, ModelPSFAndDAR, ExtractSpectrum
    from antigen.recipe import ApplyExtinction, MeasureResponse, CalibrateContext

    for nr, manifest in enumerate(response_manifest):
        element = manifest['unit_id']
        instrument = manifest['unit_instrument']
        target_name = manifest['target']
        logger.info(f'Manifest for {target_name} with instrument={instrument} and unit={element}')
        if args.standard_name.lower() in target_name.lower():
            manifest_filename = f'manifest_response_{target_name}.yml'
            save_filepath = save_path / manifest_filename
            save_manifest(manifest, str(save_filepath))
            logger.info(f'Processing response for {target_name} with instrument={instrument} and unit={element}')

            # 1) Define the recipe
            steps = [BuildGenericContext(CalibrateContext),
                     PrepareData(),
                     ModelPSFAndDAR(),
                     ExtractSpectrum(),
                     ApplyExtinction(),
                     MeasureResponse()]
            recipe = Recipe(
                name="measure_response",
                spec=SPEC_CALIBRATE,
                steps=steps,
                outputs=["flux_cal", "response"],
                description="Measure an instrument response curve from reduced frames and a standard star.",
            )

            # 2) See the plan (no execution)
            print(recipe.plan())
            from antigen import config
            config_dict = config.build_config_for_element(
                manifest['unit_instrument'].lower(),
                manifest['unit_id'].upper()
            )
            # 3) Optional: auto-generate docs of the context + steps
            md = recipe.describe_markdown(manifest, args, config_dict)

            Path(args.output_folder).mkdir(parents=True, exist_ok=True)
            (Path(args.output_folder) / "measure_response.md").write_text(md, encoding="utf-8")

            # 4) Run the recipe
            state = recipe.run(manifest, args, config_dict, outdir=Path(args.output_folder))
    return None


if __name__ == '__main__':
    sys.exit(main())
