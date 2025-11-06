#!/usr/bin/env python3
# antigen_unified_pipeline.py
"""
Unified pipeline: Base reduction → Dataset → Advanced workflows
"""
import argparse
import logging
from pathlib import Path

from antigen import recipe
from antigen import datasets

def main():
    parser = argparse.ArgumentParser(description="Unified Antigen Pipeline")
    parser.add_argument('--config', required=True, help='Base configuration file')
    parser.add_argument('--input', required=True, help='Raw data directory')
    parser.add_argument('--output', required=True, help='Output base directory')
    
    # Base reduction options
    parser.add_argument('--base-recipe', default='base_reduction.yml', 
                       help='Base reduction recipe file')
    
    # Advanced workflow options  
    parser.add_argument('--workflows', nargs='+', 
                       choices=['make_cubes', 'extract_spectrum', 'flux_calibrate'],
                       default=['make_cubes'],
                       help='Advanced workflows to run')
    
    # Workflow-specific parameters (these get passed through)
    parser.add_argument('--pixel-size', type=float, default=1.0)
    parser.add_argument('--interpolation-method', default='gdw')
    # ... other workflow parameters
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('antigen.unified')
    
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Phase 1: Base reduction using existing recipe system
    logger.info("Phase 1: Running base reduction")
    base_output = output_path / "reduced"
    
    # Load and run base reduction recipe
    base_recipe_loader = recipe.Recipe(args.base_recipe)
    base_inputs = {
        'cli': {
            'config_file': args.config,
            'input_folder': args.input,
            'output_folder': str(base_output)
        }
    }
    base_results = base_recipe_loader.load(base_inputs)
    
    # Phase 2: Build dataset from reduced files
    logger.info("Phase 2: Building dataset from reduced files") 
    
    # Extract reduced files list from base results
    reduced_files = base_results.get('reduced_files', [])
    if not reduced_files:
        # Fallback: find FITS files in output directory
        reduced_files = list(base_output.glob("*.fits"))
        reduced_files = [f.name for f in reduced_files]
    
    dataset = datasets.build_dataset_from_reduced_files(
        base_folder=str(base_output),
        filenames=reduced_files,
        config_file=args.config
    )
    
    # Phase 3: Run advanced workflows
    for workflow_name in args.workflows:
        logger.info(f"Phase 3: Running {workflow_name} workflow")
        
        workflow_output = output_path / workflow_name
        
        # Load workflow recipe
        workflow_recipe_file = f"{workflow_name}.yml"
        workflow_recipe_loader = recipe.Recipe(workflow_recipe_file)
        
        # Prepare workflow inputs by combining dataset + CLI args
        workflow_inputs = {
            'dataset': dataset,  # This contains reduced file paths, instrument info, etc.
            'cli': {
                'output_folder': str(workflow_output),
                'pixel_size': args.pixel_size,
                'interpolation_method': args.interpolation_method,
                # Add other CLI parameters as needed
            },
            'config': {
                # Config parameters come from the dataset or config file
            }
        }
        
        # Run the workflow
        workflow_results = workflow_recipe_loader.load(workflow_inputs)
        logger.info(f"Completed {workflow_name}: {workflow_results}")
    
    logger.info("Unified pipeline completed successfully")

if __name__ == "__main__":
    main()
