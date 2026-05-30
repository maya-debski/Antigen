# Examples

This page collects short, copy-pasteable examples for common Antigen tasks.

Prerequisites
- Install Antigen and dev extras for the CLI and examples:

```bash
pip install .[dev]
```

CLI: Prepare small dataset
- Run the prepare dataset script on a directory. Use --help to see options.

```bash
antigen_prepare_dataset.py --help
antigen_prepare_dataset.py -i /path/to/instrument/night/ -o ./out -c YYYYMMDD -m INSTRUMENT -j INSTRUMENT_ELEMENT -v
```
- Tips
  - If it cannot find a particular type of calibration file, use the appropriate flag (--arc_label, --bias_label, etc.) to specify the object name in the file header.
  - Logs will indicate which files are placed in the resulting file tree.

CLI: Perform base reduction a small dataset
- Run the base reduction script on a directory. Use --help to see options.

```bash
antigen_base_reduction.py --help
antigen_base_reduction.py -i /path/to/prepared/dataset/ -o ./reduced -c YYYYMMDD -r -m INSTRUMENT  -v
```

- Tips
  - Logs will indicate which frames are used and where products are written.

CLI: Measure response with a standard CALSPEC star
- Run the measure response script on reduced standard star data. Use --help to see options.

```bash
antigen_measure_response.py --help
antigen_measure_respons.py -r /path/to/reduced/data -o ./response -s STANDARD_STAR_IDENTIFIER -v 
```

- Tips
  - If your observed standard star is not within the CALSPEC library, you can skip this step and use a pre-determined response.
  - Logs will indicate which frames are used and where products are written.

CLI: Make data cubes
- Run the make cubes script to create final data products. Use --help to see options.

```bash
antigen_make_cubes.py --help
antigen_make_cubes.py -r /path/to/reduced/data -o ./cubes -v 
```

- Tips
  - Logs will indicate which frames are used and where products are written.

See also
- [API reference](../api/index.md)
- [CLI Usage](./cli.md)
- [Configuration](./configuration.md)
