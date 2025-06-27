# GWAS Neale Lab Data Processing Pipeline

This repository contains tools for processing GWAS (Genome-Wide Association Studies) summary statistics from the Neale lab UK Biobank data. The pipeline provides automated data download, format conversion, and coordinate liftover for downstream polygenic risk score (PRS) analysis.

## Quick Start

### 1. Check Requirements
Run the requirements checker to ensure your system has all necessary dependencies:

```bash
./check_requirements.sh
```

### 2. Install Dependencies
Using Pixi (recommended):
```bash
pixi install
pixi shell
```

Or manually install Python packages:
```bash
pip install polars numpy
```

### 3. Download GWAS Data
```bash
# Download by phenotype name (recommended)
./download.sh --phenotype "Calcium" --output calcium
./download.sh --phenotype "Food weight" --output food_weight

# Download by phenotype code (alternative)
./download.sh --phenotype_code 100024_irnt --output calcium

# See available options
./download.sh --help
```

### 4. Process data for PRSFNN (liftover to GRCh38, remove ambiguous SNPs, etc.)
```bash
# Complete processing pipeline
./convert.sh input_file.tsv.bgz output_file.tsv

# Advanced usage with custom parameters
pixi run python convert.py -s input_file.tsv.bgz -o output_file.tsv -d /path/to/dbsnp_directory
```

## Core Features

### 1. Automated Data Download (`download.sh`)
- **Phenotype-based selection**: Search using human-readable descriptions
- **Automatic filtering**: Only downloads both_sexes and irnt-normalized phenotypes

### 2. Data Processing Pipeline (`convert.py`)
- **Coordinate Systems**: Handles both GRCh37/hg19 and GRCh38/hg38 coordinates with automatic liftover
- **SNP ID Conversion**: Converts between rsIDs and chromosome:position:allele format
- **Quality Control**: Applies filtering steps for QC

## Requirements

### System Requirements
- **Environment Manager**: Pixi (recommended) or manual Python setup

### Python Packages
- `polars`: For efficient data processing
- `numpy`: For numerical operations

### Data Resources
- **dbSNP Database**: Parquet files for rsID ↔ coordinate conversion
  - You can recreate this with [this repo](https://github.com/weinstockj/dbSNP_to_parquet).
  - Required for processing files with rsIDs lacking coordinate information
  - variant_list_{POP}.tsv; a list of variants in the LD reference panel. You can recreate with [this repo](https://github.com/weinstockj/LD_REFERENCE_PANEL). 

## Usage

### Command Line Arguments

```bash
python convert.py [OPTIONS]
```

**Required Arguments:**
- `-s, --sumstats`: Path to input GWAS summary statistics file
- `-o, --output`: Path for output reformatted file
- `-d, --dbsnp_direc`: Directory containing dbSNP parquet files

**Optional Arguments:**
- `-ib, --initial_build`: Input coordinate system (`hg19` or `hg38`, default: `hg38`)
- `-cf, --chain_file`: Path to liftOver chain file 
- `-lt, --lift_tool`: Path to liftOver executable (default: system path)
- `-vf, --variant_file`: Path to variant list file 

### Example Usage

```bash
# Process hg19 coordinates (automatic liftover to hg38)
python convert.py \
  -s old_gwas_data.tsv.bgz \
  -o processed_gwas.tsv \
  -d ~/weinstocklab/projects/resources/dbSNP_to_parquet
  -vf ~/weinstocklab/projects/resources/LD_REFERENCE_PANEL/variant_list_EUR.tsv \
  -lt liftOver \
  -ib hg19
```

## Data Processing Details

### Quality Control Filters
- Removes variants with missing essential data (BETA, SE, P)
- Filters out low-confidence variants
- Removes ambiguous SNPs (A/T and C/G pairs)
- Restricts to autosomal chromosomes (1-22)
- Applies minor allele frequency threshold (≥1% when available)
- Removes infinite or NaN values

### Coordinate Conversion
- **Input**: GRCh37/hg19 or GRCh38/hg38
- **Output**: Always GRCh38/hg38
- **Method**: UCSC liftOver for coordinate-based variants, dbSNP lookup for rsIDs
- **Chain File**: Uses `hg19ToHg38.over.chain`

## File Structure

```
GWAS_neale_lab/
├── download.sh                          # Enhanced phenotype download script
├── convert.py                          # Main data processing script
├── convert.sh                          # Pipeline wrapper script
├── neale_lab_manifest_20180731.tsv     # Phenotype metadata
├── check_requirements.sh               # System requirements checker
├── pixi.toml                           # Environment configuration
├── README.md                           # This file
└── phenotype_directories/              # Organized output folders
    ├── calcium/
    ├── food_weight/
    └── ...
```

## Logging

All processing steps are logged with timestamps. Log files are automatically created with the same basename as your output file (e.g., `output.tsv` → `output.log`).

The log includes:
- Processing statistics (variant counts, filtering results)
- Error messages and warnings
- Performance timing information

### Getting Help

1. Run the requirements checker: `./check_requirements.sh`
2. Use `--help` flag with any script for detailed usage information
3. Check log files for detailed error messages
4. Verify your input file format matches expected GWAS summary statistics

## Contact 

For questions or issues, please open an issue on the GitHub repository or contact Josh Weinstock + April Kim. 

## Acknowledgements

Most of the convert.py script was written by Seraj Grimes (with edits from April and Josh).
