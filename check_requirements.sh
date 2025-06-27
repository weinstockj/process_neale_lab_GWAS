#!/bin/bash

# check_requirements.sh - Verify system requirements for GWAS processing pipeline

set -e

echo "Checking system requirements for GWAS processing pipeline..."
echo "================================================================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track if all requirements are met
ALL_GOOD=true

# Function to check if command exists
check_command() {
    local cmd="$1"
    local name="${2:-$cmd}"
    if command -v "$cmd" &> /dev/null; then
        echo -e "${GREEN}✓${NC} $name found: $(command -v "$cmd")"
        return 0
    else
        echo -e "${RED}✗${NC} $name not found"
        ALL_GOOD=false
        return 1
    fi
}

# Function to check Python package (preferring pixi if available)
check_python_package() {
    local package="$1"
    local python_cmd="python"
    
    # Use pixi python if pixi is available and pixi.toml exists
    if command -v pixi &> /dev/null && [ -f "pixi.toml" ]; then
        python_cmd="pixi run python"
    fi
    
    if $python_cmd -c "import $package" &> /dev/null; then
        local version=$($python_cmd -c "import $package; print(getattr($package, '__version__', 'unknown'))" 2>/dev/null || echo "unknown")
        echo -e "${GREEN}✓${NC} Python package '$package' found (version: $version)"
        return 0
    else
        echo -e "${RED}✗${NC} Python package '$package' not found"
        ALL_GOOD=false
        return 1
    fi
}

# Check Python
echo
echo "Checking Python..."
echo "-------------------"

# Determine which python command to use
PYTHON_CMD="python"
PYTHON_CONTEXT="system"

# Check if pixi is available and pixi.toml exists
if command -v pixi &> /dev/null && [ -f "pixi.toml" ]; then
    PYTHON_CMD="pixi run python"
    PYTHON_CONTEXT="pixi environment"
    echo "Using Python from pixi environment"
fi

if $PYTHON_CMD --version &> /dev/null; then
    PYTHON_VERSION=$($PYTHON_CMD --version 2>&1 | cut -d' ' -f2)
    echo -e "${GREEN}✓${NC} Python found in $PYTHON_CONTEXT: $PYTHON_VERSION"
    
    # Check if Python version is >= 3.13
    if $PYTHON_CMD -c "import sys; exit(0 if sys.version_info >= (3, 13) else 1)" 2>/dev/null; then
        echo -e "${GREEN}✓${NC} Python version meets requirements (>= 3.13)"
    else
        echo -e "${YELLOW}⚠${NC} Python version should be >= 3.13 (found: $PYTHON_VERSION)"
    fi
else
    echo -e "${RED}✗${NC} Python not found. Please install Python >= 3.13"
    if [ "$PYTHON_CONTEXT" = "pixi environment" ]; then
        echo "  Try running: pixi install"
    fi
    ALL_GOOD=false
fi

# Check Python packages
echo
echo "Checking Python packages..."
echo "----------------------------"
check_python_package "polars"
check_python_package "numpy"

# Check for standard library modules (should always be available)
echo
echo "Checking Python standard library modules..."
echo "--------------------------------------------"
for module in subprocess argparse re gzip os logging; do
    check_python_package "$module"
done

# Check Pixi (environment manager)
echo
echo "Checking environment manager..."
echo "-------------------------------"
if check_command pixi; then
    echo "  Pixi is the recommended environment manager for this project."
    if [ -f "pixi.toml" ]; then
        echo -e "${GREEN}✓${NC} pixi.toml configuration file found"
    else
        echo -e "${YELLOW}⚠${NC} pixi.toml not found in current directory"
    fi
else
    echo -e "${YELLOW}⚠${NC} Pixi not found. Consider installing Pixi for easier dependency management."
    echo "  Install from: https://pixi.sh/"
fi

# Check for convert.py script
echo
echo "Checking main script..."
echo "-----------------------"
if [ -f "convert.py" ]; then
    echo -e "${GREEN}✓${NC} convert.py found"
    
    # Check if script is executable
    if [ -x "convert.py" ]; then
        echo -e "${GREEN}✓${NC} convert.py is executable"
    else
        echo -e "${YELLOW}⚠${NC} convert.py is not executable (this is fine for 'python convert.py')"
    fi
else
    echo -e "${RED}✗${NC} convert.py not found in current directory"
    ALL_GOOD=false
fi

# Check for dbSNP directory (optional but recommended)
echo
echo "Checking optional resources..."
echo "------------------------------"
DBSNP_DIR="$HOME/weinstocklab/projects/resources/dbSNP_to_parquet/output"
if [ -d "$DBSNP_DIR" ]; then
    echo -e "${GREEN}✓${NC} dbSNP parquet directory found: $DBSNP_DIR"
    PARQUET_COUNT=$(find "$DBSNP_DIR" -name "*.parquet" 2>/dev/null | wc -l)
    echo "  Found $PARQUET_COUNT parquet files"
else
    echo -e "${YELLOW}⚠${NC} dbSNP parquet directory not found at: $DBSNP_DIR"
    echo "  This is needed for rsID to coordinate conversion"
fi

# Summary
echo
echo "Summary"
echo "======="
if [ "$ALL_GOOD" = true ]; then
    echo -e "${GREEN}✓ All core requirements are satisfied!${NC}"
    echo
    echo "You can run the conversion script with:"
    echo "  python convert.py -s <input_file> -o <output_file> -d <dbsnp_directory>"
    echo
    echo "Or with Pixi (recommended):"
    echo "  pixi run python convert.py -s <input_file> -o <output_file> -d <dbsnp_directory>"
    exit 0
else
    echo -e "${RED}✗ Some requirements are missing. Please install the missing components.${NC}"
    echo
    echo "To install dependencies with Pixi:"
    echo "  pixi install"
    echo
    echo "To install Python packages manually:"
    echo "  pip install polars numpy"
    exit 1
fi