# scTenifoldKnk

Optimized version of scTenifoldKnk for in-silico gene knockout experiments using single-cell RNA-seq data.

## What's New in Version 2.0.0

**Maintainer**: Zaoqu Liu (liuzaoqu@163.com)  
**Repository**: https://github.com/Zaoqu-Liu/scTenifoldKnk

### Key Improvements

- **3-4x Performance Boost**: For typical datasets (500-1000 genes)
  - 500 genes: 2.2 min → 0.6 min (4x faster)
  - 1000 genes: 8.6 min → 2.5 min (3.5x faster)
  
- **Perfect Accuracy**: Correlation > 0.999 with original implementation

- **Smart Parallelization**: Automatically adapts to dataset size
  - Small datasets (<200 genes): Sequential processing to avoid overhead
  - Large datasets (>200 genes): Multi-core parallel processing
  
- **Optimized Algorithms**:
  - Efficient matrix operations
  - Vectorized computations
  - Memory-optimized sparse matrix handling

### Performance Comparison

| Data Size | Original | Optimized (v2.0) | Speedup | Accuracy |
|-----------|----------|------------------|---------|----------|
| 100 genes × 3000 cells | 4.4s | 2.9s | 1.5x | 0.9999 |
| 500 genes × 2000 cells | 133s | 34s | **4.0x** | 1.0000 |
| 1000 genes × 2500 cells | 517s | 150s | **3.5x** | 1.0000 |

Tested on Apple Silicon (M-series) with 10 cores.

## Installation

```r
# Install from GitHub
library(remotes)
install_github('Zaoqu-Liu/scTenifoldKnk')
```

## Quick Start

```r
library(scTenifoldKnk)

# Load example data
scRNAseq <- system.file("single-cell/example.csv", package="scTenifoldKnk")
scRNAseq <- read.csv(scRNAseq, row.names = 1)

# Run analysis
result <- scTenifoldKnk(
  countMatrix = scRNAseq, 
  gKO = 'G100',
  qc_minLSize = 0
)

# View results
head(result$diffRegulation)
```

## Usage

### Basic Analysis

```r
result <- scTenifoldKnk(
  countMatrix = your_data,  # Gene x Cell matrix (raw counts)
  gKO = 'GENE_NAME',        # Gene to knockout
  qc_minLSize = 1000        # Minimum library size
)
```

### Adjust Parameters

```r
# Fast mode (for exploration)
result <- scTenifoldKnk(
  countMatrix = data,
  gKO = 'GENE',
  nc_nNet = 5,              # Fewer networks
  nc_nCells = 300           # Fewer cells per network
)

# Standard mode (default, balanced)
result <- scTenifoldKnk(
  countMatrix = data,
  gKO = 'GENE'
)

# Robust mode (more accurate, slower)
result <- scTenifoldKnk(
  countMatrix = data,
  gKO = 'GENE',
  nc_nNet = 20,             # More networks
  nc_nCells = 700           # More cells per network
)
```

### Control Parallelization

```r
# Auto-detect cores (default)
result <- scTenifoldKnk(data, gKO = 'GENE')

# Specify cores
result <- scTenifoldKnk(data, gKO = 'GENE', nCores = 4)

# Sequential processing
result <- scTenifoldKnk(data, gKO = 'GENE', nCores = 1)
```

## Output

Returns a list with three components:

1. **tensorNetworks**: Wild-type (WT) and knockout (KO) network matrices
2. **manifoldAlignment**: Low-dimensional gene embeddings
3. **diffRegulation**: Differential regulation statistics
   - `gene`: Gene names
   - `distance`: Euclidean distance between WT and KO
   - `Z`: Z-score
   - `FC`: Fold change
   - `p.value`: P-value
   - `p.adj`: FDR-adjusted p-value

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gKO` | - | Gene to knockout (required) |
| `qc_mtThreshold` | 0.1 | Max mitochondrial read ratio |
| `qc_minLSize` | 1000 | Min library size |
| `nc_nNet` | 10 | Number of networks |
| `nc_nCells` | 500 | Cells per network |
| `nc_nComp` | 3 | Number of PCA components |
| `nc_q` | 0.9 | Edge filtering quantile |
| `nCores` | NULL | CPU cores (NULL = auto-detect) |

## Requirements

- R ≥ 3.6
- Dependencies: Matrix, RSpectra, MASS, parallel

## Citation

**Original Method**:  
Osorio D, et al. scTenifoldKnk: An Efficient Virtual Knockout Tool for Gene Function Prediction via Single-Cell Gene Regulatory Network Perturbation. Patterns. 2020.

**Optimized Version**:  
Liu Z. Performance Optimization of scTenifoldKnk. GitHub: https://github.com/Zaoqu-Liu/scTenifoldKnk, 2024.

## License

GPL (≥ 2)

## Credits

**Original Authors**: Daniel Osorio, James Cai, and contributors from Texas A&M University

**Optimization**: Zaoqu Liu (ORCID: 0000-0002-0452-742X)

## Links

- Original Repository: https://github.com/cailab-tamu/scTenifoldKnk
- Optimized Repository: https://github.com/Zaoqu-Liu/scTenifoldKnk
- Issues: https://github.com/Zaoqu-Liu/scTenifoldKnk/issues
