# scTenifoldKnk

A workflow for in-silico gene knockout experiments using single-cell RNA-seq data.

## Version 2.0.0

**Maintainer**: Zaoqu Liu (liuzaoqu@163.com)  
**Repository**: https://github.com/Zaoqu-Liu/scTenifoldKnk

### Updates in This Version

This version includes performance optimizations while maintaining compatibility with the original implementation.

- **Performance improvements** for medium to large datasets:
  - 500 genes: ~4x faster (benchmarked: 133s → 34s)
  - 1000 genes: ~3.5x faster (benchmarked: 517s → 150s)
  - Note: Small datasets (<200 genes) show minimal speedup due to parallelization overhead
  
- **Result consistency**: High correlation (>0.999) with original version on test datasets

- **Adaptive parallelization**: Automatically selects sequential or parallel mode based on data size

- **Implementation changes**:
  - Optimized matrix operations
  - Vectorized computations where possible
  - Multi-core processing for network construction

### Benchmark Results

Results from tests on Apple Silicon (M-series, 10 cores) using synthetic data:

| Data Size | v1.0.1 | v2.0.0 | Speedup | Correlation |
|-----------|--------|--------|---------|-------------|
| 100 genes × 3000 cells | 4.4s | 2.9s | 1.5x | 0.9999 |
| 500 genes × 2000 cells | 133s | 34s | 4.0x | 1.0000 |
| 1000 genes × 2500 cells | 517s | 150s | 3.5x | 1.0000 |

*Note: Performance may vary depending on hardware, data characteristics, and system load.*

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
- Dependencies: Matrix, RSpectra, MASS, parallel, Rcpp, RcppEigen
- C++ compiler for building from source

## Limitations

- Performance improvements are most noticeable for datasets with >500 genes
- Very small datasets may not benefit from parallelization
- Results may differ slightly from v1.0.1 due to different random subsampling (use same seed for reproducibility)

## Citation

**Original Method**:  
Osorio D, et al. scTenifoldKnk: An Efficient Virtual Knockout Tool for Gene Function Prediction via Single-Cell Gene Regulatory Network Perturbation. Patterns. 2020.

**This Version**:  
If you use this optimized version, please also cite:  
Liu Z. scTenifoldKnk v2.0.0 with performance improvements. GitHub: https://github.com/Zaoqu-Liu/scTenifoldKnk, 2024.

## License

GPL (≥ 2)

## Credits

**Original Authors**: Daniel Osorio, James Cai, and contributors from Texas A&M University

**Performance improvements**: Zaoqu Liu (ORCID: 0000-0002-0452-742X)

## Links

- Original Repository: https://github.com/cailab-tamu/scTenifoldKnk
- Optimized Repository: https://github.com/Zaoqu-Liu/scTenifoldKnk
- Issues: https://github.com/Zaoqu-Liu/scTenifoldKnk/issues
