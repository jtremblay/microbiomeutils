# MicrobiomeUtils

## Summary
A lightweight Python utility to compute beta-diversity distance matrices, PCoA ordinations,
interactive 3D PCoA plots (EMPeror) and taxonomic summaries directly from a tab-separated
feature table.

The goal is to avoid the `.biom` and `.qza` formats whose use was enforced in later QIIME
releases, and the complex QIIME installation that goes with them.

Alpha-diversity metrics are not implemented, as they are already efficiently computed by the
Rarefaction Tool Kit (RTK) package, which I recommend: https://github.com/hildebra/Rarefaction/

## Installation
Requires Python 3.11.

```sh
git clone https://github.com/jtremblay/microbiomeutils.git
cd microbiomeutils
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` pins the versions the tool was tested with (notably scikit-bio 0.5.9).
Installing newer versions of scikit-bio or emperor is not yet supported.

## Docker / Apptainer
A container image with all dependencies is published on the GitHub Container Registry as
`ghcr.io/jtremblay/microbiomeutils` (tags: `<version>`, e.g. `0.9.5`; `latest` = latest release;
`main` = current `main` branch). The image runs `microbiomeutils.py` directly, so pass the
subcommand and its options as arguments.

### Docker
The image's working directory is `/data`: mount the directory holding your files there and use
relative paths.

```sh
docker run --rm -v "$PWD:/data" ghcr.io/jtremblay/microbiomeutils:0.9.5 \
    betadiv -i feature_table.tsv -m bray-curtis > bc_res.tsv
```

For `emperor`, also add `-u "$(id -u):$(id -g)"` so that the output folder is created with your
user ID.

### Apptainer (e.g. Alliance / Digital Research Alliance of Canada clusters)
Build the `.sif` image once, inside an interactive job (`salloc`) rather than on a login node, and
store it in `/project` (not in `$SLURM_TMPDIR`, which is deleted at the end of the job):

```sh
module load apptainer
export APPTAINER_CACHEDIR=$SLURM_TMPDIR
apptainer build microbiomeutils_0.9.5.sif docker://ghcr.io/jtremblay/microbiomeutils:0.9.5
```

Then run it:

```sh
apptainer run -C -B /project -B /scratch -W $SLURM_TMPDIR microbiomeutils_0.9.5.sif \
    taxsum -i feature_table.tsv -l 6 > taxonomy_L6.tsv
```

With Apptainer, the working directory is your current directory on the host (not `/data`), so
relative paths refer to it. Run the command from a directory under `/project` or `/scratch`: with
`-C`, host directories (including `$HOME`) are only visible in the container if bound with `-B`.

## Usage
```sh
./microbiomeutils.py -h
./microbiomeutils.py <subcommand> -h
```

| Subcommand | Input | Output |
|---|---|---|
| `betadiv` | feature table (+ tree for UniFrac) | distance matrix |
| `pcoa` | distance matrix | PCoA ordination results |
| `emperor` | PCoA ordination results + mapping file | interactive 3D plot (HTML folder) |
| `taxsum` | feature table with taxonomy | taxonomic summary table |

`betadiv`, `pcoa` and `taxsum` write their results to standard output; redirect it to a file.
Progress messages are written to standard error.

## Input formats

### Feature table
A tab-separated file with one feature (OTU, ASV, gene, ...) per row and one sample per column.
See `./data/feature_table.tsv` for an example.

- **Header line:** the first non-comment line holds the sample names. It may start with
  `#FEATURE_ID`, `#FEATURE ID`, `#FEATUREID`, `#FEATURE`, `#OTU ID`, `contig_id` or `gene_id`.
- **Comment lines:** lines starting with `#` before the header (for example QIIME's
  `# Constructed from biom file`) are ignored, as are blank lines.
- **Counts:** integers or decimal numbers.
- **Taxonomy (optional, required by `taxsum`):** the last column, named `taxonomy`,
  `Consensus Lineage` or `OTU Metadata`, with ranks separated by `;`
  (e.g. `k__Bacteria;p__Firmicutes;c__Bacilli`).

```
#FEATURE_ID	S1	S2	S3	taxonomy
f1	10	0	3	k__Bacteria;p__Firmicutes;c__Bacilli
f2	4	7	0	k__Bacteria;p__Proteobacteria
```

### Tree (UniFrac only)
A Newick tree whose tip labels match the feature IDs of the feature table
(see `./data/tree.fasttree`).

### Mapping file (`emperor` only)
A tab-separated sample metadata file with a `#SampleID` column holding the sample names
(see `./data/mapping_file.tsv`). Its other columns can be used to color the plot.

## Examples

### Bray-Curtis beta diversity
```sh
./microbiomeutils.py betadiv -i ./data/feature_table.tsv -m bray-curtis > ./data/bc_res.tsv
./microbiomeutils.py pcoa -i ./data/bc_res.tsv > ./data/bc_res_coords.tsv
```

### Weighted and unweighted UniFrac beta diversity
```sh
./microbiomeutils.py betadiv -i ./data/feature_table.tsv -m weighted-unifrac \
    --infile-tree ./data/tree.fasttree > ./data/wuf_res.tsv
./microbiomeutils.py pcoa -i ./data/wuf_res.tsv > ./data/wuf_res_coords.tsv

./microbiomeutils.py betadiv -i ./data/feature_table.tsv -m unweighted-unifrac \
    --infile-tree ./data/tree.fasttree > ./data/uwuf_res.tsv
./microbiomeutils.py pcoa -i ./data/uwuf_res.tsv > ./data/uwuf_res_coords.tsv
```

### Interactive 3D PCoA plot with EMPeror
```sh
./microbiomeutils.py emperor -i ./data/wuf_res_coords.tsv -m ./data/mapping_file.tsv \
    -o ./data/weighted_unifrac_3d_plot
```
Open `./data/weighted_unifrac_3d_plot/index.html` in a web browser.

### Taxonomic summaries
`-l` sets the taxonomic level (1 = first rank, e.g. kingdom; default: 3).
Lineages shorter than the requested level are padded with `Other`.
`-t relative` gives relative abundances (each sample sums to 1) instead of absolute counts
(default: `absolute`).

```sh
for level in 1 2 3 4 5 6 7; do
    ./microbiomeutils.py taxsum -i ./data/feature_table.tsv -l $level > ./data/taxonomy_L$level.tsv
done

./microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 6 -t relative > ./data/taxonomy_L6_relative.tsv
```

## Running the tests
```sh
pip install -r requirements-dev.txt
python3 -m pytest tests/
```

## License
MIT; see [LICENSE](LICENSE).

## Citation
If you use microbiomeutils in your work, please cite:

Tremblay, Julien. microbiomeutils 0.9.5: Microbiome utilities.
https://github.com/jtremblay/microbiomeutils
