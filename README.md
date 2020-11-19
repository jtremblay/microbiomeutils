# microbiomeutils

## Summary
A Python utility to generate beta-diversity distance matrix + pcoa
and taxonomic summary from a feature table. The stated goal of this utility
is to avoid using .biom tables whose usage was enforced in later QIIME releases.
I ended up realizing that myself and collaborators pretty much never use these 
.biom tables anyways. So why bother generating them in the first place? This utility also provide a great alternative to QIIME1 or QIIME2 without going through a complex installation process or dealing with .qza files.
Functions for alpha diversity metrics generation were not implemented in this utility because it
can be already efficiently acomplished using the Rarefaction Tool Kit (RTK) package - https://github.com/hildebra/Rarefaction/ - which I recommend.

## Installation
This utility was written in Python 3.9.0 using scikit-bio v0.5.6. Once Python 3.9.0 is installed, run ```pip install "scikit-bio==0.5.6"``` or ```pip install scikit-bio```. Also install emperor: ```pip install emperor```.
Then clone this repo : ```git clone https://github.com/jtremblay/microbiomeutils.git``` and run ```microbiomeutils.py``` with the appropriate arguments as described below:

## Help
```microbiomeutils.py -h```

```microbiomeutils.py betadiv -h```

```microbiomeutils.py pcoa -h```

```microbiomeutils.py emperor -h```

```microbiomeutils.py taxsum -h```

## Data format
Takes in input a tab-separated value (tsv) feature table with just one header line that can start with : '#FEATURE_ID', '#OTU ID', '#FEATURE ID', '#FEATURE' or '#FEATUREID' string. Script won't work with multi-line header.
See example input file in ```./data/feature_table.tsv```

## Computing Bray-Curtis beta-diversity
 ```microbiomeutils.py betadiv -i ./data/feature_table.tsv -m bray-curtis > ./data/bc_res.tsv ```

 ```microbiomeutils.py pcoa -i ./data/bc_res.tsv > ./data/bc_res_coords.tsv ```

## Computing Weighted-Unifrac beta-diversity
 ```microbiomeutils.py betadiv -i ./data/feature_table.tsv -m weighted-unifrac --infile-tree ./data/tree.fasttree > ./data/wuf_res.tsv ```

 ```microbiomeutils.py pcoa -i ./data/wuf_res.tsv > ./data/wuf_res_coords.tsv ```

## Computing taxonomic summaries
 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 1 > data/taxonomy_L1.tsv ```

 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 2 > data/taxonomy_L2.tsv ```

 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 3 > data/taxonomy_L3.tsv ```

 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 4 > data/taxonomy_L4.tsv ```

 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 5 > data/taxonomy_L5.tsv ```

 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 6 > data/taxonomy_L6.tsv ```

 ```microbiomeutils.py taxsum -i ./data/feature_table.tsv -l 7 > data/taxonomy_L7.tsv ```

## Generating 3d PCoA plots with emperor
 ```microbiomeutils.py emperor -i ./data/wuf_res_coords.tsv -m ./data/mapping_file.tsv -o ./data/weighted_unifrac_3d_plot```

# Citation
If you use microbiomeutils in your work, please cite:

Tremblay, Julien

microbiomeutils 0.9 : Microbiome utilities

https://github.com/jtremblay/microbiomeutils

