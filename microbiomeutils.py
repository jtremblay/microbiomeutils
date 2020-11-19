#!/usr/bin/env python
 
"""A Python script to generate beta-diversity distance matrix + pcoa
and taxonomic summary from a feature table. The stated goal of this utility
is to avoid using .biom tables whose usage was enforced in later QIIME releases.
I ended up realizing that myself and collaborators pretty much never use these 
.biom tables anyways. So why bother generating them in the first place?

Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
"""
 
import os
import sys
import argparse
import re

cwd = os.getcwd()
sys.path.append(cwd + '/nrc')
from nrc.utils import *
from skbio.diversity import beta_diversity
from skbio import TreeNode
from emperor import Emperor
from skbio import OrdinationResults
from numpy import array, concatenate, repeat, zeros, nan
import pandas as pd
import textwrap

def betadiv(infile_feature_table, metric="bc", infile_tree=False):

    sample_names, feature_ids, data, lineages = getOtuTable(infile_feature_table)
    data = array(data).transpose()

    if metric == "bray-curtis":
        bc_dm = beta_diversity("braycurtis", data, sample_names)
        #print(bc_dm)
        bc_dm.write(sys.stdout, delimiter='\t')

    elif metric == "weighted-unifrac":
        # For wuf, include tree. put validate=False in beta_diversity() because even if tree is rooted, returns that tree is not rooted...
        tree = TreeNode.read(infile_tree)
        wuf_dm = beta_diversity("weighted_unifrac", data, sample_names, tree=tree, otu_ids=feature_ids, validate=False) 
        #print(wuf_dm)
        wuf_dm.write(sys.stdout, delimiter='\t')
    else:
        raise ValueError('Only support for bray-curtis and weighted-unifrac metrics. ' + 'Metric ' + metric + 'is not supported.')

def taxsum(infile_feature_table, sumtype, level):
    # Then perform taxonomic summary open(feature_table_fp, 'U')
    if sumtype == "absolute":
        feature_table = parse_feature_table(open(infile_feature_table, 'r'))
        #feature_table_rel = convert_feature_table_relative(feature_table)
        summary, header = make_summary(feature_table, int(level), None, None)
        write_summarize_taxa(summary, header, sys.stdout, delimiter=";", transposed_output=False)
    elif sumtype == "relative":
        feature_table = parse_feature_table(open(infile_feature_table, 'r'))
        feature_table_rel = convert_feature_table_relative(feature_table)
        summary, header = make_summary(feature_table_rel, int(level), None, None)
        write_summarize_taxa(summary, header, sys.stdout, delimiter=";", transposed_output=False)

def main(arguments):
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''additional information: output of taxsum, betadiv and pcoa are written in STDOUT
If you use microbiomeutils in your work, please cite:

    Tremblay, Julien
    microbiomeutils 0.9 : Microbiome utilities
    https://github.com/jtremblay/microbiomeutils
    
Thank you.'''))
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest="command")
    parser_bd = subparsers.add_parser('betadiv')
    parser_bd.add_argument('-i', '--infile-feature-table', help="Input file", type=argparse.FileType('r'))
    parser_bd.add_argument("-m", "--metric", help="Diversity metric (default: bray-curtis)", choices=["bray-curtis", "weighted-unifrac"], default="bray-curtis")
    parser_bd.add_argument("-t", "--infile-tree", help="Tree file (for weighted uniFrac)", type=argparse.FileType('r'))
                                                                                     
    #parser_bd.set_defaults(func=betadiv)
    parser_ts = subparsers.add_parser('taxsum')
    parser_ts.add_argument('-i', '--infile-feature-table', help="Input file", type=argparse.FileType('r'))
    parser_ts.add_argument("-t", "--sumtype", help="Summary type (default: absolute)", choices=["absolute", "relative"], default="absolute")
    parser_ts.add_argument("-l", "--level", help="Level <int> 1 to 7", choices=["1","2","3","4","5","6","7", "8"], default="3")
    #parser_bd.set_defaults(func=taxsum)
    parser_ts = subparsers.add_parser('pcoa')
    parser_ts.add_argument('-i', '--infile-distance-matrix', help="Input file", type=argparse.FileType('r'))
    
    parser_ts = subparsers.add_parser('emperor')
    parser_ts.add_argument('-i', '--infile-coords', help="Input file", type=argparse.FileType('r'))
    parser_ts.add_argument('-m', '--mapping-file', help="Mapping file", type=argparse.FileType('r'))
    parser_ts.add_argument('-o', '--outdir', help="Output directory")
    
    args = parser.parse_args(arguments)
    
    if args.command == 'betadiv':
        infile_feature_table = os.path.abspath(args.infile_feature_table.name)
        sys.stderr.write("[betadiv]\n")
        if args.infile_tree is None and args.metric == "weighted-unifrac":
            raise ValueError('weighted-unifrac needs a tree supplied. --infile-tree needed')
        
        if args.metric == "bray-curtis":
            betadiv(infile_feature_table, args.metric)
        else:
            betadiv(infile_feature_table, args.metric, args.infile_tree.name)
    
    elif args.command == 'taxsum':
        infile_feature_table = os.path.abspath(args.infile_feature_table.name)
        sys.stderr.write("[taxsum]\n")
        taxsum(infile_feature_table, args.sumtype, args.level)
    
    elif args.command == 'pcoa':
        sys.stderr.write("[pcoa]\n")
        infile_distance_matrix = os.path.abspath(args.infile_distance_matrix.name)
        ord_res = do_pcoa(infile_distance_matrix)
    
    elif args.command == 'emperor':
        sys.stderr.write("[emperor]\n")
        metadata = pd.read_csv(args.mapping_file, sep='\t', index_col='#SampleID')
        ordination = OrdinationResults.read(args.infile_coords)

        # the remote argument refers to where the support files will be located
        # relative to the plot itself i.e. index.html.
        emp = Emperor(ordination, metadata, remote='.')
        output_folder = args.outdir # new folder where data will be saved

        # create an output directory
        os.makedirs(output_folder, exist_ok=True)

        with open(os.path.join(output_folder, 'index.html'), 'w') as f:
            f.write(emp.make_emperor(standalone=True))
            emp.copy_support_files(output_folder)

     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


