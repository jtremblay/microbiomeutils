#!/usr/bin/env python

__author__ = "Julien Tremblay"
__copyright__ = "Copyright 2020, NRC tools" 
__credits__ = ["Julien Tremblay"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Julien Tremblay"
__email__ = "julien.tremblay@nrc-cnrc.gc.ca"
__status__ = "Release"


"""Contains general utility code in support of nrc_tools python scripts.
Code included in this file was mostly taken from the old Qiime 1.4 release and adapted
for Python 3.9.0. The objective of creating these utilities was to avoid using the .biom
file format which is not needed in many case.
"""
import sys
import os
from numpy import array, concatenate, repeat, zeros, nan, asarray, float, where, isnan, float64
from re import compile, sub
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix

def parse_feature_table(lines,count_map_f=int):
    """parses otu table (sample ID x OTU ID map)

    Returns tuple: sample_ids, otu_ids, matrix of OTUs(rows) x samples(cols),
    and lineages from infile.
    """
    feature_table = []
    otu_ids = []
    metadata = []
    sample_ids = []
    # iterate over lines in the OTU table -- keep track of line number 
    # to support legacy (Qiime 1.2.0 and earlier) OTU tables
    for i, line in enumerate(lines):
        line = line.strip()
        
        if line:
            if i == 0 and line.startswith(('#FEATURE_ID', '#OTU ID', '#FEATURE ID', '#FEATURE', '#FEATUREID')) and not sample_ids:
                # we've got a legacy OTU table
                try:
                    sample_ids, has_metadata = process_feature_table_sample_ids(line.strip().split('\t')[1:])
                except ValueError:
                    raise ValueError(
                     "Error parsing sample IDs in OTU table. Appears to be a"+\
                     " legacy OTU table. Sample ID line:\n %s" % line)
            elif not line.startswith('#'):
                if not sample_ids:
                    # current line is the first non-space, non-comment line 
                    # in OTU table, so contains the sample IDs
                    try:
                        sample_ids, has_metadata = process_feature_table_sample_ids(
                         line.strip().split('\t')[1:])
                    except ValueError:
                        raise ValueError(
                         "Error parsing sample IDs in OTU table."+\
                         " Sample ID line:\n %s" % line)
                else:
                    # current line is OTU line in OTU table
                    fields = line.split('\t')
                    # grab the OTU ID
                    otu_id = fields[0].strip()
                    otu_ids.append(otu_id)
                    if has_metadata:
                        # if there is OTU metadata the last column gets appended
                        # to the metadata list
                        # added in a try/except to handle OTU tables containing
                        # floating numbers
                        try:
                            feature_table.append(array(list(map(count_map_f, fields[1:-1]))))
                        except ValueError:
                            feature_table.append(array(list(map(float, fields[1:-1]))))
                            
                        metadata.append(map(str.strip, fields[-1].split(';')))
                    else:
                        # otherwise all columns are appended to feature_table
                        # added in a try/except to handle OTU tables containing
                        # floating numbers
                        try:
                            feature_table.append(array(list(map(count_map_f,fields[1:]))))
                        except ValueError:
                            feature_table.append(array(list(map(float, fields[1:]))))
                        
    return sample_ids, otu_ids, (feature_table), metadata


#####
def getOtuTable(otu_source):
    """Returns parsed OTU table from putative OTU source."""
    
    #if we have a string starting with #, assume it's an OTU file,
    #otherwise assume it's a path
    # if 4-tuple, just return it
    if type(otu_source) == type((1,3,4,44)):
        return otu_source
    if hasattr(otu_source, 'startswith') and otu_source.startswith('#'):
        try:
            return parse_feature_table(StringIO(otu_source))
        except (TypeError, ValueError) as e:
            raise OtuMissingError(
                "Tried to read OTUs from string starting with # but got "+e)
    else:
        try:
            otu_file = open(otu_source, 'U')
        except (TypeError, IOError):
            raise OtuMissingError(
                "Couldn't read FEATURE table file at path: %s" % otu_source)
        result = parse_feature_table(otu_file)
        otu_file.close()
        return result


def getTree(tree_source):
    """Returns parsed tree from putative tree source"""
    if isinstance(tree_source, PhyloNode):
        tree = tree_source    #accept tree object directly for tests
    elif tree_source:
        try:
            f = open(tree_source, 'U')
        except (TypeError, IOError):
            raise TreeMissingError(
                "Couldn't read tree file at path: %s" % tree_source)
        tree = parse_newick(f, PhyloNode)
        f.close()
    else:
        raise TreeMissingError( str(self.Name) + \
            " is a phylogenetic metric, but no tree was supplied.")
    return tree


def process_feature_table_sample_ids(sample_id_fields):
    """ process the sample IDs line of an OTU table """
    if len(sample_id_fields) == 0:
            raise ValueError(
             'Error parsing sample ID line in OTU table. Fields are %s' \
             % ' '.join(sample_id_fields))
            
    # Detect if a metadata column is included as the last column. This
    # field will be named either 'Consensus Lineage' or 'OTU Metadata',
    # but we don't care about case or spaces.
    last_column_header = sample_id_fields[-1].strip().replace(' ','').lower()
    if last_column_header in ['consensuslineage', 'otumetadata', 'taxonomy']:
        has_metadata = True
        sample_ids = sample_id_fields[:-1]
    else:
        has_metadata = False
        sample_ids = sample_id_fields
    
    # Return the list of sample IDs and boolean indicating if a metadata
    # column is included.
    return sample_ids, has_metadata

###
#Taxonomic summary
def make_summary(feature_table, level, upper_percentage, lower_percentage): 
    """Returns taxonomy summary data

    header is a list of:
    [(Taxon),sample1,sample2,...]

    taxonomy_summary is a list of lists of:
    [[(taxon1),count,count,...],[(taxon2),count,count,...]...]
    """
    header = ['Taxon']
    header.extend(feature_table[0]) # sample ids

    counts_by_consensus, sample_map = sum_counts_by_consensus(feature_table, level)

    total_counts = float(sum([sum(i) for i in counts_by_consensus.values()]))
    taxonomy_summary = []
    for consensus, otu_counts in sorted(counts_by_consensus.items()):
        if lower_percentage!=None and \
                                otu_counts.sum()/total_counts>lower_percentage:
            continue
        elif upper_percentage!=None and \
                                otu_counts.sum()/total_counts<upper_percentage:
            continue
        new_row = [(consensus)]
        new_row.extend(otu_counts)
        taxonomy_summary.append(new_row)

    return taxonomy_summary, header

def sum_counts_by_consensus(feature_table, level, missing_name='Other'):
    """Returns a dict keyed by consensus, valued by otu counts

    otu counts are summed together if they have the same consensus

    if the consensus string doesn't reach to level, missing_name is appended on
    until the taxonomy string is of length level
    """
    result = {}
    sample_map = dict([(s,i) for i,s in enumerate(feature_table[0])])
    for counts, consensus in zip(feature_table[2], feature_table[3]):
        consensus = list(consensus)
        n_ranks = len(consensus)
        if n_ranks > level:
            consensus = consensus[:level]
        elif n_ranks < level:
            consensus.extend([missing_name for i in range(level - n_ranks)])
        else:
            # consensus is the correct number of levels
            pass

        consensus = tuple(consensus)
        if consensus in result:
            result[consensus] += counts
        else:
            result[consensus] = counts.copy()

    return result, sample_map

#def add_summary_mapping(feature_table, mapping, level): 
#    """Returns sample summary of sample counts by taxon
#    
#    Summary is keyed by sample_id, valued by otu counts for each taxon
#    Taxon order is a list of taxons where idx n corresponds to otu count idx n
#    """
#    counts_by_consensus, sample_map = sum_counts_by_consensus(feature_table, level)
#    
#    summary = defaultdict(list)
#    for row in mapping:
#        # grab otu idx if the sample exists, otherwise ignore it
#        sample_id = row[0]
#        if sample_id not in sample_map:
#            continue
#        otu_idx = sample_map[sample_id]
#
#        for consensus, counts in sorted(counts_by_consensus.items()):
#            summary[sample_id].append(counts[otu_idx])
#
#    taxon_order = sorted(counts_by_consensus.keys())
#
#    return summary, taxon_order

def write_summarize_taxa(summary, header, output_fp, delimiter=';', transposed_output=False):
    """ """
    # Fixing headers
    pattern = compile('\W')
    header = [sub(pattern, '.', label) for label in header]

    if transposed_output:
         # transposing the summary
         summary = [[r[col] for r in summary] for col in range(len(summary[0]))]
         # adding the first column into the new summary matrix
         for i in range(1,len(summary)):
             summary[i] = [([header[i]])] + summary[i]
         # replacing header and trimming summary
         header = ['SampleID'] + [delimiter.join(taxon) for taxon in summary[0]]
         summary = summary[1:]
    
    #of = open(output_fp,'w')
    for line in format_summarize_taxa(summary, header, delimiter):
        #of.write(line)
        sys.stdout.write(line)
    #of.close()

def format_summarize_taxa(summary, header, delimiter=';'):
    """Formats a summarized taxonomy table for output"""
    yield "%s\n" % '\t'.join(header)
    for row in summary:
        # taxon is tuple, join together for foo;bar;foobar
        taxon = row[0]
        line = [delimiter.join(taxon)]

        # add on otu counts
        line.extend(map(str, row[1:]))

        yield "%s\n" % '\t'.join(line)

def convert_feature_table_relative(feature_table):
    """Convert the OTU table to relative abundances

    this method works on a parsed OTU table
    """
    sample_ids, otu_ids, otu_counts, consensus = feature_table
    otu_counts = asarray(otu_counts, float)
    otu_counts = otu_counts / otu_counts.sum(axis=0)
    otu_counts = where(isnan(otu_counts), 0.0, otu_counts)
    return (sample_ids, otu_ids, otu_counts, consensus)

def parse_distmat(infile):
    """Parser for distance matrix file (e.g. UniFrac dist matrix).

    The examples I have of this file are just sample x sample tab-delimited
    text, so easiest way to handle is just to convert into a numpy array
    plus a list of field names. TODO: make sure square matrix is provided.
    """

    header = None
    result = []
    lines = open(infile, 'r')
    for line in lines:
        if line[0] == '\t': #is header
            header = list(map(str.strip, line.split('\t')[1:]))
            #sys.stderr.write(str(header) + "\n");
        else:
            tmp = list(map(float, line.split('\t')[1:]))
            result.append(tmp)
    
    result = asarray(result, dtype=float)
    #print(result)
    return(header, result)

def do_pcoa(infile):
    samples, distmtx = parse_distmat(infile)
    # coords, each row is an axis
    distmtx = DistanceMatrix(distmtx, ids=samples)
    ord_res = pcoa(distmtx)
    coords = ord_res.samples
    eigvals = ord_res.eigvals
    pcnts = ord_res.proportion_explained
    
    #Write results to output
    ord_res.write(sys.stdout)

