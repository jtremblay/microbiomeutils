__author__ = "Julien Tremblay"
__copyright__ = "Copyright 2026" 
__credits__ = ["Julien Tremblay"]
__license__ = "MIT"
__version__ = "0.9.6"
__maintainer__ = "Julien Tremblay"
__email__ = "jtremblay514@gmail.com"
__status__ = "Release"


"""Contains general utility code in support of nrc_tools python scripts.
Code included in this file was mostly taken from the old Qiime 1.4 release and adapted
for Python 3.9.0. The objective of creating these utilities was to avoid using the .biom
file format which is not needed in many case.
"""
import sys
import os
from io import StringIO
from numpy import array, concatenate, repeat, zeros, nan, asarray, where, isnan, float64, around
from re import compile, sub
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix


class FeatureMissingError(Exception):
    pass


def parse_feature_table(lines,count_map_f=int):
    """parses feature table (sample ID x feature ID map)

    Returns tuple: sample_ids, feature_ids, matrix of features(rows) x samples(cols),
    and lineages from infile.
    """
    feature_table = []
    feature_ids = []
    metadata = []
    sample_ids = []
    # iterate over lines in the feature table -- keep track of line number 
    # to support legacy (Qiime 1.2.0 and earlier) feature tables
    for i, line in enumerate(lines):
        
        if line.strip():
            if line.startswith('\t') and not sample_ids:
                line = line.lstrip()
                try:
                    sample_ids, has_metadata = process_feature_table_sample_ids(line.strip().split('\t')[0:])
                except ValueError:
                    raise ValueError(
                     "Error parsing sample IDs in feature table. Appears to be a"+\
                     " legacy feature table. Sample ID line:\n %s" % line)
                
                sys.stderr.write("sample_ids: " + str(sample_ids) + "\n")

            elif line.startswith(('#FEATURE_ID', '#OTU ID', '#FEATURE ID', '#FEATURE', '#FEATUREID', 'contig_id', 'gene_id')) and not sample_ids:
                # we've got a legacy feature table
                try:
                    sample_ids, has_metadata = process_feature_table_sample_ids(line.strip().split('\t')[1:])
                except ValueError:
                    raise ValueError(
                     "Error parsing sample IDs in feature table. Appears to be a"+\
                     " legacy feature table. Sample ID line:\n %s" % line)
            elif not line.startswith('#'):
                if not sample_ids:
                    # current line is the first non-space, non-comment line 
                    # in feature table, so contains the sample IDs
                    try:
                        sample_ids, has_metadata = process_feature_table_sample_ids(
                         line.strip().split('\t')[1:])
                    except ValueError:
                        raise ValueError(
                         "Error parsing sample IDs in feature table."+\
                         " Sample ID line:\n %s" % line)
                else:
                    # current line is feature line in feature table
                    fields = line.split('\t')
                    # grab the feature ID
                    feature_id = fields[0].strip()
                    feature_ids.append(feature_id)
                    if has_metadata:
                        # if there is feature metadata the last column gets appended
                        # to the metadata list
                        # added in a try/except to handle feature tables containing
                        # floating numbers
                        try:
                            feature_table.append(array(list(map(count_map_f, fields[1:-1]))))
                        except ValueError:
                            feature_table.append(array(list(map(float, fields[1:-1]))))
                            
                        metadata.append([x.strip() for x in fields[-1].split(';')])
                    else:
                        # otherwise all columns are appended to feature_table
                        # added in a try/except to handle feature tables containing
                        # floating numbers
                        try:
                            feature_table.append(array(list(map(count_map_f,fields[1:]))))
                        except ValueError:
                            feature_table.append(array(list(map(float, fields[1:]))))
                        
    return sample_ids, feature_ids, (feature_table), metadata


def getFeatureTable(feature_source):
    """Returns parsed feature table from putative feature source."""
    
    #if we have a string starting with #, assume it's a feature file,
    #otherwise assume it's a path
    # if 4-tuple, just return it
    if type(feature_source) == type((1,3,4,44)):
        return feature_source
    if hasattr(feature_source, 'startswith') and feature_source.startswith('#'):
        try:
            return parse_feature_table(StringIO(feature_source))
        except (TypeError, ValueError) as e:
            raise FeatureMissingError(
                "Tried to read features from string starting with # but got " + str(e))
    else:
        try:
            feature_file = open(feature_source, 'r')
        except (TypeError, IOError):
            raise FeatureMissingError(
                "Couldn't read FEATURE table file at path: %s" % feature_source)
        result = parse_feature_table(feature_file)
        feature_file.close()
        return result


def process_feature_table_sample_ids(sample_id_fields):
    """ process the sample IDs line of a feature table """
    if len(sample_id_fields) == 0:
            raise ValueError(
             'Error parsing sample ID line in feature table. Fields are %s' \
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

# Taxonomic summary
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
    for consensus, feature_counts in sorted(counts_by_consensus.items()):
        if lower_percentage!=None and \
                                feature_counts.sum()/total_counts>lower_percentage:
            continue
        elif upper_percentage!=None and \
                                feature_counts.sum()/total_counts<upper_percentage:
            continue
        new_row = [(consensus)]
        new_row.extend(feature_counts)
        taxonomy_summary.append(new_row)

    return taxonomy_summary, header

def sum_counts_by_consensus(feature_table, level, missing_name='Other'):
    """Returns a dict keyed by consensus, valued by feature counts

    feature counts are summed together if they have the same consensus

    if the consensus string doesn't reach to level, missing_name is appended on
    until the taxonomy string is of length level
    """
    result = {}
    sample_map = dict([(s,i) for i,s in enumerate(feature_table[0])])
    for counts, consensus in zip(feature_table[2], feature_table[3]):
        counts = counts.astype(float)
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

def write_summarize_taxa(summary, header, output_fp, delimiter=';', transposed_output=False):
    """ """
    if transposed_output:
         # transposing the summary
         summary = [[r[col] for r in summary] for col in range(len(summary[0]))]
         # adding the first column into the new summary matrix
         for i in range(1,len(summary)):
             summary[i] = [([header[i]])] + summary[i]
         # replacing header and trimming summary
         header = ['SampleID'] + [delimiter.join(taxon) for taxon in summary[0]]
         summary = summary[1:]
    
    for line in format_summarize_taxa(summary, header, delimiter):
        sys.stdout.write(line)

def format_summarize_taxa(summary, header, delimiter=';'):
    """Formats a summarized taxonomy table for output"""
    yield "%s\n" % '\t'.join(header)
    for row in summary:
        # taxon is tuple, join together for foo;bar;foobar
        taxon = row[0]
        line = [delimiter.join(taxon)]

        # add on feature counts
        line.extend(map(str, around(row[1:], 3)))

        yield "%s\n" % '\t'.join(line)

def convert_feature_table_relative(feature_table):
    """Convert the feature table to relative abundances

    this method works on a parsed feature table
    """
    sample_ids, feature_ids, feature_counts, consensus = feature_table
    feature_counts = asarray(feature_counts, float)
    feature_counts = feature_counts / feature_counts.sum(axis=0)
    feature_counts = where(isnan(feature_counts), 0.0, feature_counts)
    return (sample_ids, feature_ids, feature_counts, consensus)

def do_pcoa(infile):
    """Run PCoA on a tab-separated (lsmat) distance matrix and write the
    ordination results to stdout.

    DistanceMatrix.read() checks that the matrix is square, symmetric and
    that row labels match the header in the same order.
    """
    distmtx = DistanceMatrix.read(infile, format='lsmat')
    ord_res = pcoa(distmtx)
    ord_res.write(sys.stdout)
