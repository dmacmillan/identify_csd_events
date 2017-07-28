## Standard modules
import argparse
import os
import sys
import logging
import string
from functools import partial
from copy import deepcopy

## External modules
import pysam

## Functions
def reverse_complement(sequence, translation_table):
    return string.translate(sequence[::-1].upper(), translation_table)

# Search a string for any signals
def pas_search(string, signals, slen, mlen):
    result = {}
    string = string.upper()
    for i in xrange(slen - mlen + 1):
        substring = string[i:i + mlen]
        try:
            result[i] = signals[substring]
        except KeyError:
            continue
    return result

# Fetch a sequence from a pysam variant file
def fetch_seq(refseq, chrom, start, end):
    start = max(start, 0)
    try:
        end = min(end, refseq.get_reference_length(chrom))
    except KeyError:
        return None
    try:
        return refseq.fetch(chrom, start, end)
    except ValueError:
        return None

def get_alt_percent(variant, threshold, bc_lookup):
    num = variant.samples[0]['DP']
    den = variant.samples[1]['DP']
    result = round(100*(float(num) / den), 2)
    return result

# Check if jacusa variant is a pass or not
def check_jacusa_pass(variant):
    if variant.filter.keys() != ['*']:
        return False
    return True

def get_all_substrings(string, mlength):
    return [[i, string[i:mlength+i].upper()] for i in xrange(len(string) - mlength + 1)]

# Return true or false if the variant samples meet the coverage threshold
# threshold = int = minimum coverage threshold
def check_coverage(variant, threshold):
    for i in xrange(len(variant.samples)):
        if (variant.samples[i]['DP'] < threshold):
            return False
    return True

# Check if a variant matches one in dbsnp
# and return true if match otherwise return false
# dbsnp = pysam.cbcf.VariantFile
def check_dbsnp(variant, dbsnp):
    attrs = ('chrom', 'pos', 'start', 'stop')
    matches = dbsnp.fetch(variant.chrom, variant.start, variant.stop)
    for match in matches:
        if all([getattr(variant, attr) == getattr(match, attr) for attr in attrs]):
            return True
    return False

def classify(ref, alt, motifs, event_codes):
    # Either a shift_up, shift_down, or destroy event
    if (ref in motifs):
        logging.debug('ref {} in motifs'.format(ref))
        # Either shift_up or shift_down
        if (alt in motifs):
            logging.debug('alt {} in motifs'.format(alt))
            # Ref is stronger meaning shift_down event
            if (motifs[ref] < motifs[alt]):
                logging.debug('ref {} < alt {} shift_down'.format(ref, alt))
                return event_codes['shift_down']
            # Alt is stronger meaning shift_up event
            elif (motifs[alt] < motifs[ref]):
                logging.debug('alt {} < ref {}, shift_up'.format(alt, ref))
                return event_codes['shift_up']
        # Destroy event
        else:
            logging.debug('alt {} NOT in motifs, destroy'.format(alt))
            return event_codes['destroy']
    # Either a creation or no event
    else:
        logging.debug('ref {} NOT in motifs'.format(ref))
        # Creation event
        if (alt in motifs):
            logging.debug('alt {} in motifs, create'.format(alt))
            return event_codes['create']
        else:
            logging.debug('alt {} NOT in motifs, NO EVENT'.format(alt))
            return None

def change_character(text, index, char):
    return text[:index] + char + text[index + 1:]

def identify(variant, ref, motifs, event_codes, window=25, motif_len=6):
    results = {
        'is_only_motif': True,
        'ambiguous': False,
        'no_pas': False,
        'events': None
    }
    # Fetch initial sequence with where mutation could affect PAS
    seq = fetch_seq(
        ref,
        variant.chrom,
        variant.start - (motif_len - 1) - window,
        variant.stop + (motif_len - 1) + window
    )
    logging.debug('seq: {}'.format(seq))
    pas_refs = get_all_substrings(seq, motif_len)
    logging.debug('pas_refs left window: {}'.format(pas_refs[:window]))
    logging.debug('pas_refs regular: {}'.format(pas_refs[window : window + motif_len]))
    logging.debug('pas_refs right window: {}'.format(pas_refs[window + motif_len:]))
    # Copy the reference PAS
    pas_alts = deepcopy(pas_refs)
    # Apply the mutation to each reference PAS
    for i in xrange(window, window + motif_len):
        logging.debug('before: {}'.format(pas_alts[i][1]))
        logging.debug('Changing: change_character({}, {}, {})'.format(pas_alts[i][1], (motif_len + window - 1) - i, variant.alts[0]))
        pas_alts[i][1] = change_character(pas_alts[i][1], (motif_len + window - 1) - i, variant.alts[0])
        logging.debug('after: {}'.format(pas_alts[i][1]))
    logging.debug('pas_alts left window: {}'.format(pas_alts[:window]))
    logging.debug('pas_alts regular: {}'.format(pas_alts[window : window + motif_len]))
    logging.debug('pas_alts right window: {}'.format(pas_alts[window + motif_len:]))
    # If there are no PAS in either list then there are no events
    if not any([x[1] in motifs for x in pas_refs[window : window + motif_len] + pas_alts[window : window + motif_len]]):
        logging.debug('No PAS')
        results['no_pas'] = True
        return results
    # Classify events if any
    events = []
    logging.debug('motifs: {}'.format(motifs))
    for i in xrange(window, window + motif_len):
        logging.debug('Classify {} -> {}'.format(pas_refs[i][1], pas_alts[i][1]))
        event = classify(pas_refs[i][1], pas_alts[i][1], motifs, event_codes)
        if event is not None:
            logging.debug('event: {}'.format(event))
            # Convert local coordinate to genomic
            j = variant.pos - motif_len - window + 1 + i
            events.append({
                'pos': j,
                'event_code': event,
                'pas_ref': pas_refs[i][1],
                'pas_alt': pas_alts[i][1]
            })
    logging.debug('events: {}'.format(events))
    # Multiple events causes ambiguity
    nevents = len(events)
    if nevents > 1:
        logging.debug('multiple events: {}'.format(events))
        results['ambiguous'] = True
        results['events'] = events
        return results
    # If no events return nothing
    elif nevents == 0:
        return results
    # Check window for additional motifs
    outer_indices = range(window) + range(motif_len + window, motif_len + (window * 2) - 1)
    logging.debug('outer_indices: {}'.format(outer_indices))
    if any([pas_refs[x][1] in motifs for x in outer_indices]):
        results['is_only_motif'] = False
    results['events'] = events
    return results

def check_alts(variant):
    if len(variant.alts) > 1:
        return False
    return True

def get_strength(candidate, motifs):
    try:
        return motifs[candidate]
    except KeyError:
        return 'NAN'

def write_result(result, output_columns, delim='\t'):
    return (delim).join([
        str(result[x]) for x in output_columns
    ]) + '\n'

## Argument parser
parser = argparse.ArgumentParser(
    description = 'Given VCF output return a file displaying the ' \
    'create, shift, and destroy events for polyadenylation signals'
)

parser.add_argument(
    'variants',
    help = 'VCF output, must be tabix-indexed'
)
parser.add_argument(
    'cohort',
    help = 'Name of the cohort'
)
parser.add_argument(
    'pid',
    help = 'The patient ID'
)
parser.add_argument(
    'vcf_type',
    choices = (
        'jacusa',
        'strelka'
    ),
    help = 'Type of VCF output'
)
parser.add_argument(
    '-psw',
    '--pas_search_window',
    type = int,
    default = 25,
    help = 'The window surrounding each variant to search for additional polyadenylation signals ' \
    'in the patient\'s normal DNA. Default = 50'
)
parser.add_argument(
    '-reg',
    '--regions',
    default = '/projects/dmacmillanprj2/notebooks/PAS_SNP_Analysis/utr3.75.gff',
    help = 'A GFF file for which only variants within defined regions will be considered. ' \
    'Default = "/projects/dmacmillanprj2/notebooks/PAS_SNP_Analysis/utr3.75.gff"'
)
parser.add_argument(
    '-ref',
    '--reference',
    default = '/projects/dmacmillanprj2/notebooks/PAS_SNP_Analysis/patrick_homo_sapiens_assembly19.fasta',
    help = 'Indexed reference genome in FASTA format. Default = "/projects/dmacmillanprj2/notebooks/' \
    'PAS_SNP_Analysis/patrick_homo_sapiens_assembly19.fasta"'
)
parser.add_argument(
    '-d',
    '--delimiter',
    default = '\t',
    help = 'Delimiter to use in output. Default = "\\t"'
)
parser.add_argument(
    '-ev',
    '--exclude_variants',
    default = '/projects/dmacmillanprj2/notebooks/PAS_SNP_Analysis/common_all_20170403.vcf.gz',
    help = 'An indexed VCF file containing known variants to exclude from consideration'
)
parser.add_argument(
    '-mc',
    '--minimum_coverage',
    type = int,
    default = 10,
    help = 'The minimum level of coverage for a variant. Default = 10'
)
parser.add_argument(
    '-map',
    '--min_alt_percent',
    type = int,
    default = 20,
    help = 'For Jacusa output, the minimum percentage of alternative base ' \
    'present in the RNA-Seq pileup. For example setting this value to 60 ' \
    'would only accept variant calls if at least 60%% of bases at that position are ' \
    'one of the alternative bases. Default = 20'
)
parser.add_argument(
    '-l',
    '--log_level',
    default = 'warning',
    choices = [
        'debug',
        'info',
        'warning',
        'error',
        'critical'
    ],
    help = 'Set the logging level. Default = "warning"'
)
parser.add_argument(
    '-n',
    '--name',
    default = 'analysis',
    help = 'Name of output file. Default = "analysis"'
)
parser.add_argument(
    '-o',
    '--outdir',
    default = os.getcwd(),
    help = 'Path to output to. Default = "{}"'.format(os.getcwd())
)

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    try:
        os.makedirs(args.outdir)
    except OSError:
        pass

# Logging
logging.basicConfig(
    filename=os.path.join(args.outdir, '{}.log'.format(args.name)),
    level=getattr(logging, args.log_level.upper()),
    filemode='w',
    format='%(asctime)s %(message)s'
)

## Global variables
# Values indicate strength, 1 being the strongest
CANDIDATE_PASS = {
    'AATAAA': 1,
    'ATTAAA': 2,
    'AGTAAA': 3,
    'TATAAA': 4,
    'CATAAA': 5,
    'GATAAA': 6,
    'AATATA': 7,
    'AATACA': 8,
    'AATAGA': 9,
    'AAAAAG': 10,
    'ACTAAA': 11,
    'AAGAAA': 12,
    'AATGAA': 13,
    'TTTAAA': 14,
    'AAAACA': 15,
    'GGGGCT': 16
}

# Define motif length
motif_len = 6

# Strength lookup for PAS
REVERSE_CANDIDATE_PASS = {v:k for k,v in CANDIDATE_PASS.items()}

# Translation table for reverse complementing
translation_table = string.maketrans('TCGA','AGCT')

# Reverse complement function
rev_comp = partial(reverse_complement, translation_table = translation_table)

# Reverse complement lookup for PAS
REV_COMP_CANDIDATE_PASS = {rev_comp(k):v for k,v in CANDIDATE_PASS.items()}

# Output columns
output_columns = (
    'chrom',
    'pos',
    'gene',
    'strand',
    'cohort',
    'pid',
    'ref',
    'alt',
    'effect',
    'pas_ref',
    'pas_alt',
    'pas_ref_strength',
    'pas_alt_strength',
    'pas_start',
    'pas_end',
    'wgs_tumor_ref_depth',
    'wgs_tumor_alt_depth',
    'rna_tumor_ref_depth',
    'rna_tumor_alt_depth',
    'pass',
    'ambiguous',
    'interfered',
    'is_only_motif',
    'previously_known'
)

# Define the regions in which to consider variants
regions_file = open(args.regions, 'r')
regions_iterator = pysam.tabix_iterator(regions_file, parser=pysam.asGTF())

# Define variants
variants = pysam.VariantFile(args.variants)

# BC bases
bc_lookup = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}

# Event codes
event_codes = {
    'create': 0,
    'shift_up': 1,
    'shift_down': 2,
    'destroy': 3
}

# Reverse code lookup
codes_event = {v:k for k,v in event_codes.items()}

# Load reference sequence
reference = pysam.FastaFile(args.reference)

# Keep track of variants that did not pass all filters
not_passed = 0
# Keep track of variants that did not have enough coverage
low_coverage = 0
# Track the number of variants with more than one alternate base
multiple_alts = 0
# Track the number of variants with no PAS signals
no_pas = 0
# Track the number of variants with no events
no_event = 0
# Track the number of variants that are ambiguous
# Ambiguity comes from a single variant causing multiple effects
ambiguous = 0
# Track the number of variants that are interfered
# Interference comes from multiple variations within the same motif
interfere = 0
# Track the number of variants that are not the only motif
is_only_motif = 0
# Track the number of variants that match dbsnp calls
previously_known = 0
# Track the number of variants within regions
in_regions = 0


# Total number of variants
total_variants = len([x for x in pysam.VariantFile(args.variants)])

# If the input is not jacusa there is only one sample
sample_index = 0
# Pick index based on VCF input
if args.vcf_type == 'jacusa':
    sample_index = 1

# Define output file
outpath = os.path.join(args.outdir, '{}.info'.format(args.name))
outfile = open(outpath, 'w')
outfile.write((args.delimiter).join(output_columns) + '\n')

# Load known variants
known_variants = pysam.VariantFile(args.exclude_variants)

# Track variants which have been analyzed
analyzed = set()

last_result = None
# Iterate through defined regions:
iterator = regions_iterator
for region in iterator:
    gene = region.attributes.split('::')[1]
    motifs = CANDIDATE_PASS
    if region.strand == '-':
        motifs = REV_COMP_CANDIDATE_PASS
    subvariants = variants.fetch(region.contig, region.start, region.end)
    if not subvariants:
        continue
    for variant in subvariants:
        logging.debug('variant: {}'.format(str(variant).strip()))
        filters = {
            'pass' : True,
            'low_coverage' : False,
            'interfered' : False,
            'ambiguous' : False,
            'is_only_motif' : False,
            'previously_known' : False
        }
        identity = (variant.chrom, variant.pos, variant.ref, variant.alts)
        if identity in analyzed:
            continue
        else:
            analyzed.add(identity)
        in_regions += 1
        if not check_jacusa_pass(variant):
            # Variant did not pass
            not_passed += 1
            #continue
            filters['pass'] = False
        if check_dbsnp(variant, known_variants):
            previously_known += 1
            #continue
            filters['previously_known'] = True
        if not check_coverage(variant, args.minimum_coverage):
            # Variant does not have sufficient coverage
            low_coverage += 1
            #continue
            filters['low_coverage'] = True
        if not check_alts(variant):
            multiple_alts += 1
            continue
        result = identify(variant, reference, motifs, event_codes, window=args.pas_search_window)
        if result['no_pas']:
            no_pas += 1
            continue
        elif result['ambiguous']:
            ambiguous += 1
            #continue
            filters['ambiguous'] = result['ambiguous']
        elif not result['is_only_motif']:
            is_only_motif += 1
            #continue
            filters['is_only_motif'] = True
        if not result['events']:
            no_event += 1
            continue
        result['events'] = result['events'][0]
        if not last_result:
            last_result = result
        elif (variant.chrom == last_result['chrom']) and (abs(variant.pos - last_result['pos']) < motif_len):
            filters['interfered'] = True
            last_result = result
        logging.debug('result: {}'.format(result))
        result['wgs_tumor_ref_depth'] = variant.samples[0]['BC'][bc_lookup[variant.ref]]
        result['wgs_tumor_alt_depth'] = variant.samples[0]['BC'][bc_lookup[variant.alts[0]]]
        result['rna_tumor_ref_depth'] = variant.samples[1]['BC'][bc_lookup[variant.ref]]
        result['rna_tumor_alt_depth'] = variant.samples[1]['BC'][bc_lookup[variant.alts[0]]]
        result['chrom'] = variant.chrom
        # 1-based coord
        result['pos'] = variant.pos
        result['gene'] = gene
        result['strand'] = region.strand
        result['pid'] = args.pid
        result['cohort'] = args.cohort
        result['ref'] = variant.ref
        result['alt'] = variant.alts[0]
        result['effect'] = codes_event[result['events']['event_code']]
        result['pas_ref'] = result['events']['pas_ref']
        result['pas_alt'] = result['events']['pas_alt']
        result['pas_ref_strength'] = get_strength(result['pas_ref'], motifs)
        result['pas_alt_strength'] = get_strength(result['pas_alt'], motifs)
        # 1-based coords
        result['pas_start'] = result['events']['pos']
        result['pas_end'] = result['events']['pos'] + motif_len - 1
        result['pass'] = filters['pass']
        result['ambiguous'] = filters['ambiguous']
        result['interfered'] = filters['interfered']
        result['is_only_motif'] = filters['is_only_motif']
        result['previously_known'] = filters['previously_known']
        logging.debug('filters: {}'.format(filters))
        outfile.write(
            write_result(result, output_columns, args.delimiter)
        )

regions_file.close()
outfile.close()

# Stats file will record various stats about the run
#stats_path = os.path.join(args.outdir, '{}.stats'.format(args.name))
#stats_columns = (
#    'total_variants',
#    'in_regions',
#    'not_passed',
#    'previously_known',
#    'low_coverage',
#    'low_alt',
#    'multiple_alts',
#    'no_pas',
#    'ambiguous',
#    'interfere',
#    'is_only_motif',
#    'no_event',
#    'breakdown'
#)
## This variable will just be a comma-separated string showing how many
## variants are remaining after each filter
#diffs = [total_variants, in_regions]
#for i in xrange(2, len(stats_columns) - 1):
#    diff = diffs[i - 1] - eval(stats_columns[i])
#    diffs.append(diff)
#breakdown = (',').join([str(x) for x in diffs])
#with open(stats_path, 'w') as o:
#    o.write(
#        (args.delimiter).join(stats_columns) + '\n'
#    )
#    o.write(
#        (args.delimiter).join(
#            [str(eval(x)) for x in stats_columns]
#        ) + '\n'
#    )
