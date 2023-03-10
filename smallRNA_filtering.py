import sys # we need sys to access file input parameters
import time # let's track start time and duration
import json # create a report at the end in JSON format
import gzip
from collections import defaultdict

# returns phred score for a given character using Illumina's Quality Score Encoding
def convertQualityCharacterToScore(character):
    return ord(character)-33


# returns the lowest score for a given quality score string
def getMinimumQualityScore(quality_string):
    return min([convertQualityCharacterToScore(quality_character) for quality_character in quality_string])


# expects a fastq.gz file and creates two different counting dictionaries for high and low quality reads
# based on the given quality threshold (default is 20). Uses defaultdict for sequence counting
def getLowAndHighQualityReads(input_file, quality_threshold):
    sequence2count_high = defaultdict(int)
    sequence2count_low = defaultdict(int)
    error_line_count = 0
    
    # we open the compressed file and read it as a text.
    with gzip.open(input_file, 'rt') as f:
        for line in f: # instead of reading whole file into memory, we will proceed line by line
            if line.startswith('@'): # we found the header line
                sequence = f.readline().strip() # read sequence line
                quality_header = f.readline().strip() # read 3rd line, cutadapt converts this line to + symbol
                quality_string = f.readline().strip() # read phred score string line
                
                minimum_read_score = getMinimumQualityScore(quality_string)
                
                # initialize or increment sequence counts 
                if  minimum_read_score >= quality_threshold: # high-quality read sequence counting
                    sequence2count_high[sequence] += 1 # we increment the count of the sequence in dictionary
                else:
                    sequence2count_low[sequence] += 1 # we increment the count of the sequence in dictionary
            else: # this part is for catching errors in FASTQ file
                error_line_count += 1
                    
                
    return sequence2count_high, sequence2count_low, error_line_count


# main function to manage filtering and create fasta file
def getHighQualityReadsAndCreateReport(fastq_gz_path, quality_threshold, fasta_path):
    report_dict = {} # we will keep stats here instead of printing
    report_dict[f'fastq.gz file used for the analysis is'] = fastq_gz_path
    
    sequence2count_high, sequence2count_low, error_line_count = getLowAndHighQualityReads(fastq_gz_path, quality_threshold)
    
    # report high- and low-quality reads based on the threshold
    # instead of printing, add to the report_dict
    sum_high_quality_reads = sum(sequence2count_high.values())
    sum_low_quality_reads = sum(sequence2count_low.values())
    total_reads = sum_high_quality_reads + sum_low_quality_reads
    perc_high_quality_reads = sum_high_quality_reads/(sum_high_quality_reads+sum_low_quality_reads)

    report_dict['total reads'] = total_reads
    report_dict['high-quality reads'] = sum_high_quality_reads
    report_dict['low-quality reads'] = sum_low_quality_reads
    report_dict['% high-quality reads'] = 100*perc_high_quality_reads


    # get common sequences in both high- and low- quality reads and merge read counts into a new dictionary
    high_quality_sequences = sequence2count_high.keys()
    low_quality_sequences = sequence2count_low.keys()
    both_high_and_low_quality_sequences = set(high_quality_sequences) & set(low_quality_sequences)
    sum_common_low_quality_reads = sum([sequence2count_low[sequence] for sequence in both_high_and_low_quality_sequences])

    report_dict['common low-quality reads'] = sum_common_low_quality_reads

    
    sequence2count = sequence2count_high.copy() # use copy instead of assignment, you can read about copying
    for sequence in both_high_and_low_quality_sequences: # operate on common sequences in high and low quality
        sequence2count[sequence] += sequence2count_low[sequence] # we add counts from low-quality dictionary
    perc_merged_reads = sum(sequence2count.values())/total_reads

    report_dict['reads after merge'] = sum(sequence2count.values())
    report_dict['% reads after merge'] = 100*perc_merged_reads
    print(report_dict)

    # convert sequence:count dictionary to a fasta file
    fasta_reads = [f'>{sequence}:{sequence2count[sequence]}\n{sequence}' for sequence in sequence2count]
    print(fasta_reads[:5])

    # save fasta file
    with open(fasta_path, 'w') as f:
        f.write('\n'.join(fasta_reads))
    
    return report_dict


# this is how we run Python files, you can read why it's structured like that
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print('You need to specify <fastq.gz file path> and optionally <quality score threshold> and <fasta file path>')
        print('This script takes a smallRNA NGS fastq.gz format file and filters anything below given threshold')
        print('Default quality score is 20. Default fasta file name is created by replacing the "fastq.gz" with "fa"')
        print('Reads passing filter are saved as a sequence:count-header fasta file.')
        sys.exit()
    else:
        fastq_gz_path = sys.argv[1] # get fastq.gz file path for processing, no error checking here
        fasta_path = None
        quality_threshold = None
        
        try: # if quality score threshold is supplied, use it. else, set it to 20
            quality_threshold = int(sys.argv[2])
            try: # if fasta file path is supplied, use it. else, replace "fastq.gz" with "fa" in file_path 
                fasta_path = int(sys.argv[3])
            except:
                # assume fastq.gz is only at the end of the path; we will cover it later
                fasta_path = fastq_gz_path.replace('.fastq.gz', '.fa')
        except:
            quality_threshold = 20
            fasta_path = fastq_gz_path.replace('.fastq.gz', '.fa')
    json_report_path = fasta_path.replace('.fa', '.filtering.report.json')
    parameters_dict = {'quality_threshold':quality_threshold, 
                       'fastq.gz_path':fastq_gz_path, 
                       'fasta_path':fasta_path,
                       'json_report_path':json_report_path,
                       'start_time':time.time()}

    print(parameters_dict)
        
    # get low and high quality reads
    report_dict = getHighQualityReadsAndCreateReport(fastq_gz_path, quality_threshold, fasta_path)
    parameters_dict['duration'] = int(time.time()-parameters_dict['start_time'])
    
    # save parameters and report in JSON format to the folder fasta file is saved to 
    output_dict = {'parameters':parameters_dict, 'report':report_dict}

    with open(json_report_path, 'w') as f:
        json.dump(output_dict, f)
