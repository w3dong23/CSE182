# this file will download fg38 ref genome sequence by chromosome number
# unzip files based on user input file on marker of sliver sequence
# reference: https://gist.github.com/kstreepy/a9800804c21367d5a8bde692318a18f5

import os, gzip, shutil
import requests
import pandas as pd

# this function will create dictionary of mapping chr number to chr filename
def constructDict():
    chr_name = []
    # generate chromosome names
    for i in range(1, 23, 1):
        chr_name.append(str(i))
    chr_name.append('MT')
    chr_name.append('X')
    chr_name.append('Y')

    # intialize chromosome dict
    chr_seq_dict = {}
    for name in chr_name:
        chr_seq_dict[name] = 'hg38.chr{}.fa.gz'.format(name)
    return chr_name, chr_seq_dict

# construct name and dictionary for mapping chromosome number to filename 
chr_name, chr_seq_dict = constructDict()


# this function will download reference human genome file to current directory
def downloadGenome():
    for name in chr_name:
        url = 'https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{}.fa.gz'.format(name)
        r = requests.get(url, allow_redirects=True)
        # same the ref genome in fasta format
        filename =  chr_seq_dict[name]
        open(filename, 'wb').write(r.content)

# this function will return DNA sequence based on the marker file
# marker file format example: Chr7 1 9
def findSliverSeq(marker_file):
    markers = open(marker_file).readlines()
    result_seqs = []

    # record information about sliver seq to write to csv file later 
    chr_header = []
    start_pos = []
    end_pos = []

    for marker in markers:
        marker_list = marker.strip().split(' ')

        # record info on input file to output to df later
        chr_header.append(marker_list[0])
        start_pos.append(int(marker_list[1]))
        end_pos.append(int(marker_list[2]))

        # find actual sliver sequences 
        chr_num = marker_list[0].replace('Chr', '')
        start = int(marker_list[1])
        length = int(marker_list[2]) - start

        # fetch gzip file name from dict
        gzip_file = chr_seq_dict[chr_num]
        # unzipped file name no longer have gz extension (remove.gz at end)
        filename = gzip_file[:-3]
        # unzip corresponding file if it have not been unziped
        if (os.path.exists(gzip_file)):
            # unzip the file 
            with gzip.open(gzip_file,"rb") as f_in, open(filename,"wb") as f_out:
              shutil.copyfileobj(f_in, f_out)
            # remove zipped files 
            os.remove(gzip_file)

        # if it is already unzipped read in contents and find sequence
        # ignore header, read in sequences of chromosome
        chr_seq = open(filename).readlines()[1:]
        # our seq have 60 char per line
        # so we can calculate rough start line num 
        start_line_num = int(start / 60)
        # start_line_index refers on start index on start_line 
        start_line_index = start - start_line_num * 60 - 1

        # read in first line of sliver_seq based on start index
        if (start_line_index >= 0):
            if (length >= 60):
                sliver_seq = chr_seq[start_line_num].strip()[start_line_index:]
            # if length < 60, only read part of start line 
            else:
                sliver_seq = chr_seq[start_line_num].strip()[start_line_index:start_line_index+length]
        else:
            sliver_seq = chr_seq[start_line_num].strip()[start_line_index]
    
        # remaining = rest # of base to read
        remaining = length - len(sliver_seq)
        # update curr line number to read 
        curr_line_index = start_line_num + 1
        # while reamaining >= 60 char, read whole lines
        while(remaining >= 60):
            curr_seq = chr_seq[curr_line_index].strip()
            sliver_seq = sliver_seq + curr_seq
            curr_line_index+=1
            remaining = remaining - 60
        
        # if there is still sth to read after reading whole lines
        if (remaining > 0):
            # read the last sequences to sliver_seq
            sliver_seq = sliver_seq + chr_seq[curr_line_index].strip()[:remaining]
        result_seqs.append(sliver_seq)
    print(result_seqs)

    # append all result sequences to a text file line by line 
    # format like: 
    # ACTCTCTC
    # AACTTTTCCCC
    output = ''
    for sliver_seq in result_seqs:
        output = output + sliver_seq + '\n'
    
    # open result file and write to result file with all sliver seuqneces 
    with open('sliver.txt', 'w') as result_file:
        result_file.write(output)
    result_file.close()

    # open a csv file for storing all info about sliver sequence found 
    # cotains chromosome pos and actual sliver sequnce 
    df = pd.DataFrame()
    df['Chromosome Position'] = chr_header
    df['Start Position'] = start_pos
    df['End Position'] = end_pos
    df['Sequence'] = result_seqs
    df.to_csv('sliver.csv', index=False)

    return result_seqs


#uncomment to download genome files *YOU ONLY NEED TO DO IT ONCE*
#downloadGenome() # ONCE YOU FINISH, COMMENT IT OUT

# used to test functions 
#findSliverSeq('test/example.txt')
#findSliverSeq('test/ex2.txt')










