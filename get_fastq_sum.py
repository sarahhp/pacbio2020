"""
Merge key statistics from pacbio isoseq fastq transcript headers
"""
__author__ = "Sarah Hazell Pickering"
__date__ = "2018-09-27"

import os
import numpy
import pandas


def _getNum(attributes):
    nums = []
    for a in attributes:
        nums.append(int(a.split(sep="=")[1]))
    return nums

def _getAverage(df):
    averages = df.mean(axis=0, numeric_only=True)
    averages = averages.round()
    formatted = pandas.Series(["averages"] + averages.to_list(), index =df.columns)
    return formatted

def getIsoseqSummary(fasta, outpath="test.txt"):
    columns = ["transcript","full_length_coverage","length","num_subreads"]
    collate = []
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith(">"):
  #              breakpoint()
                trans_id,other = line.split(" ")
                trans_id = trans_id[1:] #rm ">"
                fl_coverage,length,num_subreads = _getNum(other.split(sep=";"))
                collate.append([trans_id, fl_coverage,length,num_subreads])
    result = pandas.DataFrame(collate, columns=columns)
    averages = _getAverage(result) 
#    breakpoint()
    result = result.append(averages, ignore_index=True)
    result = result.round()
    result.to_csv(outpath, sep="\t", index=False)
    return averages


def extractAverages(summary="hisat_summary.tab"):
    averages = pandas.read_table(summary).iloc[-1,1:] #select last row of summary file
    return averages
    
def countFullLength(isoseq, cupcake, outpath):
    #convert cupcake summary to cluster:[transcripts] dict
    with open(cupcake) as fin:
        rows = ( line.split('\t') for line in fin )
        c = { row[0]:(row[1].strip().split(',')) for row in rows }

    #convert isoseq summary to transcript:full_length dict
    with open(isoseq) as fin:
        rows = ( line.split('\t') for line in fin )
        i = { row[0]:row[1].strip() for row in rows } 
    collate = {}
    for cluster in c: #c.keys()?
        count = 0
        for trans in c[cluster]:
            count += int(i[trans].split('.')[0])
        collate[cluster] = count
    result = pandas.Series(collate)
    print(result)
    return result 


def addCoverageToGFF(summary, gtf, output_file):
    #convert cupcake summary to cluster:full_length_coverage dict
    with open(summary) as coverage:
        rows =  (line.split('\t') for line in coverage if not(line.startswith('#')))
        cov = { row[0]:row[1].strip() for row in rows }
    
    collate = []
    with open(gtf) as gff:
        for line in gff:
            cols = line.rstrip().split('\t')
            attrs = cols[8].split(";")
            for at in attrs:
                at = at.strip()
                if at.startswith("transcript_id"):
                    trans_id = at.split()[-1].strip('"')
            cols[8] = cols[8] + " full_length_coverage \"" + cov[trans_id] + '"'
            collate.append('\t'.join(cols) + '\n')
    with open(output_file, 'w') as out:
        out.writelines(collate)
    return collate[1:6]

if __name__ == '__main__':
    data = "new_transcripts/ASC_D07G_P7_D1_hq.fasta"
    getIsoseqSummary(data)
    #extractAverages(summary="test.txt")
    #countFullLength("test.txt", "cupcake/P7_D1/P7_D1_250bp5.collapsed.group.txt", "")
    PRE = "cupcake/ASC_D07G_P7_D1/P7_D1.collapsed"
    print(addCoverageToGFF(PRE + ".abundance.txt",PRE + ".gff", PRE + ".f.gff"))
