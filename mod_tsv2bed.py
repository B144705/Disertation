#! /usr/bin/env python3
import sys
import csv
import argparse
def parseArgs():
    parser = argparse.ArgumentParser( description='Generate methylation bedGraph file')
    parser.add_argument('-v', '--verbose', action="store_true",required=False, default = False,
            help="verbose output")
    parser.add_argument('-i', '--input', type=str, required=True,help="input methylation tsv file (default stdin)")
    #parser.add_argument('-q', '--mod',type=str,required=False,default='cpg',help="modification motif; one of cpg,gpc,dam,cpggpc")
    parser.add_argument('-g','--genome',type=str,required=False,help="Reference genome fasta")
    #parser.add_argument('-e', '--exclude',type=str,required=False,help="motif to exclude from reporting")
    #parser.add_argument('-w', '--window',type=int,required=False,default=2,
    #        help="number of nucleotides to report on either side of called nucleotide")
    parser.add_argument('--nome',action="store_true",required=False,default=False,
            help="nanoNOMe - remove calls at GCG motifs")
    parser.add_argument('--offset',type=int,required=False,default=0,
            help="nanopolish coordinate offset (1-based)")
    args = parser.parse_args()
    return args

class readQuery:
    def __init__(self,record,offset,nome):
        self.offset = offset
        self.nome = nome
        self.qname = record['read_id']
        self.rname = record['chrom']
        self.start = int(record['ref_position']) + self.offset
        self.end = self.start+1
        self.current_pos = int(record['forward_read_position'])
        self.dist=[]
        self.seq=[]
        self.prob=[]
        self.call=[]
    def update(self,record):
        prob=float(record['call_prob'])
        if record['call_code'] == '-':
            call='u'
        else:
            call=record['call_code']
        sequence = record['ref_kmer']
        record_pos = int(record['forward_read_position'])
        dist = record_pos - self.current_pos
        # update
        self.current_pos = record_pos
        self.end = int(record['ref_position']) + 1
        self.seq.append(sequence)
        self.prob.append(prob)
        self.call.append(call)
        self.dist.append(abs(dist))
    def summarizeCall(self):
        summary=""
        for ind,call in zip(self.dist,self.call):
            summary=summary+"{}{}".format(ind,call)
        return summary

    def printRead(self,extra=""):
        '''''
        if len(self.call) == 0 : return # after filtering GCG there is no data in this read
        # if GCG is first motif, adjust start to make first index 0
        if self.dist[0] != 0 :
            self.start += self.dist[0]
            self.dist[0] = 0
        # when motif is GC, positions need to shift + one
        if self.motif == "GC" :
            self.start += 1
            self.end += 1 
        '''''
        def catList(strlist,delim=""):
            return delim.join([str(x) for x in strlist])
        print("\t".join([self.rname,str(self.start),str(self.end),
            self.qname,
            self.summarizeCall(),
            catList(self.prob,","),
            catList(self.seq,",")])+extra)
def summarizeMeth(args):
    if args.input:
        in_fh = open(args.input)
    csv_reader = csv.DictReader(in_fh, delimiter='\t')
    if args.verbose: print("processing calls", file=sys.stderr)
    for record in csv_reader:
        try:
            if ((record['read_id'] != read.qname) or
                    (record['chrom'] != read.rname) or
                    (int(record['ref_position']) < read.end)):
                read.printRead()
                read = readQuery(record,args.offset,args.nome)
        except NameError:  # has not been initialized
            read = readQuery(record, args.offset, args.nome)
        except ValueError:  # header or otherwise somehow faulty
            continue
        read.update(record)
    if read.qname:
        read.printRead()
    in_fh.close()

if __name__=="__main__":
    args=parseArgs()
    summarizeMeth(args)