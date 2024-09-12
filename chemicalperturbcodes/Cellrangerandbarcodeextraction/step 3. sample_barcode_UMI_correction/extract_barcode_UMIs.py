#!/usr/bin/env python
# extract_barcode_UMIs.py
# extract barcode UMI counts based on 10X cellranger bam file
# 

import sys
import pandas as pd
import os.path
import logging
import gzip

import numpy as np
import scipy.sparse as ss
import scipy.io as sio

from re import search
from collections import Counter
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

# class
class MySamParser:
    """Class for extracting information from a sam file"""
    def __init__(self, args):
        """class constructor"""
        # input sam file
        self.sam_file = args.samfile
        # input cell barcode whitelist file
        self.whitelist_file = args.whitelistfile
        # input sample barcode file
        self.barcode_file = args.barcodefile
        # output umi counts matrix file
        self.umi_count_file_prefix = args.umicountfileprefix
        # output detailed standard sample barcode info file
        self.standard_table_file = args.stabfile
        # output detailed non-standard sample barcode info file
        self.nonstandard_table_file = args.nstabfile
        # start position of sample barcode in BFP sequence
        self.barcode_pos = 128
        # sample barcode length
        self.barcode_len = 10
        # mismatches allowed for sample barcode
        self.mismatch = 1
        # cell barcodes (dictionary: barcode --> '')
        self.cell_barcodes = None
        # sample barcodes (dictionary: id --> barcode)
        self.sample_barcodes = None
        # standard sample barcode counter: one dictionary per cell barode
        #                           |-- one dictionary per sample barcode
        #                             |-- one Counter per UMI
        self.bcount = {}
        # non-standard sample barcode counter: same structure as above
        self.nbcount = {}
        # standard sample barcode DataFrame (counter --> DataFrame)
        self.bcount_df = None
        # non-standard sample barcode DataFrame (counter --> DataFrame)
        self.nbcount_df = None

    def load_cell_barcode_whitelist(self):
        """load 10X cell barcodes"""
        # file exists?
        if not os.path.exists(self.whitelist_file):
            logging.error('Unable to load cell barcode whitelist file: {}.'.format(self.whitelist_file))
            sys.exit(1)
        # load cell barcodes
        with gzip.open(self.whitelist_file, 'rt') as fin:
            self.cell_barcodes = {x.strip():'' for x in fin.readlines()}
            logging.info('{} cell barcodes loaded.'.format(len(self.cell_barcodes)))
            #print(self.cell_barcodes[0])

    def load_sample_barcode_list(self):
        """load sample barcodes"""
        # file exists?
        if not os.path.exists(self.barcode_file):
            logging.error('Unable to load sample barcode file: {}.'.format(self.barcode_file))
            sys.exit(2)
        # load sample barcodes
        with open(self.barcode_file, 'r') as fin:
            self.sample_barcodes = {'myB{}'.format(x.split()[0]):x.split()[1] for x in fin.readlines()[1:]}
            logging.info('{} sample barcodes loaded.'.format(len(self.sample_barcodes)))
            #print(self.cell_barcodes)

    def fetch_tag_value(self, fields, ftag, ftype):
        """extract the value of a given tag from sam optional fields: TAG:TYPE:VALUE"""
        for item in fields:
            tpat = search('{}:{}:(.*)'.format(ftag,ftype), item)
            if tpat:
                #print(item)
                #print('{}:{}:(.*?)'.format(ftag,ftype))
                #print(tpat.groups())
                return tpat.groups()[0]
        return None

    def extract_read_info(self, read):
        """extract info from a read, return a dictionary"""
        elements = read.split()
        rinfo = {}
        # FLAG
        rinfo['FLAG'] = int(elements[1])
        # RNAME
        rinfo['RNAME'] = elements[2]
        # POS
        rinfo['POS'] = int(elements[3])
        # MAPQ
        rinfo['MAPQ'] = int(elements[4])
        # CIGAR
        rinfo['CIGAR'] = elements[5]
        # SEQ
        rinfo['SEQ'] = elements[9]
        # NH:i:1
        rinfo['NH'] = int(self.fetch_tag_value(elements[11:], 'NH', 'i'))
        # cell barcode
        rinfo['CB'] = self.fetch_tag_value(elements[11:], 'CB', 'Z')
        # molecular barcode
        rinfo['UB'] = self.fetch_tag_value(elements[11:], 'UB', 'Z')
        # gene id assigned, semicolon-separated
        rinfo['GX'] = self.fetch_tag_value(elements[11:], 'GX', 'Z')
        return(rinfo)

    def decode_cigar(self, rinfo):
        """decode CIGAR to get alignment info for fixing read"""
        k = 0
        prev = -1
        alninfo = []
        cigar = rinfo['CIGAR']
        while k != len(cigar):
            if cigar[k].isalpha():
                alninfo.append((int(cigar[prev+1:k]), cigar[k]))
                prev = k
            k += 1
        return alninfo

    def fix_seq(self, rinfo):
        """fix read sequence based on CIGAR so that it aligns with the reference sequence"""
        # decode CIGAR
        alninfo = self.decode_cigar(rinfo)
        # original sequence
        seq = rinfo['SEQ']
        # align read sequence
        pos = 0
        ret = []
        for alen, atype in alninfo:
            if atype == 'M':# match
                ret.append(seq[pos:pos+alen])
                pos += alen
            elif atype == 'I':# insertion
                pos += alen
            elif atype == 'D':# deletion
                ret.append('N'*alen)
            elif atype == 'N':# skip bases
                ret.append('N'*alen)
            elif atype == 'S':# soft clipping
                pos += alen
            elif atype == 'H':# hard clipping
                pos += alen
            elif atype == 'P':# padding
                pos += alen
            elif atype == '=':# sequence match
                ret.append(seq[pos:pos+alen])
                pos += alen
            elif atype == 'X':# sequence mismatch
                ret.append(seq[pos:pos+alen])
                pos += alen
            else:
                logging.error('Weird CIGAR character {}{}.'.format(alen,atype))
        return ''.join(ret)
                

    def prefilter_read(self, rinfo, read):
        """prefilter on read, decide whether we should discard it from UMI counting"""
        # -) Do NOT check whether or not mapped (assume all input reads are aligned)
        # 1) multiple mapping
        if rinfo['NH'] > 1:
            logging.warning('Multiple mapping read: {}'.format(read.rstrip()))
            return (False, 'NH')
        # 2) low quality?
        elif rinfo['MAPQ'] != 255:
            logging.warning('MAPQ != 255: {}'.format(read.rstrip()))
            return (False, 'MAPQ')
        # 3) no cell barcode
        elif rinfo['CB'] is None:
            logging.warning('Missing cell barcode: {}'.format(read.rstrip()))
            return (False, 'CB')
        # 4) not in the final cell barcode list
        elif rinfo['CB'] not in self.cell_barcodes:
            logging.warning('Not a valid cell barcode: {}'.format(read.rstrip()))
            return (False, 'CB')
        # 5) no molecular barcode
        elif rinfo['UB'] is None:
            logging.warning('Missing molecular barcode: {}'.format(read.rstrip()))
            return (False, 'UB')
        # 6) no gene id assigned
        elif rinfo['GX'] is None:
            logging.warning('No gene id assigned: {}'.format(read.rstrip()))
            return (False, 'GX')
        # 7) multiple gene ids assigned
        elif len(rinfo['GX'].split(';')) > 1:
            logging.warning('Multiple gene ids assigned: {}'.format(read.rstrip()))
            return (False, 'GX')
        ## 8) include indels in a read alignment?
        #elif 'I' in rinfo['CIGAR'] or 'D' in rinfo['CIGAR'] or 'N' in rinfo['CIGAR']:
        #    logging.warning('Indel or skips detected in read: {}'.format(read.rstrip()))
        #    return (True, 'PASS')
        # pass filter
        else:
            return (True, 'PASS')

    def extract_barcode(self, rinfo):
        """extract sample barcode from SEQ"""
        # barcode starting from the 128th position of the reference sequence, 10 bp
        pos = rinfo['POS']
        # fix the read sequence
        seq = self.fix_seq(rinfo)
        #print(rinfo['CIGAR'])
        #print(rinfo['SEQ'])
        #print(seq)
        # get barcode sequence (missing bases by 'N')
        barcode = []
        end = pos + len(seq) -1
        # read does NOT overlap with barcode bases
        if end < self.barcode_pos or pos > self.barcode_pos + self.barcode_len - 1:
            barcode.extend(['N' for k in range(self.barcode_len)])
        else:
            left = self.barcode_pos - pos
            right = end - (self.barcode_pos + self.barcode_len - 1)
            if left >= 0:
                if right >= 0:
                    barcode.append(seq[left:left+self.barcode_len])
                else:# right < 0
                    barcode.append(seq[left:])
                    barcode.extend(['N' for k in range(-right)])
            else:# left < 0
                if right >= 0:
                    barcode.extend(['N' for k in range(-left)])
                    barcode.append(seq[:self.barcode_len+left])
                else:# right < 0
                    barcode.extend(['N' for k in range(-left)])
                    barcode.append(seq)
                    barcode.extend(['N' for k in range(-right)])
        return ''.join(barcode)

    def hamming_distance(self, b1, b2):
        """calculate the Hamming distance between two barcodes"""
        return len(['' for k in range(len(b1)) if b1[k] != b2[k]])

    def process_sam(self):
        """process aligned reads in the sam file"""
        # file exists?
        if not os.path.exists(self.sam_file):
            logging.error('Unable to load sam file: {}'.format(self.sam_file))
            sys.exit(3)
        # process reads in sam file
        nFilt = Counter()
        with gzip.open(self.sam_file, 'rt') as fin:
            for read in fin:
                # extract info from the read
                rinfo = self.extract_read_info(read)
                #print(rinfo)
                # prefilter on read
                mypass, mylabel = self.prefilter_read(rinfo, read)
                nFilt[mylabel] += 1
                # not passing filter, discard!!!
                if not mypass:
                    continue
                # extract barcode
                barcode = self.extract_barcode(rinfo)
                #print(barcode)
                # sample barcode with allowed mismatch
                gid = rinfo['GX']
                ref_barcode = self.sample_barcodes[gid]
                bdiff = self.hamming_distance(ref_barcode, barcode)
                #with open('/tmp/barcodes.txt','a') as fout:
                #    fout.write('{}\t{}\t{}\n'.format(barcode,ref_barcode,bdiff))
                if bdiff <= self.mismatch:# match reference sample barcode
                    cb = rinfo['CB']
                    ub = rinfo['UB']
                    if cb not in self.bcount:
                        self.bcount[cb] = {}
                        self.bcount[cb][gid] = Counter()
                        self.bcount[cb][gid][ub] += 1
                    else:
                        if gid not in self.bcount[cb]:
                            self.bcount[cb][gid] = Counter()
                            self.bcount[cb][gid][ub] += 1
                        else:
                            self.bcount[cb][gid][ub] += 1
                else:# do NOT match reference sample barcode
                    if 'N' not in barcode:
                        cb = rinfo['CB']
                        ub = rinfo['UB']
                        if cb not in self.nbcount:
                            self.nbcount[cb] = {}
                            self.nbcount[cb][barcode] = Counter()
                            self.nbcount[cb][barcode][ub] += 1
                        else:
                            if barcode not in self.nbcount[cb]:
                                self.nbcount[cb][barcode] = Counter()
                                self.nbcount[cb][barcode][ub] += 1
                            else:
                                self.nbcount[cb][barcode][ub] += 1

    def counter_to_df(self, ct, colname='count', idxname='UB'):
        """convert Counter to DataFrame"""
        return pd.DataFrame({colname : pd.Series(ct)}).reset_index().rename(columns={'index':idxname})

    def dict_to_df(self, dt, colname='count', idxname='UB', keyname='gid'):
        """convert Dict to DataFrame"""
        df_list = []
        for k in dt:
            df = self.counter_to_df(dt[k], colname, idxname)
            df[keyname] = k
            df_list.append(df)
        # concat DataFrames
        return pd.concat(df_list, ignore_index=True)

    def count_to_df(self, count, colname='count', idxname='UB', keyname='gid', topname='CB'):
        """convert barcode count dictionary to DataFrame"""
        df_list = []
        for k in count:
            df = self.dict_to_df(count[k], colname, idxname, keyname)
            df[topname] = k
            df_list.append(df)
        # concat DataFrames
        return pd.concat(df_list, ignore_index=True)

    def bcount_to_df(self):
        """convert standard sample barcode count Dict to DataFrame"""
        self.bcount_df = self.count_to_df(self.bcount, colname='count', idxname='UB', keyname='gid', topname='CB')

    def nbcount_to_df(self):
        """convert non-standard sample barcode count Dict to DataFrame"""
        self.nbcount_df = self.count_to_df(self.nbcount, colname='count', idxname='UB', keyname='barcode', topname='CB')

    def bcount_df_to_tsv(self):
        """write standard sample barcode count DataFrame to file"""
        if self.standard_table_file is not None:
            self.bcount_df.to_csv(self.standard_table_file, sep='\t', columns=['CB','gid','UB','count'], index=False)

    def nbcount_df_to_tsv(self):
        """write non-standard sample barcode count DataFrame to file"""
        if self.nonstandard_table_file is not None:
            self.nbcount_df.to_csv(self.nonstandard_table_file, sep='\t', columns=['CB','barcode','UB','count'], index=False)

    def bcount_to_mtx(self):
        """write standard sample barcode count Dict to mtx file"""
        arlist = []
        cells = []
        barcodes = list(self.sample_barcodes.keys())
        for c in self.cell_barcodes:
            cells.append(c)
            if c not in self.bcount:# no data found for this cell barcode
                arlist.append([0 for k in range(len(self.sample_barcodes))])
            else:
                ar = []
                for s in barcodes:
                    if s not in self.bcount[c]:# no data found for this sample barcode
                        ar.append(0)
                    else:
                        ar.append(len(self.bcount[c][s]))
                arlist.append(ar)
        # convert array list to ndarray, transpose to barcode X cell
        npa = np.array(arlist).transpose()
        # convert to sparse matrix and write to file
        sio.mmwrite('{}.matrix.mtx'.format(self.umi_count_file_prefix), ss.csr_matrix(npa), field='integer')
        # write rownames (barcodes) and colnames (cells) to file
        with gzip.open('{}.barcodes.tsv.gz'.format(self.umi_count_file_prefix), 'wt') as fb, gzip.open('{}.features.tsv.gz'.format(self.umi_count_file_prefix), 'wt') as ff:
            ff.write('\n'.join(barcodes)+'\n')
            fb.write('\n'.join(cells)+'\n')

# functions
def get_arguments():
    """fetch command line arguments"""
    parser = ArgumentParser(description="""For Jeya's project: given a 10X cellranger bam, extract barcode UMI counts.""", prog='extract_barcodes_UMIs.py')
    parser.add_argument("-v", "--version", action="version", version='%(prog)s v0.1')
    parser.add_argument("-i", "--sam", nargs="?", required=True, help="sam file", metavar="sam_file", dest="samfile")
    parser.add_argument("-w", "--whitelist", nargs="?", required=True, help="cell barcode whitelist file", metavar="cell_barcode_file", dest="whitelistfile")
    parser.add_argument("-b", "--barcode", nargs="?", required=True, help="sample barcode file", metavar="sample_barcode_file", dest="barcodefile")
    parser.add_argument("-l", "--log", nargs="?", default="extract_barcode_UMIs.log", help="log file", metavar="log_file", dest="logfile")
    parser.add_argument("-o", "--umicount", nargs="?", required=True, help="output sample barcode umi counts matrix prefix", metavar="umi_count_fileprefix", dest="umicountfileprefix")
    parser.add_argument("-t", "--stable", nargs="?", default=None, help="output detailed standard sample barcode information", metavar="standard_barcode_table_file", dest="stabfile")
    parser.add_argument("-n", "--nstable", nargs="?", default=None, help="output detailed non-standard sample barcode information", metavar="non-standard_barcode_table_file", dest="nstabfile")
    return parser.parse_args()

def setup_logging(logfile, level=logging.INFO):
    """set up logging"""
    # prepare loggings
    log_formatter = logging.Formatter('%(levelname)s: %(message)s')
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # logging to file
    file_handler = logging.FileHandler(logfile)
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)

    # logging to stdout
    #console_handler = logging.StreamHandler(sys.stdout)
    #console_handler.setFormatter(log_formatter)
    #root_logger.addHandler(console_handler)

    # another common way to set up logging, but of less flexibility
    #logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    #logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    return root_logger

def main():
    """call me to get started!"""
    # load arguments
    args = get_arguments()
    # set up logging
    root_logger = setup_logging(args.logfile, level=logging.INFO)
    #root_logger = setup_logging(args.logfile, level=logging.DEBUG)
    #root_logger = setup_logging(args.logfile, level=logging.INFO)
    # create a MySamParser object
    d = MySamParser(args)
    # load cell barcodes
    d.load_cell_barcode_whitelist()
    # load sample barcodes
    d.load_sample_barcode_list()
    # process sam entries
    d.process_sam()
    ##print(d.bcount)
    ##print(d.dict_to_df(d.bcount['GTTGCTCAGCCATTCA-1'], colname='count', idxname='UB', keyname='gid'))
    ##print(d.count_to_df(d.bcount, colname='count', idxname='UB', keyname='gid', topname='CB'))
    # convert standard sample count Dict to DataFrame
    d.bcount_to_df()
    # convert non-standard sample count Dict to DataFrame
    d.nbcount_to_df()
    # write standard sample count DataFrame to file
    d.bcount_df_to_tsv()
    # write non-standard sample count DataFrame to file
    d.nbcount_df_to_tsv()
    # write standard sample count matrix to file
    d.bcount_to_mtx()

# main
if __name__ == '__main__':
    main()

