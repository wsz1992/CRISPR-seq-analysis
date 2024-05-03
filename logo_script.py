import functools
import io
import glob
import logomaker
import operator
import os
import subprocess
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from Bio import motifs
from Bio.Seq import Seq
from cutadapt.adapters import (FrontAdapter, BackAdapter, LinkedAdapter)
from dnaio import Sequence
from functools import partial
from multiprocessing import Pool
from . import cmdutil
import argparse
import logging
from lib import commands


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    args = commands.parse_args()
    args.func(args)

parser = argparse.ArgumentParser(description='Analysis pipeline.')
subparser = parser.add_subparsers(
    help='Sub-commands (use with -h for more info)')
PAM_SEQ = {
    'PAM1': ['GAACGGCTCGGAGATCATCATTGCG', 'TTGTTTGCCACCATG', 'GTGAGCAAGGGCGAG'],
    'PAM2': ['GAGAGTAGAGGCGGCCACGACCTG', 'TTGTTTGCCACCATG','GTGAGCAAGGGCGAG'], 
    'PAM3': ['GAGAGTAGAGGCGGCCACGACCTG', 'TTGTTTGCCACCATG','GTGAGCAAGGGCGAG'], 
    'PAM4': ['GGGCTTCAAGCAACTTGTAGTGGG','TTGTTTGCCACCATG','GTGAGCAAGGGCGAG'],
    'PAM5': ['GAACGGCTCGGAGATCATCATT', 'TTGTTTGCCACCATG', 'GTGAGCAAGGGCGAG'],
    'PAM6': ['GGAGAGTAGAGGCGGCCACGAC', 'TTGTTTGCCACCATG', 'GTGAGCAAGGGCGAG']
}
def get_read_info(row, pam='PAM1'):
    desinged_target = PAM_SEQ[pam][0]
    seq = row['read']
    prefix_range = ()
    suffix_ragne = ()
    target_range = ()
    read_prefix_seq = ''
    read_suffix_seq = ''
    read_target_seq = ''
    target_len = 0
    is_edit = -1
    is_3fold = 0
    #=====================
    #5' 3' adapter 公用参数
    parameters= {'max_error_rate': 0.2,
                 #与params.PAM_SEQ[0]相关
                 'min_overlap': len(desinged_target),
                 'read_wildcards': False,
                 'adapter_wildcards': False,
                 'indels': False}
    
    #这里不考虑 anchor 问题，如果需要 anchor 则用 PrefixAdapter
    front_adapter = FrontAdapter(PAM_SEQ[pam][1], **parameters)
    #这里不考虑 anchor 问题，如果需要 anchor 则用 SuffixAdapter
    back_adapter = BackAdapter(PAM_SEQ[pam][2], **parameters)

    linked_adapter = LinkedAdapter(front_adapter, back_adapter, 
                                   front_required=True, back_required=True, 
                                   name='target_region_recognition')
    #匹配target区域前后的位置
    r = linked_adapter.match_to(seq)
    #如果有匹配到
    if r is not None:
        prefix_range = (r.front_match.rstart, r.front_match.rstop, 
                        r.front_match.errors)
        suffix_ragne = (r.back_match.rstart + r.front_match.rstop, 
                        r.back_match.rstop + r.front_match.rstop, 
                        r.back_match.errors)
        target_range = (r.front_match.rstop, 
                        r.back_match.rstart + r.front_match.rstop)
        read_prefix_seq = seq[prefix_range[0]:prefix_range[1]]
        read_suffix_seq = seq[suffix_ragne[0]:suffix_ragne[1]]
        read_target_seq = seq[target_range[0]:target_range[1]]
        target_len = len(read_target_seq)
        #判断与设计的gRNA 序列是否一致
        if not desinged_target in read_target_seq:
            is_edit = 1
        #判断是否 3 倍变化
        if (target_len % 3 == 0) and ('N' not in read_target_seq):
            is_3fold = 1
    return (read_prefix_seq, read_suffix_seq, read_target_seq, target_len, 
            is_edit, is_3fold)

def get_GCG_3p(row, pam='PAM1'):
    if pam in ['PAM1', 'PAM2']:
        if pam == 'PAM1':
            end = -10
            length = 7
            vlid_triple = 'GCG'
        else:
            end = -11
            length = 8
            vlid_triple = 'CTG'
        r = row[end:]
        is_valid_triple = 0
        if vlid_triple == r[:3] and len(r[3:])==length:
            is_valid_triple = 1
        return r[3:],is_valid_triple
    elif pam in ['PAM3', 'PAM4']:
        if pam == 'PAM3':
            end = 8
            length = 5
            vlid_triple = 'GGA'
        else:
            end = 8
            length = 5
            vlid_triple = 'GGG'
        r = row[:end]
        is_valid_triple = 0
        if vlid_triple == r[5:]:
            is_valid_triple = 1
        return r[:5], is_valid_triple
    else:
        if pam == 'PAM5':
            end = -13
            length = 10
            vlid_triple = 'ATT'
        else:
            end = -14
            length = 11
            vlid_triple = 'GAC'
        r = row[end:]
        is_valid_triple = 0
        if vlid_triple == r[:3] and len(r[3:])==length:
            is_valid_triple = 1
        return r[3:],is_valid_triple
    
def get_Seq_logo(row):
    seq = [Seq(row['7mer'])]
    seqs = seq * row['counts']
    return seqs

def get_info(data, pam='PAM1'):
    # para_i = 0 if pam == 'PAM1' else 1
    data[['read_prefix_seq', 'read_suffix_seq', 'read_target_seq', 
          'target_len', 'is_edit', 'is_3fold']] = data.apply(
                        lambda x: pd.Series(get_read_info(x, pam)), axis=1)
    return data

def parallelize_dataframe(df, func, pam, processes=1):
    num_partitions = processes
    num_cores = processes
    df_split = np.array_split(df, num_partitions)
    pool = Pool(num_cores)
    #固定函数，由参数传递给偏函数两个值,固定 PAM 以后 patial_func 只会有 data 一个参数了
    partial_func = partial(get_info, pam=pam)
    df = pd.concat(pool.map(partial_func, df_split))
    pool.close()
    pool.join()
    return df

def do_weblogo(row, out_path, processes):
    file = row['file']
    pam = row['pam']
    out_file = '{}/triple_valid/{}.pkl'.format(out_path, row['sample'])
    out_file1 = '{}/triple_valid/{}.csv'.format(out_path, row['sample'])

    if os.path.exists(out_file):
        print('File {} exits'.format(out_file))
        df_gcg = pd.read_pickle(out_file)
    else:
        df_unique = pd.read_pickle(file)
        print(pam)
        # if pam == 'PAM2':
        #     df_unique = df_unique[df_unique['read'].str.contains('GAGAGT')]
        # else:
        #     df_unique = df_unique[df_unique['read'].str.contains('GAACGG')]
        df_info = parallelize_dataframe(df_unique, get_info, pam, processes)
        print('Processing 3fold valid reads: {}'.format(out_file))
        df_3fold = df_info[df_info['is_3fold'] == 1].copy()    
        df_3fold[['7mer','is_GCG']] = df_3fold.read_target_seq.apply(
                                        lambda x: pd.Series(get_GCG_3p(x, pam)))
        df_gcg = df_3fold[df_3fold['is_GCG'] == 1].iloc[:, 1:]
        df_gcg.to_pickle(out_file)
    #weblogo
    print('Processing seqlogo')
    df_gcg.reset_index(drop=True).to_csv(out_file1, sep='\t', index=False)
    s_seqs = df_gcg.apply(lambda x: get_Seq_logo(x), axis=1)
    instances = functools.reduce(operator.iconcat, list(s_seqs), [])
    kargs = {'color_scheme':'color_custom',
             'symbols0': 'T',
             'symbols1': 'A',
             'symbols2': 'C',
             'symbols3': 'G',
             'color0': '#FD8D3C',
             'color1': '#74C476',
             'color2': '#9170C0',
             'color3': '#6BAED6',
            }
    print(len(s_seqs))
    if len(s_seqs)>0:
        with open('{}/triple_valid/{}.fa'.format(out_path,row['sample']),'w')as f:
            for i in instances:
                f.write(str(i)+'\n')
        m = motifs.create(instances)
        cmd='weblogo --format pdf  --color Orange T \'T\' --color SpringGreen A \'A\' --color DarkOrchid C \'C\' --color CornflowerBlue G \'G\' < {}/triple_valid/{}.fa > {}/logos/{}.pdf'.format(out_path,row['sample'], out_path,row['sample'])
        print(cmd)
        os.system(cmd)    
def cmd_weblogo(args):
    '''Plot weblogo.'''
    cmdutil.do_mkdir('{}/triple_valid'.format(args.out_path))
    cmdutil.do_mkdir('{}/logos'.format(args.out_path))
    df_meta = pd.read_csv(args.sample, sep='\t')
    df_meta['file'] = df_meta['sample'].apply(
                lambda x:'{}/Unique/unique_{}_q10.pkl'.format(args.out_path, x))
    
    df_meta.apply(lambda x: do_weblogo(x, args.out_path, args.processes), 
                  axis=1)

p_weblogo = subparser.add_parser('weblogo', help=cmd_weblogo.__doc__)
p_weblogo.add_argument('-s', '--sample', help='Sample information file.')
p_weblogo.add_argument('-o', '--out-path', help='Result output path.')
p_weblogo.add_argument('-p', '--processes', default=1, type=int, 
                       help='Multiple processes.')
p_weblogo.set_defaults(func=cmd_weblogo)
