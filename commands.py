# -*- encoding: utf-8 -*-


import argparse
import glob
import os

import pandas as pd

from . import logo_script,cmdutil


parser = argparse.ArgumentParser(description='Analysis pipeline.')
subparser = parser.add_subparsers(
    help='Sub-commands (use with -h for more info)')


def cmd_preprocess(args):
    '''Run PAM pre-analysis.'''
    cmdutil.do_mkdir(args.out_path)
    cmdutil.do_mkdir('{}/Merged'.format(args.out_path))
    cmdutil.do_mkdir('{}/Unique'.format(args.out_path))

    df_meta = pd.read_csv(args.sample, sep='\t')
    for sample_name, _pam in df_meta.itertuples(index=False):
        tmp_list = glob.glob('{}/{}/*.gz'.format(args.in_path, sample_name))
        read_list = list(sorted(tmp_list))
        fastp_out = '{}/Merged/{}.fq'.format(args.out_path, sample_name)
        logo_script.do_fastp(read_list, fastp_out, args.processes)
        seqtk_out = '{}/Merged/{}_q10.fq'.format(args.out_path, sample_name)
        logo_script.do_seqtk(fastp_out, seqtk_out)
        df_unique = logo.fastq2df(seqtk_out)
        unique_out = '{}/Unique/unique_{}_q10.pkl'.format(args.out_path, 
                     sample_name)
        df_unique.to_pickle(unique_out)

p_preprocess = subparser.add_parser('preprocess', help=cmd_preprocess.__doc__)
p_preprocess.add_argument('-i', '--in-path', help='Input samples fastq path.')
p_preprocess.add_argument('-s', '--sample', help='Sample information file.')
p_preprocess.add_argument('-o', '--out-path', help='Result output path.')
p_preprocess.add_argument('-p', '--processes', default=1, type=int, 
                       help='Multiple processes.')
p_preprocess.set_defaults(func=cmd_preprocess)


def cmd_weblogo(args):
    '''Plot weblogo.'''
    cmdutil.do_mkdir('{}/triple_valid'.format(args.out_path))
    cmdutil.do_mkdir('{}/logos'.format(args.out_path))
    df_meta = pd.read_csv(args.sample, sep='\t')
    df_meta['file'] = df_meta['sample'].apply(
                lambda x:'{}/Unique/unique_{}_q10.pkl'.format(args.out_path, x))
    
    df_meta.apply(lambda x: logo.do_weblogo(x, args.out_path, args.processes), 
                  axis=1)

p_weblogo = subparser.add_parser('weblogo', help=cmd_weblogo.__doc__)
p_weblogo.add_argument('-s', '--sample', help='Sample information file.')
p_weblogo.add_argument('-o', '--out-path', help='Result output path.')
p_weblogo.add_argument('-p', '--processes', default=1, type=int, 
                       help='Multiple processes.')
p_weblogo.set_defaults(func=cmd_weblogo)


def parse_args(args=None):
    '''Parse the command line.'''
    return parser.parse_args(args=args)
