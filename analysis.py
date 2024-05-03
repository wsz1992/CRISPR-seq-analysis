# -*- encoding: utf-8 -*-


'''Command-line interface for CNV.'''


import logging
from . import commands


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    args = commands.parse_args()
    args.func(args)
