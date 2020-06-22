import os
import subprocess
import sys
import logging
import re

from genomic_tools_lib import Logging

def run(args):
    """
    1. Identify all DAP-G output in a directory tree
    2. For each phenotype, copy all of the SNP lines (from all of the
        different files) to a joined file.
    3. Execute William's perl script.
    :param args:
    :return:
    """
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    else:
        raise ValueError("Output directory already exists")



def execute_perl_script(perl_script, dir, vcf):
    cmd = ['perl', perl_script, '-dir', dir, '-vcf', vcf]
    logging.log(9, "Running reorganize Perl on {}".format(dir))
    return_state = subprocess.run(cmd, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    logging.log(8, "Finished reorganize")
    if return_state.returncode != 0:
        logging.warning(return_state.stdout)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-file_pattern", nargs='+')
    parser.add_argument("-perl_script_path")
    parser.add_argument("-output_dir")
    parser.add_argument("--parsimony", type=int, default=10)

    args = parser.parse_args()

    Logging.configure_logging(level = args.parsimony, target=sys.stdout)

