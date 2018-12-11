__author__ = "alvaro barbeira"

import logging
import sys

def configure_logging(level=5, target=sys.stderr, log_file=None, with_date=False):
    logger = logging.getLogger()
    logger.setLevel(level)

    # create console handler and set level to info
    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    if with_date:
        formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
                                      datefmt='%Y-%m-%d %H:%M:%S')
    else:
        formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if log_file:
        fileHandler = logging.FileHandler("log_file")
        fileHandler.setFormatter(formatter)
        logger.addHandler(fileHandler)