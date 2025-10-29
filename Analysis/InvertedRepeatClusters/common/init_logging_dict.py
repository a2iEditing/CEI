#!/usr/bin/env python3.8
"""
This function is for initializing a logging object.
"""

import logging
import os

def init_logging_dict(log_file):
    if not os.path.isdir(os.path.dirname(log_file)):
        os.makedirs(os.path.dirname(log_file))
    logging.basicConfig(level=logging.DEBUG,
                        format='[%(asctime)s] %(processName)s %(module)s %(levelname)-8s %(message)s',
                        datefmt='%y-%m-%d %H:%M:%S',
                        filename=log_file,
                        filemode='w')
    # define a Handler which writes WARN messages or higher to the sys.stderr
    console_logger = logging.StreamHandler()
    console_logger.setLevel(logging.WARN)
    # set a format which is simpler for console_logger use
    formatter = logging.Formatter('[%(asctime)s] %(module)s %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console_logger.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console_logger)
