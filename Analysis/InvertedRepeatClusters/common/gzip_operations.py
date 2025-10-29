#!/usr/bin/env python3.8
"""
This module reads and writes gzipped files
"""
__author__ = 'Roni'


# Function for wrapping subprocesses opening
import gzip
import os
import shutil


def gzip_and_remove_original(file_path):
    # zip file
    with open(file_path, 'rb') as f_in, gzip.open(file_path + ".gz", 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    # remove original
    os.remove(file_path)


def gunzip_file(file_path, suffix='gz'):
    # unzip zip file: remove suffix and '.' from end
    with gzip.open(file_path, 'rb') as f_in, open(file_path[:-(len(suffix) + 1)], 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return file_path[:-3]


def write_gzip_file(file_path, output):
    with gzip.open(file_path, "wt") as f:
        f.write(output)
