#!/usr/bin/env python3.8
"""
Module for working with BED files using server in Python
"""

# region Imports
import sys
import os
sys.path.append(os.path.dirname(__file__))
from command_executer import execute
# endregion

# region Consts
BEDTOOLS_PROGRAM = "bedtools"
INTERSECT_PROGRAM = "%s intersect"
SORT_PROGRAM = "%s sort"
MERGE_PROGRAM = "%s merge"
STDIN_STR = "stdin"
# endregion


def intersect_files(repeat_bed_path, regions_bed_path, extra_params_list=[], stdin=None, bedtools=BEDTOOLS_PROGRAM):
    # run command
    if stdin:
        # check for input validity
        assert repeat_bed_path == STDIN_STR or regions_bed_path == STDIN_STR
        return execute(" ".join([INTERSECT_PROGRAM % bedtools, "-a", repeat_bed_path, "-b", regions_bed_path] +
                                extra_params_list), input_str=stdin).strip().split("\n")
    return execute(" ".join([INTERSECT_PROGRAM % bedtools, "-a", repeat_bed_path, "-b", regions_bed_path] +
                     extra_params_list)).strip().split("\n")

def sort_and_merge_bed_str(bed_str, bedtools=BEDTOOLS_PROGRAM):
    return merge_bed_str(sort_bed_str(bed_str, bedtools=bedtools), bedtools=bedtools)

'''
Using mutable default argument - single mutable for all empty calls of method.
In this case it doesn't matter as we do not change the list
'''
def sort_bed_str(bed_str, extra_params_list=[], bedtools=BEDTOOLS_PROGRAM):
    return execute(" ".join([SORT_PROGRAM % bedtools, "-i", "stdin"] + extra_params_list), input_str=bed_str)

def merge_bed_str(bed_str, extra_params_list=[], bedtools=BEDTOOLS_PROGRAM):
    return execute(" ".join([MERGE_PROGRAM % bedtools, "-i", "stdin"] + extra_params_list), input_str=bed_str)

def sort_and_merge_file(bed_path, bedtools=BEDTOOLS_PROGRAM):
    return merge_bed_str(sort_file(bed_path, bedtools=bedtools), bedtools=bedtools)

def sort_file(bed_path, extra_params_list=[], bedtools=BEDTOOLS_PROGRAM):
    return execute(" ".join([SORT_PROGRAM % bedtools, "-i", bed_path] + extra_params_list))

def merge_file(bed_path, extra_params_list=[], bedtools=BEDTOOLS_PROGRAM):
    return execute(" ".join([MERGE_PROGRAM % bedtools, "-i", bed_path] + extra_params_list))