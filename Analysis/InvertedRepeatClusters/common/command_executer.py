#!/usr/bin/env python3.8
"""
Function for wrapping subprocesses opening
"""

# region Imports
import subprocess
import traceback
import logging
import sys
# endregion


# START_STEP = "Process: %(ps_name)s; Started Step: %(step)s - %(name)s"
# RUNNING_STEP = "Process: %(ps_name)s; Running Step: %(step)s"
# ERROR_STEP_MSG = "Process: %(ps_name)s; Going To Error Step: %(step)s"
# END_STEP = "Process: %(ps_name)s; Finished Step: %(step)s"
# STARTED_RUNING = "Started Running"
# FILE_FLOW_START = "Started Processing File %(file)s"
# PS_END = "Process: %(ps_name)s; Exited With Code: %(ret_code)s"
# RAN_ON_MSG = "Ran (Or Running) On %s Files So Far"
# region Consts
# endregion


def execute(command, input_str=None):
    try:
        # info
        logging.debug(RUNNING_COMMAND % {'cmd': command})
        sys.stdout.flush()
        if input_str:
            p = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               input=input_str, universal_newlines=True, text=True)
        else:
            p = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
        # return output
        return p.stdout
    except subprocess.CalledProcessError as e:
        logging.error(ERROR_MSG % {'cmd': command, 'ret_code': e.returncode})
        # print("ERROR:",  p.stderr)
        print(traceback.format_exc())
        exit(1)


    # # clear stdout
    # sys.stdout.flush()
    # # intersect
    # p = subprocess.run(" ".join([INTERSECT_PROGRAM, "-a", repeat_bed_path, "-b", regions_bed_path] + extra_params_list), shell=True,
    #     stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    # print(p, p.stdout, p.stderr)
    # # stdout, stderr = p.stdout()
    # # check for error
    # if p.stderr:
    #     print("Could not intersect files: " + p.stderr)
    #     exit(1)
    # # split and strip output - return list of str
    # return p.stdout.strip().split("\n")
    # bedtools.2.27.1 sort -i stdin |bedtools.2.27.1 merge -i stdin | gzip
    # sys.stdout.flush()
    # p = subprocess.Popen(" ".join([SORT_PROGRAM, "-i", "stdin", "|", MERGE_PROGRAM, "-i", "stdin"]), shell=True,
    #                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    # out, err = p.communicate("\n".join(repeats_bed_all) + "\n")