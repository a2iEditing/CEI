#!/usr/bin/env python3.8
"""This module finds oppositely oriented repeats by given parameters.
Using pybedtools and wrappers
"""

# region Imports
import argparse
import errno
import inspect
import logging
import os
import sys
# for logging
from datetime import datetime
from pprint import pformat
from typing import List, Any, Iterable
# pybedtools

# use this for importing external self commands
if __name__ == '__main__':
    sys.path.append(os.path.dirname(__file__))
from common.init_logging_dict import init_logging_dict
from common.files_general import write_dict_to_csv
from common.PyBedToolsCommons import *
from findInvertedRepeatsConsts import *
import csv


# endregion

def flatten_recursive(lst: List[Any]) -> Iterable[Any]:
    """Flatten a list using recursion."""
    for item in lst:
        if isinstance(item, list):
            yield from flatten_recursive(item)
        else:
            yield item

def write_output(repeats_file_path_suffix, repeats_bed_all, found_regions_file_path, region_bed_all, closest_info, merge_output,
                 zip_output, set_type):
    """
    Writes BED files output, regions and repeats
    :param set_type: type of control files (name)
    :param repeats_file_path_suffix: repeats path prefix
    :param repeats_bed_all: repeat BedTool object
    :param found_regions_file_path: regions path prefix
    :param region_bed_all: regions BedToll object
    :param merge_output: should output be merged
    :param zip_output: should output bed zipped
    :return: bedTools regions, bedtool repeats
    """

    logging.debug(SAVING_MSG % set_type)

    repeats_bed_all = repeats_bed_all.sort()
    region_bed_all = region_bed_all.sort()

    logging.info(NUM_ELEMENTS_FOUND_MSG % {'num_elements': len(repeats_bed_all), 'num_regions': len(region_bed_all),
                                           'type': set_type})

    if merge_output:
        found_regions_file_path += MERGED_BED_SUFFIX + COMPRESSED_BED_SUFFIX if zip_output else MERGED_BED_SUFFIX
        repeats_file_path_suffix += MERGED_BED_SUFFIX + COMPRESSED_BED_SUFFIX if zip_output else MERGED_BED_SUFFIX
        region_bed_all = region_bed_all.merge().saveas(found_regions_file_path)
        repeats_bed_all = repeats_bed_all.merge().saveas(repeats_file_path_suffix)
        logging.info(
            NUM_MERGED_ELEMENTS_FOUND_MSG % {'num_elements': len(repeats_bed_all), 'num_regions': len(region_bed_all),
                                             'type': set_type})
    else:
        found_regions_file_path += SORTED_BED_SUFFIX + COMPRESSED_BED_SUFFIX if zip_output else SORTED_BED_SUFFIX
        repeats_file_path_suffix += SORTED_BED_SUFFIX + COMPRESSED_BED_SUFFIX if zip_output else SORTED_BED_SUFFIX
        region_bed_all = region_bed_all.saveas(found_regions_file_path)
        repeats_bed_all = repeats_bed_all.saveas(repeats_file_path_suffix)

    logging.debug(SAVED_REPEATS_MSG % {'file': repeats_file_path_suffix,
                                       'type': set_type})
    logging.debug(SAVED_REGIONS_MSG % {'file': found_regions_file_path,
                                       'type': set_type})
    
    if closest_info:
            write_dict_to_csv(dict_list=closest_info, ordered_keys_list=CLOSEST_INFO_HEADER, 
                              output_file=repeats_file_path_suffix + CLOSEST_CSV_SUFFIX)
            logging.debug(SAVED_CLOSEST_INFO_MSG % {'file': repeats_file_path_suffix + CLOSEST_CSV_SUFFIX,
                                       'type': set_type})
    
    return region_bed_all, repeats_bed_all


def calc_elements_sign_distribution():
    pass


def filter_clusters_by_element_length_range(bedtool_tuples_list, max_len_range=None, min_len_range=None):
    # get dict of regions
    region_dict = {pair.first_position(): set() for pair in bedtool_tuples_list}
    # associate full repeats (including strand) with each region
    [region_dict[pair.first_position()].add(pair.second_str()) for pair in bedtool_tuples_list]
    # make into BedTool object
    region_dict = {k: BedTupleList.str_list_to_bedtool(v) for k, v in region_dict.items()}
    # filter regions to include only those with an element within range
    if max_len_range and min_len_range:
        # if range is given, make sure there is at least one element within range in cluster
        regions_with_clusters_passing_limits = [k for k, v in region_dict.items() if
                                                any([min_len_range <= len(i) <= max_len_range for i in v])]
    elif min_len_range:
        # if bottom limit is given, make sure there is at least one element with under such length
        regions_with_clusters_passing_limits = [k for k, v in region_dict.items() if
                                                any([min_len_range <= len(i) for i in v])]
    elif max_len_range:
        # if top limit is given, make sure there is at least one element with over such length
        regions_with_clusters_passing_limits = [k for k, v in region_dict.items() if
                                                any([len(i) <= max_len_range for i in v])]
    # filter tuples and return
    return [pair for pair in bedtool_tuples_list if pair.first_position() in regions_with_clusters_passing_limits]


def max_distance_between_elements(bedtool_tuples_list, max_opp_distance=None, min_opp_distance=None, opp_dist_control=False):
    # get dict of regions
    region_dict = {pair.first_position(): set() for pair in bedtool_tuples_list}
    # associate full repeats (including strand) with each region
    [region_dict[pair.first_position()].add(pair.second_str()) for pair in bedtool_tuples_list]
    # make into BedTool object
    region_dict = {k: BedTupleList.str_list_to_bedtool(v) for k, v in region_dict.items()}
    # sort each bedtool
    closest = {k: v.sort() for k, v in region_dict.items()}
    # get minimal distance between 2 opposite elements
    if opp_dist_control:
        closest_info = [[{"Region" : k.split("\t")[0]+ ":" + k.split("\t")[1] + "-" + k.split("\t")[2], 
                      "PosRepeat": i[0] + ":" + i[1] + "-" + i[2], 
                      "StrandRepeat" : i[5],
                      "PosClosestRepeat": i[6] + ":" + i[7] + "-" + i[8], 
                      "StrandClosestRepeat" : i[11],
                      "RepeatPairDistance" : i[12]} for i in v.closest(v, d=True, s=True, io=True)] for k, v in closest.items()]
        # for control elements - find closes element
        closest = {k: min([int(i[12]) for i in v.closest(v, d=True, s=True, io=True)]) for k, v in closest.items()}
    else:
        closest_info = [[{"Region" : k.split("\t")[0]+ ":" + k.split("\t")[1] + "-" + k.split("\t")[2], 
                      "PosRepeat": i[0] + ":" + i[1] + "-" + i[2], 
                      "StrandRepeat" : i[5],
                      "PosClosestRepeat": i[6] + ":" + i[7] + "-" + i[8], 
                      "StrandClosestRepeat" : i[11],
                      "RepeatPairDistance" : i[12]} for i in v.closest(v, d=True, S=True)] for k, v in closest.items()]
        closest = {k: min([int(i[12]) for i in v.closest(v, d=True, S=True)]) for k, v in closest.items()}
    # filter regions to include only those with minimal distance up to required maximum
    if max_opp_distance and min_opp_distance:
        regions_with_min_distance_under_max = [k for k, v in closest.items() if
                                               min_opp_distance <= v < max_opp_distance]
    elif min_opp_distance:
        regions_with_min_distance_under_max = [k for k, v in closest.items() if min_opp_distance <= v]
    elif max_opp_distance:
        regions_with_min_distance_under_max = [k for k, v in closest.items() if v < max_opp_distance]
    # filter tuples and return
    return [pair for pair in bedtool_tuples_list if pair.first_position() in regions_with_min_distance_under_max], closest_info


def percent_of_repeat_bases_in_region():
    pass


def get_distribution():
    pass


def filter_by_num_elements_in_cluster(intersect_out_tuples, repeat_family, min_num_elements=None, max_num_elements=None):
    # add option of filtering first by number of participants within cluster
    # both minimum and maximum number of elements per cluster
    # use dictionary to count number of element per cluster and keep those within limits
    # do this using method - better
    # add option to also get statistics in the next part - number of sign + and - per cluster along with the elements
    # initialize dictionary - only for regions with a family repeat
    count_dict = {pair.first_position(): set() for pair in intersect_out_tuples if pair.second().score == repeat_family}
    [count_dict[pair.first_position()].add(pair.second_position()) for pair in intersect_out_tuples if
     pair.second().score == repeat_family]
    # # get signs on intersect - only if repeat belongs to family
    # for pair in intersect_out_tuples:
    #     if pair.second().score == repeat_family:
    #         count_dict[pair.first_position()] =+ 1
    if min_num_elements and max_num_elements:
        return [pair for pair in intersect_out_tuples if min_num_elements <= len(count_dict[pair.first_position()]) <= max_num_elements]
    elif min_num_elements and not max_num_elements:
        return [pair for pair in intersect_out_tuples if min_num_elements <= len(count_dict[pair.first_position()])]
    elif max_num_elements and not min_num_elements:
        return [pair for pair in intersect_out_tuples if len(count_dict[pair.first_position()]) <= max_num_elements]


def find_repeats_in_regions(out, regions_with_reps_bed, rf, keep_edges, maximal_opp_distance=None,
                            minimal_opp_distance=None, opp_dist_control=False,
                            min_len_in_cluster=None, max_len_in_cluster=None):
    """
    Gets tuples of region and repeat for given intersect output, family and for regions position list
    :param opp_dist_control:
    :param max_len_in_cluster:
    :param min_len_in_cluster:
    :param minimal_opp_distance:
    :param maximal_opp_distance:
    :param out: intersect output, first regions then repeats
    :param regions_with_reps_bed: lst(str) of region positions
    :param rf: repeat family string
    :param keep_edges: should repeats overflowing regions be trimmed to region limits or kept fully
    :return: list(str) of unique found repeats, bedtool tuples of found pairs
    """
    # get pairs
    pairs = [pair for pair in out if
             pair.first_position() in regions_with_reps_bed and pair.second().score == rf]

    # filter by certain element length within cluster
    if min_len_in_cluster or max_len_in_cluster:
        logging.debug(LENGTH_IN_CLUSTER_FILTER_MSG)
        pairs = filter_clusters_by_element_length_range(bedtool_tuples_list=pairs, max_len_range=max_len_in_cluster,
                                                        min_len_range=min_len_in_cluster)
        regions_with_reps_bed = set([pair.first_position() for pair in pairs])

    # filter by minimal gap between + and -
    if maximal_opp_distance or minimal_opp_distance:
        logging.debug(MAX_DISTANCE_FILTER_MSG)
        pairs, closest_info = max_distance_between_elements(bedtool_tuples_list=pairs, max_opp_distance=maximal_opp_distance,
                                              min_opp_distance=minimal_opp_distance, opp_dist_control=opp_dist_control)
        regions_with_reps_bed = set([pair.first_position() for pair in pairs])

    # get unique
    reps_pair = set([pair.second_position() for pair in pairs])

    # intersect to remove repeat "edges" sliding out of region
    if not keep_edges:
        logging.debug(REMOVE_EDGES_MSG)
        reps_pair = set(
            [str(i).strip() for i in BedTupleList.str_list_to_bedtool(reps_pair).intersect(
                BedTupleList.str_list_to_bedtool(regions_with_reps_bed))])
    return reps_pair, regions_with_reps_bed, closest_info if maximal_opp_distance or minimal_opp_distance else None


def get_opposites_and_controls(all_region_bed_str_lst, all_repeats_bed_str_lst, intersect_out_tuples, keep_edges,
                               reduce_discrepancies, repeat_family, closest_info, opposites_key, controls_key,
                               maximal_opp_distance=None, minimal_opp_distance=None, apply_opp_distance_to_control=False,
                               min_len_in_cluster=None, max_len_in_cluster=None):
    """
    Splits an intersect BedTool object into regions with opposite signs and regions with same sign (and their repeats)
    :param apply_opp_distance_to_control: should cluster distance be also calculated for control
    :param min_len_in_cluster:
    :param max_len_in_cluster:
    :param minimal_opp_distance:
    :param maximal_opp_distance:
    :param reduce_discrepancies:
    :param all_region_bed_str_lst: list to add found regions
    :param all_repeats_bed_str_lst: list to add found repeats
    :param intersect_out_tuples: intersect output, first regions then repeats
    :param keep_edges: should repeats overflowing regions be trimmed to region limits or kept fully
    :param repeat_family: repeat family string
    :param opposites_key: RepeatSet enum of wanted opposites group
    :param controls_key: RepeatSet enum of wanted opposites group
    :return: set(str) of unique opposite and control regions
    """

    # initialize dictionary - only for regions with a family repeat
    sign_dict = {pair.first_position(): set() for pair in intersect_out_tuples if pair.second().score == repeat_family}
    # get signs on intersect - only if repeat belongs to family
    [sign_dict[pair.first_position()].add(pair.second().strand) for pair in intersect_out_tuples if
     pair.second().score == repeat_family]

    # region Opposites
    # get found regions with elements with 2 different signs (unique)
    reps_with_opp_pair = ""
    regions_with_opp_reps_bed = set(
        [k for k in sign_dict.keys() if Sign.PLUS in sign_dict[k] and Sign.MINUS in sign_dict[k]])
    if regions_with_opp_reps_bed:
        # find repeats and add to general list
        reps_with_opp_pair, regions_with_opp_reps_bed, regions_with_opp_reps_closest_info = find_repeats_in_regions(out=intersect_out_tuples,
                                                                                                                    regions_with_reps_bed=regions_with_opp_reps_bed,
                                                                                                                    rf=repeat_family, keep_edges=keep_edges,
                                                                                                                    maximal_opp_distance=maximal_opp_distance,
                                                                                                                    minimal_opp_distance=minimal_opp_distance,
                                                                                                                    min_len_in_cluster=min_len_in_cluster,
                                                                                                                    max_len_in_cluster=max_len_in_cluster)
        # add repeats and regions into overall list
        all_repeats_bed_str_lst[opposites_key].extend(reps_with_opp_pair)
        all_region_bed_str_lst[opposites_key].extend(regions_with_opp_reps_bed)
        if closest_info and regions_with_opp_reps_closest_info:
            closest_info[opposites_key].extend(regions_with_opp_reps_closest_info)
        # notify
        logging.info(
            NUM_ELEMENTS_FOUND_PER_FAMILY_MSG % {'family': repeat_family, 'num_elements': len(reps_with_opp_pair),
                                                 'num_regions': len(regions_with_opp_reps_bed),
                                                 'type': opposites_key.name})
    # endregion

    # region Controls
    # get found regions with elements with 2 different signs (unique)
    regions_with_ctrl_reps_bed = set(
        [k for k in sign_dict.keys() if (Sign.PLUS in sign_dict[k] and Sign.MINUS not in sign_dict[k]) or (
                Sign.MINUS in sign_dict[k] and Sign.PLUS not in sign_dict[k])])
    if regions_with_ctrl_reps_bed:
        # find repeats and add to general list
        # no need to check for distance from opposite - this is the controls
        # does check for element length by range in cluster
        if apply_opp_distance_to_control:
            reps_with_no_opp_pair, regions_with_ctrl_reps_bed, regions_with_ctrl_reps_closest_info = find_repeats_in_regions(out=intersect_out_tuples,
                                                                                        regions_with_reps_bed=regions_with_ctrl_reps_bed,
                                                                                        rf=repeat_family,
                                                                                        keep_edges=keep_edges,
                                                                                        maximal_opp_distance=maximal_opp_distance,
                                                                                        minimal_opp_distance=minimal_opp_distance,
                                                                                        opp_dist_control=True,
                                                                                        min_len_in_cluster=min_len_in_cluster,
                                                                                        max_len_in_cluster=max_len_in_cluster
                                                                                        )
        else:
            reps_with_no_opp_pair, regions_with_ctrl_reps_bed, regions_with_ctrl_reps_closest_info = find_repeats_in_regions(out=intersect_out_tuples,
                                                                                                                            regions_with_reps_bed=regions_with_ctrl_reps_bed,
                                                                                                                            rf=repeat_family,
                                                                                                                            keep_edges=keep_edges,
                                                                                                                            min_len_in_cluster=min_len_in_cluster,
                                                                                                                            max_len_in_cluster=max_len_in_cluster
                                                                                                                            )

        # remove repeats that have an opposite using some other region (overlapping genes or isoforms)
        if reduce_discrepancies:
            reps_with_no_opp_pair = {r for r in reps_with_no_opp_pair if r not in reps_with_opp_pair}
        # add regions and repeats into overall list
        all_region_bed_str_lst[controls_key].extend(regions_with_ctrl_reps_bed)
        all_repeats_bed_str_lst[controls_key].extend(reps_with_no_opp_pair)
        if closest_info and regions_with_ctrl_reps_closest_info:
            closest_info[controls_key].extend(regions_with_ctrl_reps_closest_info)
        # notify
        logging.info(
            NUM_ELEMENTS_FOUND_PER_FAMILY_MSG % {'family': repeat_family, 'num_elements': len(reps_with_no_opp_pair),
                                                 'num_regions': len(regions_with_ctrl_reps_bed),
                                                 'type': controls_key.name})
    # endregion
    return regions_with_opp_reps_bed, regions_with_ctrl_reps_bed


# noinspection DuplicatedCode
def find_repeat_sets_per_region(intersect_out_tuples, all_region_bed_str_lst, all_repeats_bed_str_lst, repeat_family,
                                keep_edges, repeat_elements_bedtool, closest_info, control_window_size=None, region_strand=False, window_size_left_upstream=None, window_size_right_downstream=None,
                                reduce_discrepancies=False, maximal_opp_distance=None, minimal_opp_distance=None,
                                min_len_in_cluster=None, max_len_in_cluster=None, max_num_elements=None, min_num_elements=None):
    """
    Finds regions in which given family has a set of repeats where at least one is of the opposite sign to others
    :param min_num_elements: filter clusters with less or more than the given number of elements
    :param max_num_elements:
    :param max_len_in_cluster:
    :param min_len_in_cluster:
    :param minimal_opp_distance:
    :param maximal_opp_distance:
    :param reduce_discrepancies:
    :param control_window_size:
    :param repeat_elements_bedtool:
    :param intersect_out_tuples: intersect output, first regions then repeats
    :param all_region_bed_str_lst: list to add found regions
    :param all_repeats_bed_str_lst: list to add found repeats
    :param repeat_family: repeat family string
    :param keep_edges: should repeats overflowing regions be trimmed to region limits or kept fully
    :return: None
    """
    logging.debug(RUNNING_ON_FAMILY_MSG % repeat_family)

    # filtering first by both minimum and maximum number of elements per cluster
    # TODO add option to also get statistics in the next part - number of sign + and - per cluster along with the elements
    if min_num_elements or max_num_elements:
        intersect_out_tuples = filter_by_num_elements_in_cluster(intersect_out_tuples, repeat_family=repeat_family,
                                                                 max_num_elements=max_num_elements, min_num_elements=min_num_elements)

    # get opposites and control regions
    _, regions_with_ctrl_reps_bed = get_opposites_and_controls(all_region_bed_str_lst=all_region_bed_str_lst,
                                                               all_repeats_bed_str_lst=all_repeats_bed_str_lst,
                                                               intersect_out_tuples=intersect_out_tuples,
                                                               keep_edges=keep_edges,
                                                               maximal_opp_distance=maximal_opp_distance,
                                                               minimal_opp_distance=minimal_opp_distance,
                                                               # only if no window is considered and opp distance exists
                                                               apply_opp_distance_to_control=True if control_window_size <= 0 and (
                                                                           maximal_opp_distance or minimal_opp_distance) else False,
                                                               min_len_in_cluster=min_len_in_cluster,
                                                               max_len_in_cluster=max_len_in_cluster,
                                                               reduce_discrepancies=reduce_discrepancies,
                                                               repeat_family=repeat_family,
                                                               opposites_key=RepeatSet.OPPOSITES,
                                                               controls_key=RepeatSet.CONTROLS,
                                                               closest_info = closest_info)
    
    # split control repeats with one sign only into:
    # 1. control repeats with other sign in surrounding window and
    # 2. control repeats without other sign in surrounding window
    # Do this only if window is greater than 0    
    if control_window_size > 0 or window_size_left_upstream > 0 or window_size_right_downstream > 0:
        if control_window_size > 0:
            # get BedTool object of regions and intersect with repeats around given window
            ctrl_window_intersect_out = BedTupleList.str_list_to_bedtool(regions_with_ctrl_reps_bed).window(
                repeat_elements_bedtool, w=control_window_size)
            logging.info(INTERSECT_WINDOW_CONTROL_MSG % {'type': RepeatSet.CONTROLS, 'bp': control_window_size,
                                                        'lines': len(ctrl_window_intersect_out)})
        elif region_strand:
            ctrl_window_intersect_out = BedTupleList.str_list_to_bedtool(regions_with_ctrl_reps_bed).window(
                repeat_elements_bedtool, l=window_size_left_upstream, r=window_size_right_downstream, sw=True)
            #strand %(strand)s considered
            logging.info(INTERSECT_ASYMMETRICAL_WINDOW_CONTROL_MSG % {'type': RepeatSet.CONTROLS, 
                                                                      'lbp': window_size_left_upstream,
                                                                      'rbp' : window_size_right_downstream,
                                                                      'strand' : 'IS',
                                                        'lines': len(ctrl_window_intersect_out)})
        else:
            ctrl_window_intersect_out = BedTupleList.str_list_to_bedtool(regions_with_ctrl_reps_bed).window(
                repeat_elements_bedtool, l=window_size_left_upstream, r=window_size_right_downstream)
            logging.info(INTERSECT_ASYMMETRICAL_WINDOW_CONTROL_MSG % {'type': RepeatSet.CONTROLS, 
                                                                      'lbp': window_size_left_upstream,
                                                                      'rbp' : window_size_right_downstream,
                                                                      'strand' : 'NOT',
                                                        'lines': len(ctrl_window_intersect_out)})
        # create BedTool Tuple
        # print(ctrl_window_intersect_out)
        ctrl_window_intersect_out = [BedToolTuple.bed_tool_tuple_from_bed3_6(interval, region_strand) for interval in ctrl_window_intersect_out]
        
        # get opposites and control regions in control window
        get_opposites_and_controls(all_region_bed_str_lst=all_region_bed_str_lst,
                                   all_repeats_bed_str_lst=all_repeats_bed_str_lst,
                                   intersect_out_tuples=ctrl_window_intersect_out, keep_edges=keep_edges,
                                   maximal_opp_distance=maximal_opp_distance, minimal_opp_distance=minimal_opp_distance,
                                   apply_opp_distance_to_control=True,
                                   reduce_discrepancies=reduce_discrepancies,
                                   repeat_family=repeat_family, opposites_key=RepeatSet.CONTROLS_WITH_WINDOW_OPPOSITE,
                                   controls_key=RepeatSet.CONTROLS_NO_WINDOW_OPPOSITE,
                                   closest_info = closest_info)

    # TODO - add statistics filters such as minimal gap between + and -,
    #  and sign distribution
    #  Question: why do some regions appear both as "same in window" and "opposites in window" [Answer: isoforms]


# noinspection DuplicatedCode
def main(repeat_elements_bed, genomic_regions_file, output_dir, output_files_prefix, merge_output, zip_output, log_path,
         keep_edges, window_size, region_strand, window_size_left_upstream, window_size_right_downstream, 
         maximal_opp_distance, minimal_opp_distance, reduce_discrepancies, min_len, max_len,
         min_len_in_cluster, max_len_in_cluster, min_num_elements, max_num_elements, output_mode):
    """
    Runs program flow over repeat families
    :param min_num_elements: limit number of elements per cluster
    :param max_num_elements: limit number of elements per cluster
    :param minimal_opp_distance: limit minimal distance of elements in cluster
    :param min_len_in_cluster:
    :param max_len_in_cluster:
    :param output_mode:
    :param min_len:
    :param max_len:
    :param maximal_opp_distance: limit minimal distance of elements in cluster
    :param reduce_discrepancies:
    :param window_size:
    :param repeat_elements_bed: repeats file path
    :param genomic_regions_file: regions file path
    :param output_dir:
    :param output_files_prefix: prefix to add to filenames if wanted
    :param merge_output: should output be merged
    :param zip_output: should putput be zipped
    :param log_path: logging file path
    :param keep_edges: should repeat edges be trimmed to region
    :return: None
    """
    # init logging
    init_logging_dict(log_path)
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    logging.debug(RUN_PARAMS_DEBUG_MSG % pformat(dict([(i, values[i]) for i in args])))

    # create bedtool element
    repeat_elements_bed = BedTool(repeat_elements_bed)
    # filter by minimal length if required
    if min_len or max_len:
        logging.debug(FILTER_READ_LEN_MSG % (min_len, max_len))
    if min_len and max_len:
        repeat_elements_bed = repeat_elements_bed.filter(lambda x: min_len < len(x) <= max_len).saveas()
    elif min_len:
        repeat_elements_bed = repeat_elements_bed.filter(lambda x: len(x) > min_len).saveas()
    elif max_len:
        repeat_elements_bed = repeat_elements_bed.filter(lambda x: len(x) <= max_len).saveas()

    # run first intersect command
    out = BedTool(genomic_regions_file).intersect(repeat_elements_bed, wa=True, wb=True)
    
    logging.info(INTERSECT_MSG % len(out))
    # Convert BedTool objects to tuples and export to CSV
    out_tuples = [BedToolTuple.bed_tool_tuple_from_bed12(interval, region_strand) for interval in out]

    # Export out_tuples to CSV
    csv_output_path = os.path.join(output_dir, output_files_prefix + "intersection_results.csv")
    try:
        with open(csv_output_path, 'w', newline='') as csvfile:
            # Define CSV headers based on BedToolTuple structure
            csv_headers = ['region', 'region_name', 'region_score', 'region_strand', 'region_length',
                          'repeat', 'repeat_name', 'repeat_score', 'repeat_strand', 'repeat_length', 'repeat_cropped']
            writer = csv.DictWriter(csvfile, fieldnames=csv_headers)
            writer.writeheader()
            
            # Write each BedToolTuple as a row in the CSV
            for tup in out_tuples:
                first = tup.first()
                second = tup.second()
                # Calculate the cropped repeat coordinates
                repeat_cropped_start = max(first.start, second.start)
                repeat_cropped_end = min(first.end, second.end)
                repeat_cropped = f"{second.chrom}:{repeat_cropped_start}-{repeat_cropped_end}"
                
                writer.writerow({
                    'region': f"{first.chrom}:{first.start}-{first.end}",
                    'region_name': getattr(first, 'name', '.'),
                    'region_score': getattr(first, 'score', '.'),
                    'region_strand': first.strand,
                    'region_length': first.end - first.start,
                    'repeat': f"{second.chrom}:{second.start}-{second.end}",
                    'repeat_name': getattr(second, 'name', '.'),
                    'repeat_score': second.score,
                    'repeat_strand': second.strand,
                    'repeat_length': second.end - second.start,
                    'repeat_cropped': repeat_cropped
                })
        logging.info(f"Exported intersection results to CSV: {csv_output_path}")
    except Exception as e:
        logging.error(f"Failed to export intersection results to CSV: {str(e)}")

    # get repeat families
    repeat_family_set = set()
    [repeat_family_set.add(pair.second().score) for pair in out_tuples]
    logging.info(FAMILIES_FOUND_MSG % {'num_families': len(repeat_family_set), 'families': repeat_family_set})

    # do this for each family - get lists of bed records in str format (for each set of wanted region types)
    region_bed_str_lst = {s: [] for s in RepeatSet}
    repeats_bed_str_lst = {s: [] for s in RepeatSet}
    closest_info = {s: [] for s in RepeatSet} if maximal_opp_distance or minimal_opp_distance else None
    for rf in repeat_family_set:
        find_repeat_sets_per_region(intersect_out_tuples=out_tuples, all_region_bed_str_lst=region_bed_str_lst,
                                    all_repeats_bed_str_lst=repeats_bed_str_lst, repeat_family=rf,
                                    keep_edges=keep_edges, repeat_elements_bedtool=repeat_elements_bed,
                                    control_window_size=window_size, region_strand=region_strand, window_size_left_upstream=window_size_left_upstream, window_size_right_downstream=window_size_right_downstream,
                                    maximal_opp_distance=maximal_opp_distance, minimal_opp_distance=minimal_opp_distance, min_len_in_cluster=min_len_in_cluster,max_len_in_cluster=max_len_in_cluster, reduce_discrepancies=reduce_discrepancies,
                                    min_num_elements=min_num_elements,max_num_elements=max_num_elements,
                                    closest_info = closest_info)
    
    # region OutputPrep
    # create output file directory if does not exist
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
            
    # opposite orientation
    if OutputModeOptions.All in output_mode or OutputModeOptions.Inverted in output_mode:
        write_output(
            repeats_file_path_suffix=os.path.join(output_dir, output_files_prefix + OPPOSITE_REPEATS_FILE_PREFIX),
            repeats_bed_all=BedTupleList.str_list_to_bedtool(set(repeats_bed_str_lst[RepeatSet.OPPOSITES])),
            found_regions_file_path=os.path.join(output_dir, output_files_prefix + OPPOSITE_REGIONS_FILE_PREFIX),
            region_bed_all=BedTupleList.str_list_to_bedtool(set(region_bed_str_lst[RepeatSet.OPPOSITES])),
            # closest_info=[k for i in closest_info[RepeatSet.OPPOSITES] for k in i],
            closest_info=list(flatten_recursive(closest_info[RepeatSet.OPPOSITES])) if closest_info else None,
            merge_output=merge_output, zip_output=zip_output, set_type=RepeatSet.OPPOSITES.name)
    # same orientations (control)
    if OutputModeOptions.All in output_mode or OutputModeOptions.Tandem in output_mode:
        write_output(
            repeats_file_path_suffix=os.path.join(output_dir, output_files_prefix + CONTROL_REPEATS_FILE_PREFIX),
            repeats_bed_all=BedTupleList.str_list_to_bedtool(set(repeats_bed_str_lst[RepeatSet.CONTROLS])),
            found_regions_file_path=os.path.join(output_dir, output_files_prefix + CONTROL_REGIONS_FILE_PREFIX),
            region_bed_all=BedTupleList.str_list_to_bedtool(set(region_bed_str_lst[RepeatSet.CONTROLS])),
            closest_info=list(flatten_recursive(closest_info[RepeatSet.CONTROLS])) if closest_info else None,
            merge_output=merge_output, zip_output=zip_output, set_type=RepeatSet.CONTROLS.name)
    # same orientation, with other opposite repeat in window
    if OutputModeOptions.All in output_mode or OutputModeOptions.WindowInverted in output_mode:
        write_output(
            repeats_file_path_suffix=os.path.join(output_dir,
                                                  output_files_prefix + OPPOSITE_IN_WINDOW_REPEATS_FILE_PREFIX),
            repeats_bed_all=BedTupleList.str_list_to_bedtool(
                set(repeats_bed_str_lst[RepeatSet.CONTROLS_WITH_WINDOW_OPPOSITE])),
            found_regions_file_path=os.path.join(output_dir,
                                                 output_files_prefix + OPPOSITE_IN_WINDOW_REGIONS_FILE_PREFIX),
            region_bed_all=BedTupleList.str_list_to_bedtool(
                set(region_bed_str_lst[RepeatSet.CONTROLS_WITH_WINDOW_OPPOSITE])),
            closest_info=list(flatten_recursive(closest_info[RepeatSet.CONTROLS_WITH_WINDOW_OPPOSITE])) if closest_info else None,
            merge_output=merge_output, zip_output=zip_output, set_type=RepeatSet.CONTROLS_WITH_WINDOW_OPPOSITE.name)
    # same orientations (control)
    if OutputModeOptions.All in output_mode or OutputModeOptions.WindowTandem in output_mode:
        write_output(
            repeats_file_path_suffix=os.path.join(output_dir,
                                                  output_files_prefix + CONTROL_IN_WINDOW_REPEATS_FILE_PREFIX),
            repeats_bed_all=BedTupleList.str_list_to_bedtool(
                set(repeats_bed_str_lst[RepeatSet.CONTROLS_NO_WINDOW_OPPOSITE])),
            found_regions_file_path=os.path.join(output_dir,
                                                 output_files_prefix + CONTROL_IN_WINDOW_REGIONS_FILE_PREFIX),
            region_bed_all=BedTupleList.str_list_to_bedtool(set(region_bed_str_lst[RepeatSet.CONTROLS_NO_WINDOW_OPPOSITE])),
            closest_info=list(flatten_recursive(closest_info[RepeatSet.CONTROLS_NO_WINDOW_OPPOSITE])) if closest_info else None,
            merge_output=merge_output, zip_output=zip_output, set_type=RepeatSet.CONTROLS_NO_WINDOW_OPPOSITE.name)

    # endregion

    logging.debug(DONE_MSG)


if __name__ == '__main__':
    # Create parser
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='Find repeats that are in opposite orientation at given regions, '
                                                 'by families.')
    parser.add_argument('-r', '--repeat_input', dest='repeat_input_file', action='store', metavar='repeat_input_file',
                        nargs='?', required=True,
                        help='Path to repeat BED file: chr, start, end, (score), family, strand')
    parser.add_argument('-i', '--regions_input', dest='regions_input_file', action='store',
                        metavar='regions_input_file', nargs='?', required=True,
                        help='Path to wanted regions BED file: chr, start, end, (name), (score), strand')
    parser.add_argument('-o', '--output', dest='output_dir', action='store', metavar='output_dir', nargs='?',
                        required=True, help='Path to output directory of BED files, will be created if does not exist')
    parser.add_argument('-p', '--prefix', dest='output_files_prefix', action='store', metavar='output_files_prefix',
                        default="", help="Add prefix to output files' names")
    parser.add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path',
                        default="", help="The path where the log will be written")
    parser.add_argument('--min_rep_len', dest='min_rep_len', action='store', metavar='min_rep_len', type=int,
                        default=None, help="Filter by given minimal repeat length before generating pairs (exclusive)")
    parser.add_argument('--max_rep_len', dest='max_rep_len', action='store', metavar='max_rep_len', type=int,
                        default=None, help="Filter by given maximal repeat length before generating pairs (inclusive)")
    parser.add_argument('-w', '--window_size', dest='window_size', action='store', metavar='window_size', type=int,
                        default=0, help="Window size in which to look for opposite pair for control elements (symmetrical)")
    parser.add_argument('-wl', '--window_size_left', dest='window_size_left', action='store', metavar='window_left', type=int,
                        default=0, help="Window size left in which to look for opposite pair for control elements (asymmetrical of right). Combine with -s to look for upstream")
    parser.add_argument('-wr', '--window_size_right', dest='window_size_right', action='store', metavar='window_right', type=int,
                        default=0, help="Window size right in which to look for opposite pair for control elements (asymmetrical of left). Combine with -s for downstream")
    parser.add_argument('--keep_edges', dest='keep_edges', action='store_true', default=False,
                        help='Allows repeats to extend beyond regions, if such a case occurs. Required for getting also the repeats in window outside of region in output.')
    parser.add_argument('-s', '--strand', dest='output_strand', action='store_true', default=False,
                        help='Find inverted repeats when considering strand (that is, a position with + is not equal to a position with -, both in repeats and regions). Consider strand when using window (position A on + strand and on - strand will be considered seperately)')
    parser.add_argument('--min_cluster_size', dest='min_cluster_size', action='store', metavar='min_cluster_size', type=int,
                        default=None, help="Keep only clusters that contain at least this number of elements (inclusive)")
    parser.add_argument('--max_cluster_size', dest='max_cluster_size', action='store', metavar='max_cluster_size', type=int,
                        default=None, help="Keep only clusters that contain at most this number of elements (inclusive)")
    parser.add_argument('-d', '--maximal_distance', dest='maximal_distance', action='store', metavar='maximal_distance',
                        type=int, default=None,
                        help="Maximal allowed distance between opposite elements. Minimal possible distance is considered, including overlap (exclusive). "
                             "For same orientation clusters, the closest element is considered without regards to orientation.")
    parser.add_argument('--minimal_distance', dest='minimal_distance', action='store', metavar='minimal_distance',
                        type=int, default=None,
                        help="Minimal allowed distance between opposite elements. Minimal possible distance is considered, including overlap (inclusive). "
                             "For same orientation clusters, the closest element is considered without regards to orientation. USE 0 TO OUTPUT FOR EACH  PAIR DISTANCE INFORMATION.")
    parser.add_argument('--minimal_length_within_cluster', dest='min_len_in_cluster', action='store', metavar='min_len_in_cluster',
                        type=int, default=None, help="Require each cluster to contain an element at least this long (inclusive)")
    parser.add_argument('--maximal_length_within_cluster', dest='max_len_in_cluster', action='store', metavar='max_len_in_cluster',
                        type=int, default=None, help="Require each cluster to contain an element at up to this length (inclusive)")
    parser.add_argument('--reduce_discrepancies', dest='reduce_discrepancies', action='store_true', default=False,
                        help='Repeats that do not have a pair in one region but do have an opposite pair in some another '
                             'region will only be considered as oppositely oriented (will not be included in control)')
    parser.add_argument('--merge', dest='merge_output_file', action='store_true', default=False,
                        help='Should found repeats output file be merged?')
    parser.add_argument('--zip', dest='zip_output_file', action='store_true', default=False,
                        help='Should output files be gzipped?')
    parser.add_argument('--output_mode', dest='output_mode', choices=OutputModeOptionsDict.keys(), nargs='+',
                        default=[OutputModeOptions.All.name],
                        help='Which files to output (will not affect calculation speed)')

    options = parser.parse_args()
    
    assert not options.window_size or (not options.window_size_left and not options.window_size_right), "Cannot define both symmetrical (-w) and asymmetrical windows (-wl, -wr)"

    main(repeat_elements_bed=options.repeat_input_file, genomic_regions_file=options.regions_input_file,
         output_dir=options.output_dir, output_files_prefix=options.output_files_prefix,
         merge_output=options.merge_output_file, zip_output=options.zip_output_file,
         log_path=options.log_path if options.log_path else LOG_DIR % {'out_dir': options.output_dir,
                                                                       'suffix': options.output_files_prefix,
                                                                       'time': datetime.today().isoformat()},
         keep_edges=options.keep_edges, window_size=options.window_size, region_strand=options.output_strand, 
         window_size_left_upstream = options.window_size_left, window_size_right_downstream = options.window_size_right,
         maximal_opp_distance=options.maximal_distance, minimal_opp_distance=options.minimal_distance,
         min_len_in_cluster=options.min_len_in_cluster, max_len_in_cluster=options.max_len_in_cluster,
         reduce_discrepancies=options.reduce_discrepancies,
         min_len=options.min_rep_len, max_len=options.max_rep_len,
         min_num_elements=options.min_cluster_size,max_num_elements=options.max_cluster_size,
         output_mode=[OutputModeOptionsDict[i] for i in options.output_mode])
