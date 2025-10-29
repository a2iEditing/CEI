#!/usr/bin/env python3.8
"""
This module wraps some of the pybedtools objects to allow for working with BED12 in an ordered manner
"""
__author__ = 'Roni'

from pybedtools.cbedtools import create_interval_from_list
from pybedtools import BedTool

class BedToolTuple:
    def __init__(self, intervalA, intervalB, include_strand=False):
        self.intervalA = intervalA
        self.intervalB = intervalB
        self.include_strand = include_strand
        # self.original = bed12interval

    def __repr__(self):
        return str((self.intervalA, self.intervalB))

    @classmethod
    def bed_tool_tuple_from_bed12(cls, bed12interval, region_strand=False):
        return BedToolTuple(create_interval_from_list(bed12interval.fields[:6]),
                            create_interval_from_list(bed12interval.fields[6:]), 
                            region_strand)

    @classmethod
    def bed_tool_tuple_from_bed3_6(cls, bed3_6interval, region_strand=False):
        return BedToolTuple(create_interval_from_list(bed3_6interval.fields[:3]),
                            create_interval_from_list(bed3_6interval.fields[3:]),
                            region_strand) if not region_strand else BedToolTuple(create_interval_from_list(bed3_6interval.fields[:6]),
                            create_interval_from_list(bed3_6interval.fields[6:12]),
                            region_strand)

    # def bed12(self):
    #     return self.original

    def first(self):
        return self.intervalA

    def second(self):
        return self.intervalB

    def first_position(self):
        if self.include_strand: 
            return "\t".join(list(self.intervalA[:3]) + ["0", "0", self.intervalA.strand])
        return "\t".join(self.intervalA[:3])

    def second_position(self):
        if self.include_strand: 
            return "\t".join(list(self.intervalB[:3]) + ["0", "0", self.intervalB.strand])
        return "\t".join(self.intervalB[:3])

    def first_str(self):
        return "\t".join(self.intervalA) 

    def second_str(self):
        return "\t".join(self.intervalB)

    # def second_position_trimmed(self):
    #     """
    #     trim intervalB by intervalA limits
    #     :return:
    #     """
    #     # if intervalB is included within intervalA return it
    #     if self.intervalA.start <= self.intervalB.start and self.intervalA.end >= self.intervalB.end:
    #         return self.second_position()
    #     # if intervalB is to right of intervalA
    #     elif self.intervalA.start <= self.intervalB.start and self.intervalA.end <= self.intervalB.end:
    #         trimmed_interval = self.intervalB
    #         trimmed_interval.end = self.intervalA.end
    #         "\t".join(trimmed_interval[:3])
    #     # if intervalB is to left of intervalA
    #     elif self.intervalA.start >= self.intervalB.start and self.intervalA.end >= self.intervalB.end:
    #         trimmed_interval = self.intervalB
    #         trimmed_interval.start = self.intervalA.start
    #         "\t".join(trimmed_interval[:3])
    #     # if intervalB includes intervalA within it
    #     elif self.intervalA.start >= self.intervalB.start and self.intervalA.end <= self.intervalB.end:
    #         return self.first_position()
    #     else:
    #         print(self)
    #         raise Exception("Intervals not intersecting")


    # def overlap_interval(self):
    #
    #
    #     try:
    #         if not (self.intervalA.start <= self.intervalB.start and self.intervalA.stop >= self.intervalB.stop):
    #             a=BedToolTuple(self.intervalA,
    #                                 BedTool(self.second_position(), from_string=True).intersect(
    #                                     BedTool(self.first_position(), from_string=True)))
    #             return BedToolTuple(self.intervalA,
    #                                 BedTool(self.second_position(), from_string=True).intersect(
    #                                     BedTool(self.first_position(), from_string=True)))
    #     except NotImplementedError:
    #         return self


from enum import Enum


class Sign(str, Enum):
    PLUS = "+"
    MINUS = "-"

    @classmethod
    def sign(cls, sign_str):
        if sign_str == Sign.MINUS:
            return Sign.MINUS
        elif sign_str == Sign.PLUS:
            return Sign.PLUS
        raise ValueError("Unrecognized strand sign: %s" % sign_str)

class BedTupleList:
    @staticmethod
    def bed_file(tuple_list):
        return "\n".join([str(o.original).strip() for o in tuple_list])

    @staticmethod
    def bed_file_first(tuple_list):
        return "\n".join([str(o.intervalA).strip() for o in tuple_list])

    @staticmethod
    def bed_file_second(tuple_list):
        return "\n".join([str(o.intervalB).strip() for o in tuple_list])

    @staticmethod
    def bedtool_first(tuple_list):
        return BedTool(BedTupleList.bed_file_first(tuple_list), from_string=True)

    @staticmethod
    def bedtool_second(tuple_list):
        return BedTool(BedTupleList.bed_file_second(tuple_list), from_string=True)


    @staticmethod
    def str_list_to_bed(rec_str):
        return "\n".join(rec_str)

    @staticmethod
    def str_list_to_bedtool(rec_str):
        return BedTool(BedTupleList.str_list_to_bed(rec_str), from_string=True)