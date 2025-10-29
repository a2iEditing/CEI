#!/usr/bin/env python3.8
"""
This module attempts to wrap BED files and work with this data structure
"""


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

    # def __str__(self):
    #     return self.value

class BedRecord:
    def __init__(self, record_lst):
        self._chr = record_lst[0]
        self._start = record_lst[1]
        self._end = record_lst[2]
        self._name = record_lst[3]
        self._score = record_lst[4]
        self._strand = Sign.sign(record_lst[5])
        self._interval = "_".join([self.chr, self.start, self.end])

    def __repr__(self):
        return "<" + " ".join([self.chr, self.start, self.end, self.name, self.score, self.strand]) + ">"

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return self.chr == other.chr and self.start == other.start and self.end == other.start and \
               self.name == other.name and self.score == other.score and self.strand == other.strand

    # python-ic constructor overloading using class methods
    @classmethod
    def bed_record_from_str(cls, bed_str, bed_sep="\t"):
        return BedRecord(bed_str.split(bed_sep))

    # region Properties
    @property
    def chr(self):
        return self._chr
    @property
    def start(self):
        return self._start
    @property
    def end(self):
        return self._end
    @property
    def name(self):
        return self._name
    @property
    def score(self):
        return self._score
    @property
    def strand(self):
        return self._strand
    # endregion

    def bed(self):
        return "\t".join([self.chr, self.start, self.end, self.name, self.score, self.strand])

    def interval(self):
        return self._interval

    def sign(self):
        return self.strand

class BedTuple:
    def __init__(self, concat6bed_str, bed_sep="\t"):
        self.bedA = BedRecord(concat6bed_str.split(bed_sep)[:6])
        self.bedB = BedRecord(concat6bed_str.split(bed_sep)[6:])

    def __repr__(self):
        return str([self.bedA, self.bedB])

    def bed(self):
        return " | ".join([self.bedA.bed(), self.bedB.bed()])

    def first(self):
        return self.bedA

    def second(self):
        return self.bedB

class BedFile:
    def __init__(self, bed_record_list=None):
        self._bed_list = bed_record_list if bed_record_list else []

    # python-ic constructor overloading using class methods
    @classmethod
    def bed_file_from_str_list(cls, bed_str_list):
        return BedFile([BedRecord.bed_record_from_str(line) for line in bed_str_list])

    @property
    def bed_list(self):
        return self._bed_list

    def add_record(self, value):
        self._bed_list.append(value)

    def add_record_from_list(self, bed_list):
        self.add_record(BedRecord(bed_list))

    def add_record_from_str(self, bed_str):
        self.add_record(BedRecord.bed_record_from_str(bed_str))

    def add_record_list(self, rec_list):
        self._bed_list.extend(rec_list)

    def __add__(self, other):
        return BedFile(self.bed_list.extend(other.bed_list))

    def keep_chr(self, chromosome):
        return BedFile([b for b in self.bed_list if b.name == chromosome])

    def keep_name(self, name):
        return BedFile([b for b in self.bed_list if b.name == name])

    def keep_score(self, score):
        return BedFile([b for b in self.bed_list if b.score == score])

    def keep_strand(self, strand):
        return BedFile([b for b in self.bed_list if b.strand == strand])

    def bed_str_list(self):
        return [b.bed() for b in self.bed_list]

    def bed(self):
        return "\n".join(self.bed_str_list())

    # def sort(self):
    #     # TODO make sure this works... not sure
    #     self._bed_list.sort(key=lambda bed: (bed.chr, int(bed.start)))