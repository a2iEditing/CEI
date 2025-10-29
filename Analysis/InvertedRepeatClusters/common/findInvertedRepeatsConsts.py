# region Consts
# logging consts
LOG_DIR = os.path.join("%(out_dir)s", "InvertedRepeatsLogs",
                       "%(suffix)sFindInvertedRepeats_%(time)s.log")
RUN_PARAMS_DEBUG_MSG = "Starting Run with Params: %s"
INTERSECT_MSG = "Intersected regions with repeats, resulting in %s lines"
FAMILIES_FOUND_MSG = "Found %(num_families)s families to run on: %(families)s"
RUNNING_ON_FAMILY_MSG = "Running on %s repeats"
NUM_ELEMENTS_FOUND_PER_FAMILY_MSG = "Found %(num_elements)s %(type)s elements in %(num_regions)s regions for %(family)s repeats"
NUM_ELEMENTS_FOUND_MSG = "Found a total of %(num_elements)s elements in %(num_regions)s regions"
NUM_MERGED_ELEMENTS_FOUND_MSG = "Found a total of %(num_elements)s merged elements in %(num_regions)s merged regions"
SAVED_REPEATS_MSG = "Saved repeats BED file into %s"
SAVED_REGIONS_MSG = "Saved regions BED file into %s"
DONE_MSG = "DONE"
# output prefixes
REPEATS_PREFIX = "RepeatsIn%sOrientationAtRegions"
REGIONS_PREFIX = "RegionsWith%sOrientationRepeat"
OPPOSITE_REPEATS_FILE_PREFIX = REPEATS_PREFIX % "Inverted"
OPPOSITE_REGIONS_FILE_PREFIX = REGIONS_PREFIX % "Inverted"
CONTROL_REPEATS_FILE_PREFIX = REPEATS_PREFIX % "Tandem"
CONTROL_REGIONS_FILE_PREFIX = REGIONS_PREFIX % "Tandem"
# output suffixes
SORTED_BED_SUFFIX = ".sorted.bed"
MERGED_BED_SUFFIX = ".sorted.merged.bed"
COMPRESSED_BED_SUFFIX = ".gz"
# endregion

RepeatSet = Enum('RepeatSet', 'OPPOSITES CONTROLS')
