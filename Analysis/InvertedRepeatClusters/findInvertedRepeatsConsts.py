# region Imports
import os
from enum import Enum

# endregion

# region Consts
# region Logging
LOG_DIR = os.path.join("%(out_dir)s", "InvertedRepeatsLogs",
                       "%(suffix)sFindInvertedRepeats_%(time)s.log")
RUN_PARAMS_DEBUG_MSG = "Starting Run with Params: %s"
FILTER_READ_LEN_MSG = "Filtering repeats to include only repeats within %s to %s length"
INTERSECT_MSG = "Intersected regions with repeats, resulting in %s lines"
INTERSECT_WINDOW_CONTROL_MSG = "Intersected %(type)s regions with repeats around %(bp)s bp window, resulting in %(lines)s lines"
INTERSECT_ASYMMETRICAL_WINDOW_CONTROL_MSG = "Intersected %(type)s regions with repeats around %(lbp)s bp window upstream and %(rbp)s downstream, strand %(strand)s considered, resulting in %(lines)s lines"
FAMILIES_FOUND_MSG = "Found %(num_families)s families to run on: %(families)s"
RUNNING_ON_FAMILY_MSG = "Running on %s repeats"
SAVING_MSG = "Saving %s"
REMOVE_EDGES_MSG = "Removing repeats edges that extend out of region"
MAX_DISTANCE_FILTER_MSG = "Removing regions exceeding maximal or minimal distance allowed between closest opposite/control repeats"
LENGTH_IN_CLUSTER_FILTER_MSG = "Removing clusters not including an element within the length range specified"
NUM_ELEMENTS_FOUND_PER_FAMILY_MSG = "Found %(num_elements)s %(type)s elements in %(num_regions)s regions for %(family)s repeats"
NUM_ELEMENTS_FOUND_MSG = "Found a total of %(num_elements)s %(type)s elements in %(num_regions)s regions"
NUM_MERGED_ELEMENTS_FOUND_MSG = "Found a total of %(num_elements)s merged %(type)s elements in %(num_regions)s merged regions"
SAVED_REPEATS_MSG = "Saved %(type)s repeats BED file into %(file)s"
SAVED_REGIONS_MSG = "Saved %(type)s regions BED file into %(file)s"
SAVED_CLOSEST_INFO_MSG = "Saved  %(type)s closest pair information into %(file)s"
DONE_MSG = "DONE"
# endregion
# region Output
# region Prefixes
REPEATS_PREFIX = "RepeatsIn%sOrientationAtRegions"
REGIONS_PREFIX = "RegionsWith%sOrientationRepeat"
OPPOSITE_REPEATS_FILE_PREFIX = REPEATS_PREFIX % "Inverted"
OPPOSITE_REGIONS_FILE_PREFIX = REGIONS_PREFIX % "Inverted"
CONTROL_REPEATS_FILE_PREFIX = REPEATS_PREFIX % "Tandem"
CONTROL_REGIONS_FILE_PREFIX = REGIONS_PREFIX % "Tandem"
OPPOSITE_IN_WINDOW_REPEATS_FILE_PREFIX = REPEATS_PREFIX % "WindowInverted"
OPPOSITE_IN_WINDOW_REGIONS_FILE_PREFIX = REGIONS_PREFIX % "WindowInverted"
CONTROL_IN_WINDOW_REPEATS_FILE_PREFIX = REPEATS_PREFIX % "WindowTandem"
CONTROL_IN_WINDOW_REGIONS_FILE_PREFIX = REGIONS_PREFIX % "WindowTandem"
# endregion
# region Suffixes
SORTED_BED_SUFFIX = ".sorted.bed"
MERGED_BED_SUFFIX = ".sorted.merged.bed"
COMPRESSED_BED_SUFFIX = ".gz"
CLOSEST_CSV_SUFFIX = ".ClosestDistance.csv"
# endregion
# endregion
# endregion

# region Enums
RepeatSet = Enum('RepeatSet', 'OPPOSITES CONTROLS CONTROLS_WITH_WINDOW_OPPOSITE CONTROLS_NO_WINDOW_OPPOSITE')
OutputModeOptions = Enum('OutputModeOptions', 'All Inverted WindowInverted WindowTandem Tandem')
OutputModeOptionsDict = {i.name : i for i in OutputModeOptions}
# endregion

# region Headers
CLOSEST_INFO_HEADER = ["Region", "PosRepeat", "StrandRepeat", "PosClosestRepeat", "StrandClosestRepeat", "RepeatPairDistance"]
# endregion
