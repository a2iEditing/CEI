import csv
import os

def write_dict_to_csv(dict_list, ordered_keys_list, output_file, append = False):
    '''
    Writes list of dictionaries to CSV file
    :param dict_to_write: dictionary
    :param ordered_keys_list: wanted order of keys
    :param output_file: output path
    :param append: should the file be overwritten and not appended to?
    :return: None
    '''
    
    # if file exists and should be appended to
    if append and os.path.isfile(output_file):
        with open(output_file, 'a') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list)
            writer.writerows(dict_list)
    else:
        # else write file
        with open(output_file, 'w') as csvfile:
            # header treats group columns as one field, as that is the key : sums implementation
            # this also ensures our wanted order (to create a BED file)
            writer = csv.DictWriter(csvfile, lineterminator="\n", fieldnames=ordered_keys_list)
            # write header
            writer.writeheader()
            writer.writerows(dict_list)