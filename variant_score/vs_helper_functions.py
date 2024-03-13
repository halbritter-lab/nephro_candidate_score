from datetime import datetime
import os

def get_latest_file(prefix, directory, extension=None):
    """Return the latest file with given prefix in given directory. 
    If 'extension' is not specified, it will be recognized automatically, else, only filenames with the
    given extension are returned.
    'extension' can be specified in the format 'vcf', 'csv', 'vcf.gz', etc. without preceding '.'.
    """

    # list all files in the current directory
    files = os.listdir(directory)
    
    # list all files that match the prefix
    matching_files = [file for file in files if file.startswith(prefix)]

    # list only files with no additional string after the prefix (except date)
    exactly_matching_files = [i.split(prefix)[1] for i in matching_files if len(i.split(prefix)[1].split("_")) <= 2 and i.split(prefix)[1].split("_")[0] == ""]

    # exrtract dates
    date_strings = [i.split("_")[1].split(".")[0] for i in exactly_matching_files]

    if len(date_strings) == 0:
        raise FileNotFoundError(f"No files found with the prefix: {prefix}")

    # convert the date strings to datetime objects
    date_objects = [datetime.strptime(date, '%Y-%m-%d') for date in date_strings]

    # find the most recent date
    most_recent_date = max(date_objects)

    # convert most recent date to string
    most_recent_date_string = most_recent_date.strftime('%Y-%m-%d')
    # most_recent_date_string
    # file_dates = [datetime.strptime(i, date_formats) for i in date_strings]

    # get latest file with date
    latest_prefix_with_date = f"{prefix}_{most_recent_date_string}"

    # get all files with given prefix and latest date
    latest_files = [i for i in files if i.startswith(latest_prefix_with_date)]

    # get latest file matching prefix and extension
    if len(latest_files) == 1:
        latest_file = latest_files[0]
    elif len(latest_files) > 1 and extension == None:
        raise ValueError("There are multiple latest files with different extensions. Please specify extension.")
    elif len(latest_files) > 1 and extension is not None:
        latest_file = f"{latest_prefix_with_date}.{extension}"
    else:
        latest_file = None
        
    if extension is not None and not latest_file.endswith(extension):
        raise FileNotFoundError(f"No files found with given extension: .{extension}")

    return latest_file