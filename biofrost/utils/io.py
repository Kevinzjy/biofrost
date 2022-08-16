"""Functions for dealing with input/output files"""
import os
import gzip

__all__ = [
    "is_gzipped",
    "check_file",
    "check_dir",
]

class FileConflictError(Exception):
    """Directory name conflict with existed files"""
    pass


def is_gzipped(fname):
    """Determine whether file is gzipped

    From https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed

    Args:
        fname (str): Path of input file

    Returns:
        bool: True if file is gzipped, otherwise False
    """
    with gzip.open(fname, 'r') as f:
        try:
            f.read(1)
        except OSError:
            return False
    return True


def check_file(file_name):
    """Check is a file exist

    Args:
        file_name (str): query file

    Returns:
        file_path: absolute path of input file
    """
    if os.path.exists(file_name) and os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(f'{file_name} not found')


def check_dir(dir_name):
    """Check if a directory exist, create if not exist

    Args:
        dir_name (str): query directory

    Returns:
        dir_path: absolution path of directory
    """
    if os.path.exists(dir_name):
        if os.path.isdir(dir_name):
            pass
        else:
            raise FileConflictError(f'Directory: {dir_name} conflict with existed files')
    else:
        os.mkdir(dir_name)
    return os.path.abspath(dir_name)