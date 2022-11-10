import os
import sys
import time
import logging
from multiprocessing.managers import BaseManager


class ProgressBar(object):
    """A simple progress bar with date stamp"""

    def __init__(self, width=50, total=100):
        """Init the ProgressBar object

        Paramters
        ---------
        width : integer, optional (default=50)
            The width of progress bar
        total : integer, optional (default=100)
            The total jobs of progress bar
        """
        self.last_x = -1
        self.width = width
        self.total = total
        self.processed = 0

    def update(self, x, force=False):
        """Update progress bar

        Args:
            x (int): The percentage of progress in [0, 100], if x equals input from the last time, the progress bar will not be updated
        """
        x = min(100, x)
        assert 0 <= x <= 100
        if self.last_x == int(x) and force is not True:
            return
        self.last_x = int(x)
        p = int(self.width * (x / 100.0))
        time_stamp = time.strftime("[%a %Y-%m-%d %H:%M:%S]", time.localtime())
        sys.stderr.write('\r%s [%-5s] [%s]' % (time_stamp, str(int(x)) + '%', '#' * p + '.' * (self.width - p)))
        sys.stderr.flush()

    def update_total(self, x):
        """Update total number of jobs

        Args:
            x (int): total number
        """
        self.total = int(x)

    def process(self, x):
        """Update progress bar using processed jobs

        Args:
            x (int): The number of jobs recently finished
        """
        self.processed += x
        self.update(100 * self.processed / self.total)

    def self_update(self):
        self.update(max(100, self.last_x + 1))

    def close(self):
        """Finish progress bar and print new line"""
        self.update(100)
        sys.stderr.write('\n')


class ProgManager(BaseManager):
    pass


ProgManager.register('Prog', ProgressBar)


def get_logger(logger_name='logger', fname=None, verbosity=False):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    level = logging.DEBUG if verbosity else logging.INFO

    # LOG format
    fmt = "%(asctime)-15s [%(levelname)-5s] %(message)s"
    datefmt = "[%a %Y-%m-%d %H:%M:%S]"
    formatter = logging.Formatter(fmt, datefmt)

    # LOG file
    if fname is not None:
        fh = open(fname, 'w')
        fh.close()
        file_handler = logging.FileHandler(fname)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # LOG console
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def find_logger_basefilename(logger):
    log_file = None
    parent = logger.__dict__['parent']
    if parent.__class__.__name__ == 'RootLogger':
        for h in logger.__dict__['handlers']:
            if h.__class__.__name__ == 'FileHandler':
                log_file = h.baseFilename
    else:
        log_file = find_logger_basefilename(parent)
    return log_file


def error_callback(e):
    """Error callback functions for multiprocessing Pool

    Args:
        e (Exception): Exception caught in processes
    """
    sys.stderr.write(f"Error in worker process: {e.__cause__}")