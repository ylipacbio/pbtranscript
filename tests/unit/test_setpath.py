#!/usr/bin/python
from os import path
import ConfigParser

"""Define test data path for pbtranscript."""

THIS_DIR = path.dirname(path.abspath(__file__))
ROOT_DIR = path.dirname(THIS_DIR)
NOSE_CFG = path.join(THIS_DIR, "nose.cfg")

def _get_siv_data_std_dir():
    """Get siv data/stdout directories for pbtranscript unit tests
    """
    nosecfg = ConfigParser.SafeConfigParser()
    nosecfg.readfp(open(NOSE_CFG), 'r')
    if nosecfg.has_section('data'):
        data_dir = path.abspath(nosecfg.get('data', 'SIV_DATA_DIR'))
        std_dir = path.abspath(nosecfg.get('data', 'SIV_STD_DIR'))
        return data_dir, std_dir
    else:
        msg = "Unable to find section [DATA] option [dataDir]" + \
              "and [stdDir] in config file {f}.".format(f=NOSE_CFG)
        raise KeyError(msg)

OUT_DIR  = path.join(ROOT_DIR, "out")
DATA_DIR = path.join(ROOT_DIR, "data")
STD_DIR  = path.join(ROOT_DIR, "stdout")
SIV_DATA_DIR, SIV_STD_DIR = _get_siv_data_std_dir()

