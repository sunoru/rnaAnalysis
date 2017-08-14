#!/usr/bin/python2

import os
import datetime
import json
import requests
import re
import logging
from logging import info, warn, error

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
LOGS_DIR = os.path.join(BASE_DIR, "logs")
DATA_DIR = os.path.join(BASE_DIR, "data")
TEST_DIR = os.path.join(BASE_DIR, "testfiles")

logging.basicConfig(filename=os.path.join(LOGS_DIR, "%s.log" % datetime.date.today().strftime("%Y%m%d")), level=logging.INFO)

def fetch_raw(url):
    x = requests.get(url)
    return x.content if x.ok else None

def findall(pattern, string, flags=0):
    return re.findall(pattern, string, flags)

def data_file(filename):
    return os.path.join(DATA_DIR, filename)
