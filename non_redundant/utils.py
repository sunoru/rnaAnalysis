#!/usr/bin/python2

import os
import datetime
import json
import requests
import re
import logging
from logging import info, warn, error

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
logging.basicConfig(filename="logs/%s.log" % datetime.date.today().strftime("%Y%m%d"), level=logging.INFO)

def fetch_raw(url):
    x = requests.get(url)
    return x.content

def findall(pattern, string, flags=0):
    return re.findall(pattern, string, flags)

def data_file(filename):
    return os.path.join(DATA_DIR, filename)
