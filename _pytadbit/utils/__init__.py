import datetime
from time           import time


def printime(msg):
    print(msg +
          (' ' * (79 - len(msg.replace('\n', '')))) +
          '[' +
          str(datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
          ']')
