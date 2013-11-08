# -*- coding:utf-8 -*-

import gzip

class gzopen(object):
   def __init__(self, fname):
      f = open(fname)
      magic_number = f.read(2)
      f.seek(0)
      if magic_number == '\x1f\x8b':
         self.f = gzip.GzipFile(fileobj=f)
      else:
         self.f = f
   def __exit__(self, type, value, traceback):
      try:
         self.f.fileobj.close()
      except AttributeError:
         pass
      finally:
         self.f.close()
   def __enter__(self):
      return self
   def __getattr__(self, name):
      return getattr(self.f, name)
   def __iter__(self):
      return iter(self.f)
