"""
CFFI builder for Liquid-DSP
"""

import sys, platform

if not (sys.platform == "win32" and platform.architecture()[0] == "64bit"):
    raise Exception("Built only for Windows x64 for now...")

import os
from cffi import FFI

curr_dir = os.path.abspath(os.path.dirname(__file__))

ffibuilder = FFI()

ffibuilder.set_source(
    "liquidcffi",
    r"""
        #include "liquid.pre.nc.h"
    """,
    libraries=["libliquid"],
    library_dirs=[os.path.normpath(os.path.join(curr_dir, "lib/win64/"))],
    include_dirs=[os.path.normpath(os.path.join(curr_dir, "include/"))],
    #extra_objects=[]
)

ffibuilder.cdef(open(os.path.normpath(os.path.join(curr_dir, "include/liquid.pre.nc.h"))).read())

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
