import os
import sys
import platform

from setuptools import setup, find_packages

if not (sys.platform == "win32" and platform.architecture()[0] == "64bit"):
    raise Exception("Built only for Windows x64 for now...")

os.chdir(os.path.dirname(sys.argv[0]) or ".")

setup(
    name="liquidcffi",
    version="0.1",
    description="Python bindings for Liquid-DSP using CFFI",
    long_description=open("README.md", "rt").read(),
    url="https://github.com/ruilvo",
    author="Rui Oliveira",
    author_email="ruilvo@ua.pt",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
        "License :: OSI Approved :: MIT License",
    ],
    packages=find_packages(),
    install_requires=["cffi>=1.0.0"],
    setup_requires=["cffi>=1.0.0"],
    cffi_modules=[
        "./liquidcffi/build_liquidcffi.py:ffibuilder",
    ],
    data_files=[('', ["liquidcffi/lib/win64/libliquid.dll"])],
)
