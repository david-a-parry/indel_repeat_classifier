try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import re
v_file = "indel_repeat_classifier/version.py"
v_line = open(v_file, "rt").read()
v_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
match = re.search(v_re, v_line, re.M)
if match:
    verstr = match.group(1)
else:
    raise RuntimeError("Unable to find version string in {}.".format(v_file))

setup(
    name="indel_repeat_classifier",
    packages=["indel_repeat_classifier"],
    scripts=[
        "bin/short_repeats_from_fasta", "bin/filter_repeat_csv_by_regions"
    ],
    version=verstr,
    description="Module for classifying repeat context of indels from VCFs.",
    author="David A. Parry",
    author_email="david.parry@ed.ac.uk",
    url='https://git.ecdf.ed.ac.uk/RER_deletions_paper/indel_repeat_classifier',
    download_url=
    'https://git.ecdf.ed.ac.uk/RER_deletions_paper/indel_repeat_classifier/archive/{}.tar.gz'
    .format(verstr),
    license='MIT',
    install_requires=['pyfaidx'],
    tests_require=['nose2', 'pysam'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
