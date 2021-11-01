from setuptools import setup, find_packages

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name = 'charcoal-bio',
    version = "0.1",
    description="a tool for decontaminating genomes",
    url="https://github.com/dib-lab/charcoal",
    author="C. Titus Brown and Taylor Reiter",
    author_email="titus@idyll.org,tereiter@ucdavis.edu",
    license="BSD 3-clause",
    packages = find_packages(),
    classifiers = CLASSIFIERS,
    entry_points = {'console_scripts': [
        'charcoal  = charcoal.__main__:main'
        ]
    },
    include_package_data=True,
    package_data = { "charcoal": ["Snakefile", "*.yml", "*.ipynb"] },
    setup_requires = [ "setuptools>=38.6.0",
                       'setuptools_scm', 'setuptools_scm_git_archive' ],
    use_scm_version = {"write_to": "charcoal/version.py"},
    install_requires = ['snakemake==6.10.0', 'click>=7,<8']
)
