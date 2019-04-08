"""
Setup module for standard Pypi installation.
"""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biobb_structure_checking",
    version="0.0.1",
    author="Biobb developers",
    author_email="josep.gelpi@bsc.es",
    description="BioBB_structure_checking performs MDWeb structure checking set as a command line utility.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="Bioinformatics Workflows BioExcel Compatibility",
    url="https://github.com/bioexcel/biobb_analysis",
    project_urls={
        "Documentation": "http://BioBB_structure_checking.readthedocs.io/en/latest/",
        "Bioexcel": "https://bioexcel.eu/"
    },
    packages=setuptools.find_packages(exclude=['docs', 'test']),
    install_requires=['biobb_structure_manager==0.0.6'],
    python_requires='==3.6.*',
    classifiers=(
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ),
)
