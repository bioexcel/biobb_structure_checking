"""
Setup module for standard Pypi installation.
"""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biobb_structure_checking",
    version="3.15.6",

    author="Biobb developers",
    author_email="gelpi@ub.edu",
    description="BioBB_structure_checking performs MDWeb structure checking tool set as a command line utility.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="Bioinformatics Workflows BioExcel Compatibility",
    url="https://github.com/bioexcel/biobb_structure_checking",
    project_urls={
        "Documentation": "http://BioBB_structure_checking.readthedocs.io/en/latest/",
        "Bioexcel": "https://bioexcel.eu/"
    },
    packages=setuptools.find_packages(exclude=['docs', 'test']),
    package_data={'biobb_structure_checking': ['py.typed']},
    include_package_data=True,
    install_requires=['psutil', 'numpy', 'biopython>=1.79'],
    python_requires='>=3.8',
    extras_require={},
    entry_points={
        "console_scripts": [
            "check_structure = biobb_structure_checking.check_structure:main"
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
    ],
)
