import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cctyper", 
    version="1.8.0",
    author="Jakob Russel",
    author_email="russel2620@gmail.com",
    description="CRISPRCasTyper: Automatic detection and subtyping of CRISPR-Cas operons",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Russel88/CRISPRCasTyper",
    download_url="https://github.com/Russel88/CRISPRCasTyper/archive/v1.8.0.tar.gz",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 4 - Beta"],
    python_requires='>=3.10',
    install_requires=[
        "numpy >= 1.26",
        "pandas >= 2.2",
        "scipy >= 1.12",
        "biopython >= 1.85",
        "multiprocess >= 0.70.18",
        "scikit-learn >= 1.4",
        "xgboost >= 2.0",
        "tqdm >= 4.66",
        "drawsvg >= 2.4.1",
        "setuptools"],
    scripts=['bin/cctyper',
             'bin/repeatType',
             'bin/repeatTrain']
)
