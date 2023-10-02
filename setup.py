import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tfunct",
    version="1.0.1.dev",
    author="Ernesto Giron Echeverry",
    author_email="e.giron.e@gmail.com",
    description="Library for estimating grain yield using temperature response functions",
    keywords="PRFT, WEFT, TPF, VPD, wheat, iPAR, photoperiod, RUE, NDVI, crop modeling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/egiron/TemperatureFunct",
    packages=setuptools.find_packages(),
    include_package_data=True,
    #packages=['tfunct'],
    # scripts=['bin/tfunct'],
    install_requires=[
        'numpy >=1.22.4',
        'numba >=0.51.2',
        'pandas >=1.5.3',
        'scikit-learn >=1.2.2',
        'scipy >=1.10.1',
        'tqdm >=4.64.0',
        'seaborn >=0.11.0',
        'Shapely >=1.7.1',
        'ipython >=7.21.0',
        'duckdb >=0.8.1',
        'pyarrow',
        'click'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Free for non-commercial use",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "License :: Free for non-commercial use",
        "Environment :: Console",
        "Framework :: Jupyter",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research"
    ],
)