# based on https://realpython.com/pypi-publish-python-package
# How to upload:
#  - change package version in `setup.py` and `__init__.py`
#  - `python setup.py sdist`
#  - `twine upload dist/opengenomebrowser-tools-?.tar.gz`
import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / 'README.md').read_text()

setup(
    name='orthofinder-tools',
    version='0.0.2',
    description='Annotate orthogenes and create Roary-like plots',
    long_description=README,
    long_description_content_type='text/markdown',
    url='https://github.com/MrTomRod/orthofinder-tools/',
    author='Thomas Roder',
    author_email='roder.thomas@gmail.com',
    license='MIT',
    packages=['orthofinder_tools'],
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'fire',
        'matplotlib',
        'seaborn'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    entry_points={
        'console_scripts': [
            'annotate_orthogroups=orthofinder_tools.annotate_orthogroups:main',
            'orthofinder_plots=orthofinder_tools.orthofinder_plots:main',
        ],
    },
)
