#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = [ ]

setup(
    author="Edvin Fako",
    author_email='edvin.fako@basf.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A heuristic approach of generating molecules and intermediates for at heterogeneous interfaces.",
    entry_points={
        'console_scripts': [
            'autoadsorbate=autoadsorbate.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='autoadsorbate',
    name='autoadsorbate',
    packages=find_packages(include=['autoadsorbate', 'autoadsorbate.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/fakoed/autoadsorbate',
    version='0.2.0',
    zip_safe=False,
)
