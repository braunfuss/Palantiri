#!/usr/bin/env python

from setuptools import setup
from setuptools.command.install import install


class CustomInstallCommand(install):
    def run(self):
        install.run(self)


setup(
    cmdclass={
        'install': CustomInstallCommand,
    },

    name='palantiri',

    description='A python based seismological backprojection tool.\
                 In Quenya, palantíri means "far-seeing"',

    version='0.4',

    author='Andreas Steinberg',

    author_email='andreas.steinberg@ifg.uni-kiel.de',

    packages=[
        'palantiri',
        'palantiri.apps',
        'palantiri.cluster',
        'palantiri.common',
        'palantiri.process',
        'palantiri.skeleton',
        'palantiri.tools',
        'palantiri.data',
        'palantiri.waveform',
    ],
    python_requires='>=3.5',
    entry_points={
        'console_scripts': [
            'bat = palantiri.apps.arraytool:main',
            'palantiri_plot = palantiri.apps.plot:main',
            'palantiri_down = palantiri.apps.palantiri_down:main',
            'palantiri_eventsearch = palantiri.tools.eventsearch:main',
            'palantiri_cluster = palantiri.cluster.callcluster:main',
            'palantiri_create = palantiri.tools.create:main',
            'palantiri_process = palantiri.process.main:main',

        ]
    },
    package_dir={'palantiri': 'src'},

    package_data={
        'palantiri': [
            'data/*nd',
            'skeleton/example.config']},
    data_files=[],

    license='GPLv3',

    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Application Frameworks',
        ],

    keywords=[
        'seismology, waveform analysis, earthquake modelling, geophysics,'
        ' backprojection'],
    )
