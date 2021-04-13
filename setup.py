from setuptools import setup
import setuptools
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))
operatingSystem = sys.platform

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.md')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()


installFolder = os.path.abspath(os.path.join(os.path.split(setuptools.__file__)[0],'..'))
pythonPath =  os.path.relpath(os.path.join(installFolder,'MJOLNIR'),sys.base_prefix)

setup(
    name='MJOLNIR',

    version='1.1.18',
    description=('Neutron Scattering software suite.'),
    long_description=long_description,
    author='Jakob Lass',
    author_email='lass.jakob@gmail.com',
    url='https://github.com/jakob-lass/MJOLNIR',
    license='MPL-2.0',
    data_files = [(pythonPath, ["LICENSE.txt"]),((os.path.join(pythonPath),['MJOLNIR/Normalization_1.calib'])),((os.path.join(pythonPath),['MJOLNIR/CalibrationFlatCone.csv'])),((os.path.join(pythonPath),['MJOLNIR/CalibrationMultiFLEXX.csv'])),
                ((os.path.join(pythonPath,'Geometry'),['MJOLNIR/Geometry/detsequence.dat']))],#,(pythonPath+'/CommandLineScripts/',['MJOLNIR/CommandLineScripts/.settings'])],
    packages=['MJOLNIR','MJOLNIR/Data','MJOLNIR/Geometry','MJOLNIR/Statistics','MJOLNIR/CommandLineScripts'],
    #scripts=['MJOLNIR/CommandLineScripts/MJOLNIRCalibrationInspector','MJOLNIR/CommandLineScripts/MJOLNIRHistory','MJOLNIR/CommandLineScripts/MJOLNIRConvert',
    #'MJOLNIR/CommandLineScripts/MJOLNIR3DView'],
    entry_points = {
        "console_scripts": ['MJOLNIRHistory = MJOLNIR.CommandLineScripts.MJOLNIRHistory:main',
                            'MJOLNIRConvert = MJOLNIR.CommandLineScripts.MJOLNIRConvert:main',
                            'MJOLNIRCalibrationInspector = MJOLNIR.CommandLineScripts.MJOLNIRCalibrationInspector:main',
                            'MJOLNIR3DView = MJOLNIR.CommandLineScripts.MJOLNIR3DView:main']
        },
    python_requires='>=3.5' if not operatingSystem == 'darwin' else '>=3.6',
    install_requires=['matplotlib>=3,<3.4','numpy>=1.14','h5py>=2.5','scipy','datetime','pytest>=4.6','pyperclip','decorator','pandas','future',
                    'pip>=20','ufit>=1.4.0','pyqtgraph','regex'], # ,'ufit','sip','PyQt5-sip','PyQt5<=5.12'
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'],
    )
