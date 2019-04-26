from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.md')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()


pythonPath = 'lib/python{}.{}/site-packages/MJOLNIR/'.format(*sys.version_info[:2])

setup(
    name='MJOLNIR',

    version='1.1.0',
    description=('Neutron Scattering software suite.'),
    long_description=long_description,
    author='Jakob Lass',
    author_email='lass.jakob@gmail.com',
    url='https://github.com/jakob-lass/MJOLNIR',
    license='MPL-2.0',
    data_files = [(pythonPath, ["LICENSE.txt"])],#,(pythonPath+'/CommandLineScripts/',['MJOLNIR/CommandLineScripts/.settings'])],
    packages=['MJOLNIR','MJOLNIR/Data','MJOLNIR/Geometry','MJOLNIR/Statistics','MJOLNIR/CommandLineScripts'],
    scripts=['MJOLNIR/CommandLineScripts/MJOLNIRCalibrationInspector','MJOLNIR/CommandLineScripts/MJOLNIRHistory','MJOLNIR/CommandLineScripts/MJOLNIRConvert',
    'MJOLNIR/CommandLineScripts/MJOLNIR3DView'],
    python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*,<3.8.*',
    install_requires=['matplotlib>=2.0.2','numpy>=1.14','h5py>=2.5','scipy','datetime','shapely','pytest','pyperclip','shapely','decorator'],
    #data_files=[('MJOLNIR/TestData',['TestData/CAMEA_Full_2.xml','TestData/EnergyNormalization_2.calib','TestData/EnergyNormalization_3.calib','TestData/EnergyNormalization_8.calib','TestData/VanNormalization.h5','TestData/VanNormalization.nxs','TestData/cameasim2018n000005.nxs'])],
    
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Visualization',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    )
