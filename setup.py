from setuptools import setup
import setuptools, os, sys, platform

_here = os.path.abspath(os.path.dirname(__file__))
operatingSystem = sys.platform

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.md')) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()

if 'linux' in platform.platform().lower():
    dependency_links = ['https://download.pytorch.org/whl/cpu']
else:
    dependency_links = ['']

installFolder = os.path.abspath(os.path.join(os.path.split(setuptools.__file__)[0],'..'))
pythonPath =  os.path.relpath(os.path.join(installFolder,'MJOLNIR'),sys.base_prefix)

setup(
    name='MJOLNIR',
    version='1.3.2.post3',
    description=('Neutron Scattering software suite.'),
    long_description=long_description,
    author='Jakob Lass',
    author_email='lass.jakob@gmail.com',
    url='https://github.com/jakob-lass/MJOLNIR',
    license='MPL-2.0',
    #data_files = [(pythonPath, ["LICENSE.txt"]),((os.path.join(pythonPath),['MJOLNIR/Normalization_1.calib', 'MJOLNIR/Normalization_2.calib', 'MJOLNIR/Normalization_3.calib', 'MJOLNIR/Normalization_4.calib',
    #'MJOLNIR/Normalization_5.calib', 'MJOLNIR/Normalization_6.calib', 'MJOLNIR/Normalization_7.calib', 'MJOLNIR/Normalization_8.calib'])),((os.path.join(pythonPath),['MJOLNIR/CalibrationFlatCone.csv'])),((os.path.join(pythonPath),['MJOLNIR/CalibrationMultiFLEXX.csv'])),((os.path.join(pythonPath),['MJOLNIR/CalibrationBambus.csv'])),
    #            ((os.path.join(pythonPath,'Geometry'),['MJOLNIR/Geometry/detsequence.dat']))],#,(pythonPath+'/CommandLineScripts/',['MJOLNIR/CommandLineScripts/.settings'])],
    include_package_data=True,
    package_data ={'license':['license.txt'], 'Normalization_1':['Normalization_1.calib'], 'Normalization_2':['Normalization_2.calib'], 'Normalization_3':['Normalization_3.calib'],
                   'Normalization_4':['Normalization_4.calib'], 'Normalization_5':['Normalization_5.calib'], 'Normalization_6':['Normalization_6.calib'], 'Normalization_7':['Normalization_7.calib'], 
                    'Normalization_8':['Normalization_8.calib'],  'CalibrationFlatCone':['CalibrationFlatCone.csv'], 'CalibrationMultiFLEXX':['CalibrationMultiFLEXX.csv'],
                    'CalibrationBambus':['CalibrationBambus.csv'], 'Geometry/detsequence.dat':['Geometry/detsequence.dat'], 'CommandLineScripts/.settings':['CommandLineScripts/.settings']},

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
    install_requires=['matplotlib>=3','numpy>=1.14,<2.0','h5py>=2.5','scipy','datetime','pytest>=4.6','pyperclip','decorator','pandas','future','sympy',
                    'pip>=20','ufit>=1.4.0','pyqtgraph','regex','torch','torchvision','torchaudio'], # ,'ufit','sip','PyQt5-sip','PyQt5<=5.12'
    dependency_links = dependency_links,
    extras_require={
        "gui": [
            "MJOLNIRGui",
        ]},
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
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11'],
    )
