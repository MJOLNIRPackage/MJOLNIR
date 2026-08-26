MJOLNIR
============
Neutron scattering software to be used at CAMEA-like Neutron Instruments.

## Installation
This package is to be install through the Python Package Manager by issuing one of the two commands

```shell
pip install MJOLNIR[qt5]
```

or 

```shell
pip install MJOLNIR[qt6]
```

The reason is that MJOLNIR supports voth pyqt5 and pyqt6 but require one of them. This is done for python version compatibility. If MJOLNIR is installed without the qt version specified the package will inform the user upon first import of the package/modules requiring pyqt.

For full documentation see http://mjolnir.readthedocs.io/en/latest/ or https://www.psi.ch/en/sinq/camea/data-treatment. 

