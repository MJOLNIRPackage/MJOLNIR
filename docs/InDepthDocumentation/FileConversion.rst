.. _`Data file conversion`:

Data file conversion
====================

In general, it is expected that a CAMEA-like instrument is to be run during experiments in 2 different scan modes:

    - :math:`A3` scans

    - External parameter

Most often, one would rotate the back-end into a suitable :math:`A4` and :math:`E_i` position to cover the interesting physics, and then perform an :math:`A3` scan. The length of this scan depends on the symmetry of the crystal in the given scattering plane as performing a 360 degrees scan with a 90 degrees symmetry does not provide additional information. After such a scan, would rotate :math:`A4` by half a wedge coverage angle (3.75 degrees) to cover the dark angles and then perform an identical :math:`A3` scan. This could be performed with different incoming energies to expand the covered area in the energy direction.

.. Having the raw data in the H5 format, converting the data files into :math:`S(\vec{q},\omega)` is rather strraihgt forward. 



NXsqom file format
------------------

The format into which data is converted is the `NXsqom <http://download.nexusformat.org/sphinx/classes/applications/NXsqom.html>`_ format. It is a standard of the nexus files and is designed for data converted into reciprocal space. With this choice of conversion it is believed that some pre-existing data handling routines exist in other software solutions already. 


Below is a HDF converted file in the NXsqom format for a :math:`A3` scan. Here :math:`NP` is the number of scan points and :math:`NNP` is the number of unique pixels converted.

::

    cameasim2018n0000xx.nxs
    └── entry
        ├── CAMEA (NX_class = NXinstrument)
        │   ├── detector (NX_class = NXdetector)
        │   │   ├── data (32-bit integer,    NP x 104 x 452)
        │   │   └── polar_angle (32-bit floating-point,    1)
        │   │
        │   └── monochromator
        │       ├── d_spacing (32-bit floating-point,    1)
        │       ├── energy (32-bit floating-point,    1)
        │       ├── rotation_angle (32-bit floating-point,    1)
        │       └── type (String, length = 70,    1)
        │
        ├── control (NX_class = NXmonitor)
        │   ├── data (32-bit integer,    NP)
        │   ├── mode (String, length = 70,    1)
        │   ├── preset (32-bit integer, 1)
        │   └── time (32-bit floating-point,    NP)
        │
        ├── data (NX_class = NXdata)
        │   ├── data (32-bit integer,    NNP)
        │   ├── en (32-bit floating-point,    NNP)
        │   ├── monitor (32-bit integer,    NNP)
        │   ├── normalization (32-bit floating-point,    NNP)
        │   ├── qx (32-bit floating-point,    NNP)
        │   ├── qy (32-bit floating-point,    NNP)
        │   ├── qz (32-bit floating-point,    NNP)
        │   ├── rawdata (32-bit integer,    NP x 104 x 452)
        │   └── rotation_angle (32-bit floating-point,    NP)
        │
        ├── definition (String, length = 70,    1)
        │
        ├── end_time (String, length = 70,    1)
        │
        ├── proposal_id (String, length = 70,    1)
        │
        ├── proposal_user (NX_class = NXuser)
        │   └── name (String, length = 70,    1)
        │
        ├── reduction (NX_class = NXprocess)
        │   └── MJOLNIR_algorithm_convert (NX_class = NXprocess)
        │       ├── author (String, length = 70,    1)
        │       ├── date (String, length = 70,    1)
        │       ├── description (String, length = 70,    1)
        │       ├── normalization table (String, length = 70,    1)
        │       └── rawdata (String, length = 70,    1)
        │
        ├── sample (NX_class = NXsample)
        │   ├── name (String, length = 70,    1)
        │   ├── orientation_matrix (32-bit floating-point,    3 x 3)
        │   ├── plane_normal (32-bit floating-point,    3)
        │   ├── polar_angle (32-bit floating-point,    1)
        │   ├── rotation_angle (32-bit floating-point,    NP)
        │   └── unit_cell (32-bit floating-point,    6)
        │
        ├── start_time (String, length = 70,    1)
        │
        └── title (String, length = 70,    1)

