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
            ├── CAMEA
            │    ├── detector
            │    │    ├── data
            │    │    ├── polar_angle
            │    └── monochromator
            │         ├── d_spacing
            │         ├── energy
            │         ├── rotation_angle
            │         └── type
            ├── calibration
            │    ├── 1_pixels
            │    ├── 3_pixels
            │    ├── 8_pixels
            ├── control
            │    ├── data
            │    ├── mode
            │    ├── preset
            │    ├── time
            ├── data
            │    ├── data
            │    ├── en
            │    ├── monitor
            │    ├── normalization
            │    ├── qx
            │    ├── qy
            │    ├── qz
            │    └── rotation_angle
            ├── definition
            ├── end_time
            ├── proposal_id
            ├── proposal_user
            │    └── name
            ├── reduction
            │    └── MJOLNIR_algorithm_convert
            │         ├── author
            │         ├── binning
            │         ├── date
            │         ├── description
            │         └── rawdata
            ├── sample
            │    ├── name
            │    ├── orientation_matrix
            │    ├── plane_normal
            │    ├── polar_angle
            │    └── unit_cell
            ├── start_time
            └── title

