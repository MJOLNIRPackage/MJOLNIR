Convert Data File
=================

If a single file or a list of files needs to be converted, creating a complete python script might be a little too much and thus this script is provided. It converts the file(s) to the nxs format that can be loaded directly by MJOLNIR.

The script *MJOLNIRConvert* and it has the following help text::

    $ MJOLNIRConvert -h
    usage: MJOLNIRConvert [-h] [-s SAVE] [-b BINNING] [DataFile]

    Conversion tool for converting output h5 files to nxs files.

    positional arguments:
    DataFile              Data file to convert. If none provided file dialogue
                            will appear.

    optional arguments:
    -h, --help            show this help message and exit
    -s SAVE, --save SAVE  Location to which the generated file will be saved.
    -b BINNING, --binning BINNING
                            Binning performed. Default '8'

When no data files are provided a dialogue window will appear and the wanted files can be selected. This can be combined with the -s flag to specify a save location and/or the binning flag -b to specify the number 
of prismatic analyser pixels to be used per energy segment. 
