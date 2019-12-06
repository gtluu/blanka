# BLANKA

BLANKA is a command line tool used to remove noise and blank background signals from mass spectrometry data.

This repo contains the publication version of BLANKA. For the current version please see https://github.com/gtluu/blanka_v1.


## Workflow
![Blanka_Workflow](https://user-images.githubusercontent.com/46392631/51808383-c63f1480-2258-11e9-98d3-57b759b97d76.png)

## Installation
No installation is required for use of BLANKA. Simply download the scripts and place them in the desired directory. The config.ini file must be editted to include the path to a local installation of MSConvert before usage.

## Dependencies
Python 2.7.15\
argparse 1.3.0 or higher\
pyteomics 3.5.1 or higher\
numpy 1.15.4 or higher\
pandas 0.24.0 ~~0.23.4~~ or higher (the bug fix required is now a part of pandas 0.24.0) ~~(must be released Nov 2018 or later - at the time of writing, Anaconda has the Aug 2018 release of pandas, and a newer release must be manually installed from Github or another repository)~~\
subprocess\
os\
sys\
multiprocessing\
functools\
copy

## Parameters
#### Required
--sample : sample input directory with single or multiple files or sample .mzXML file\
--control : LCQ/QTOF mode = control input file path with .mzXML file extension, 
            DD mode = control sample spot name\
--instrument : instrument/protocol used for experiment ('lcq', 'qtof', 'dd', or 'ims')\
#### Optional
--dd_template : dried droplet excel sheet with sample names (same template as IDBac) (required if instrument = 'dd')\
--output : output directory for all generated files (default = sample directory)\
--cpu : number of threads used (default = max-1)\
--signal_noise_ratio : integer signal to noise ratio used for noise removal (default = 4)\
--retention_time_tolerance : retention time tolerance range in seconds (default = 0.1 sec)\
--precursor_mz_tolerance : absolute precursor m/z error tolerance in Da (default = 0.02 Da)\
--peak_mz_tolerance : absolute precursor m/z error tolerance in Da (default = 0.02 Da)\
--noise_removal_only : only perform noise removal (default = False)\
--blank_removal_only : only perform blank removal (default = False)

## Examples
Print usage information.\
```python blanka```

Perform noise and blank removal on data actinomycetes.mzXML\
```python blanka --sample E:\lcms_data\actinomycetes.mzXML --control E:\lcms_data\media_control.mzXML --instrument lcq```

Perform noise and blank removal on data found in E:\lcms_data and outputs data to E:\blanka_output\
```python blanka --sample E:\lcms_data --control E:\lcms_data\media_control.mzXML --output E:\blanka_output --instrument lcq```

Perform noise and blank removal on data found in E:\lcms_data\sample using multiple control files found in E:\lcms_data\control\
```python blanka --sample E:\lcms_data\sample --control E:\lcms_data\control --output E:\blanka_output --instrument lcq```

Convert .RAW data to .mzXML using MSConvert and perform noise and blank removal on actinomycetes.mzXML\
```python blanka --sample E:\lcms_data\actinomycetes.RAW --control E:\lcms_data\media_control.RAW --instrument lcq```

Perform noise and blank removal on data actinomycetes.mzXML with custom retention time and precursor mz tolerance\
```python blanka --sample E:\lcms_data\actinomycetes.mzXML --control E:\lcms_data\media_control.mzXML --instrument lcq --retention_time_tolerance 0.5 --peak_mz_tolerance 0.1```

Performs noise and blank removal on dried droplet maldi data using 'media_control' spots as control\
```python blanka --sample E:\maldi_data\ --control media_control --output E:\blanka_output --instrument dd```

Performs blank removal only on dried droplet maldi data using 'media_control' spots as control\
```python blanka --sample E:\maldi_data\ --control media_control --output E:\blanka_output --instrument dd --blank_removal_only True```
