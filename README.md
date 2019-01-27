# BLANKA

BLANKA is a command line tool used to remove noise and blank background signals from mass spectrometry data.


## Workflow
![Blanka_Workflow](https://user-images.githubusercontent.com/46392631/51808383-c63f1480-2258-11e9-98d3-57b759b97d76.png)

## Installation
No installation is required for use of BLANKA. Simply download the scripts and place them in the desired directory. The config.ini file must be editted to include the path to a local installation of MSConvert before usage.

## Dependencies
Python 2.7.15\
argparse 1.3.0 or higher\
pyteomics 3.5.1 or higher\
numpy 1.15.4 or higher\
pandas. 0.23.4 or higher (must be released Nov 2018 or later - at the time of writing, Anaconda has the Aug 2018 release of pandas, and a newer release must be manually installed from Github or another repository)\
subprocess\
os\
sys\
multiprocessing\
functools\
copy

## Parameters
#### Required
-s [--sample] : sample input directory with single or multiple files or sample .mzXML file\
-c [--control] : LCQ/QTOF mode = control input file path with .mzXML file extension, 
                 DD mode = control sample spot name\
-i [--instrument] : instrument/protocol used for experiment ('lcq', 'qtof', 'dd', or 'ims')\
#### Optional
--dd_template : dried droplet excel sheet with sample names (same template as IDBac) (required if instrument = 'dd')\
-o [--output] : output directory for all generated files (default = sample directory)\
--cpu : number of threads used (default = max-1)\
--signal_noise_ratio : integer signal to noise ratio used for noise removal (default = 4)\
-r [--retention_time_tolerance] : retention time tolerance range in seconds (default = 0.1 sec)\
-m [--precursor_mz_tolerance] : absolute precursor m/z error tolerance in Da (default = 0.02 Da)\
-p [--peak_mz_tolerance] : absolute precursor m/z error tolerance in Da (default = 0.02 Da)\
--noise_removal_only : only perform noise removal (default = False)\
--blank_removal_only : only perform blank removal (default = False)

## Examples
Print usage information.\
```python blanka```

Perform noise and blank removal on data actinomycetes.mzXML\
```python blanka -s E:\lcms_data\actinomycetes.mzXML -c E:\lcms_data\media_control.mzXML -i lcq```

Perform noise and blank removal on data found in E:\lcms_data and outputs data to E:\blanka_output\
```python blanka -s E:\lcms_data -c E:\lcms_data\media_control.mzXML -o E:\blanka_output -i lcq```

Perform noise and blank removal on data found in E:\lcms_data\sample using multiple control files found in E:\lcms_data\control\
```python blanka -s E:\lcms_data\sample -c E:\lcms_data\control -o E:\blanka_output -i lcq```

Convert .RAW data to .mzXML using MSConvert and perform noise and blank removal on actinomycetes.mzXML\
```python blanka -s E:\lcms_data\actinomycetes.RAW -c E:\lcms_data\media_control.RAW -i lcq```

Perform noise and blank removal on data actinomycetes.mzXML with custom retention time and precursor mz tolerance\
```python blanka -s E:\lcms_data\actinomycetes.mzXML -c E:\lcms_data\media_control.mzXML -i lcq -r 0.5 -m 0.1```

Performs noise and blank removal on dried droplet maldi data using 'media_control' spots as control\
```python blanka -s E:\maldi_data\ -c media_control -o E:\blanka_output -i dd```

Performs blank removal only on dried droplet maldi data using 'media_control' spots as control\
```python blanka -s E:\maldi_data\ -c media_control -o E:\blanka_output -i dd --blank_removal_only True```
