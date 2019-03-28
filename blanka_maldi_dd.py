import argparse, subprocess, os, copy, sys, pandas, numpy, timeit
import pyteomics.mzxml as pytmzxml
import pyteomics.mgf as pytmgf
from multiprocessing import Pool, cpu_count
from functools import partial

def mzxml_data_detection(directory):
    # scan directory for .mzXML files
    return [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(directory)
           for files in filenames if files.endswith('.mzXML')]

def raw_data_detection(args):
    # looks for data generated by Sanchez Lab Autoflex
    return [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(args['sample'])
            for files in filenames if files == 'fid']

def parse_maldi_template(args, msconvert_list):
    # parse MALDI DD names from IDBac Excel template
    template_df = pandas.read_excel(args['dd_template']).fillna(0)
    template_list = template_df.values.tolist()
    names_list = [(str(j), str(i[0]) + str(count)) for i in template_list for count, j in enumerate(i[1:], 1) if j != 0]
    return [(j, i[0] + '_' + i[1]) for i in names_list for j in msconvert_list if i[1] + '\\' in j]

def msconvert(args, msconvert_list):
    # convert raw data files detected into .mzXML format using MSConvert and default Sanchez Lab settings
    with open('config.ini', 'r') as config_file:
        msconvert_path = config_file.read().split('=')[1] + ' '
    for files in msconvert_list:
        if args['output'] == '':
            args['output'] = files[0][:files[0].find('fid')]
        msconvert = msconvert_path + files[0] + " -o " + args['output'] + " --outfile " + files[1] + \
                    ' --mzXML --32 --mz32 --inten32\
                    --filter "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>"\
                    --filter "peakPicking true 1-2"'
        print msconvert
        subprocess.call(msconvert)

def spectrum_to_dataframe(spectrum):
    # convert spectrum to dataframe
    mz_array = spectrum['m/z array'].newbyteorder('=')
    intensity_array = spectrum['intensity array'].newbyteorder('=')
    # fix byte ordering issues during multiprocessing
    # for linear processing: spectrum[''].byteswap().newbyteorder()
    return pandas.DataFrame({'m/z': mz_array, 'intensity': intensity_array}).astype({'m/z': numpy.float32,
                                                                                      'intensity': numpy.float32})
    # set all arrays as float32 to prevent different dtype error when using merge_asof

def load_control_data(args):
    # load control data files
    if args['output'] == '':
        directory = args['sample']
    else:
        directory = args['output']
    mzxml_list = [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(directory)
                  for files in filenames if files.startswith(args['control']) and files.endswith('.mzXML')]
    return [[[list(pytmzxml.read(mzxml))[0], mzxml], mzxml_list] for mzxml in mzxml_list]

def combine_control_spectra(noiseless_control_data):
    noiseless_control_dataframes = [pandas.DataFrame({'m/z': spectrum[0]['m/z array'].byteswap().newbyteorder(),
                                                 'intensity': spectrum[0]['intensity array'].byteswap().newbyteorder()})
                               .astype({'m/z': numpy.float32, 'intensity': numpy.float32})
                               for spectrum in noiseless_control_data]
    big_control_df = pandas.concat(noiseless_control_dataframes).sort_values(by=['m/z']).groupby(by=['m/z'],
                                                                                                as_index=False).sum()
    return {'m/z array': big_control_df['m/z'].values, 'intensity array': big_control_df['intensity'].values}

def noise_removal(signal_noise_ratio, spectrum):
    # remove noise using average of 5% lowest intensity peaks as noise level and user defined signal/noise ratio
    if len(spectrum['m/z array']) == 0:
        spectrum_dataframe = spectrum_to_dataframe(spectrum[0])
        spectrum_dataframe = spectrum_dataframe.drop(spectrum_dataframe[spectrum_dataframe['intensity'] == 0].index)
        if not spectrum_dataframe.empty:
            sorted_intensities = sorted(spectrum_dataframe['intensity'].values.tolist())
            lowest_sorted_intensities = sorted_intensities[:int(round((len(sorted_intensities)*0.05)))]
            if len(lowest_sorted_intensities) == 0:
                spectrum_noise = sum(lowest_sorted_intensities)/len(lowest_sorted_intensities)
                print spectrum_noise
                spectrum_dataframe = spectrum_dataframe.drop(spectrum_dataframe[spectrum_dataframe['intensity'] <= (signal_noise_ratio * spectrum_noise)].index)
                spectrum[0]['m/z array'] = spectrum_dataframe['m/z'].values
                spectrum[0]['intensity array'] = spectrum_dataframe['intensity'].values
    return [spectrum[0], spectrum[1]]

def blank_removal(peak_mz_tolerance, control_spectrum, sample_spectrum):
    # remove control spectrum peaks from sample spectrum if m/z within specified tolerance
    # returns processed sample spectrum and dictionary with removed peaks
    sample_spectrum_df = spectrum_to_dataframe(sample_spectrum[0])
    control_spectrum_df = spectrum_to_dataframe(control_spectrum)
    if not sample_spectrum_df.empty and not control_spectrum_df.empty:
        blank_subtracted_df = pandas.merge_asof(sample_spectrum_df, control_spectrum_df, on='m/z',
                                                tolerance=peak_mz_tolerance, direction='nearest').fillna(0)
        changed_peaks = blank_subtracted_df[blank_subtracted_df['intensity_y'] != 0.0]
        changed_dict = copy.deepcopy(sample_spectrum[0])
        changed_dict['m/z array'] = changed_peaks['m/z'].values
        changed_dict['intensity array'] = changed_peaks['intensity_x'].values
        blank_subtracted_df = blank_subtracted_df.drop(blank_subtracted_df[blank_subtracted_df['intensity_y']
                                                                           != 0.0].index)
        sample_spectrum[0]['m/z array'] = blank_subtracted_df['m/z'].values
        sample_spectrum[0]['intensity array'] = blank_subtracted_df['intensity_x'].values
        return [[sample_spectrum[0], sample_spectrum[1]], changed_dict]
    else:
        return [[sample_spectrum[0], sample_spectrum[1]], None]

def mgf_writer(spectrum_data_dict, output_dir, datatype):
    with open(output_dir + datatype + '_data_ms2.mgf', 'a') as mgf_file:
        if spectrum_data_dict['msLevel'] >= 2:
            mgf_file.write("BEGIN IONS" + "\n")
            mgf_file.write("TITLE=" + output_dir.split("\\")[-1].split("_blanka_")[0] + "." + spectrum_data_dict['num']
                       + "." + spectrum_data_dict['num'] + '. File:"' + output_dir.split("\\")[-1].split("_blanka_")[0]
                       + '.RAW", NativeID:"controllerType=0 controllerNumber=1 scan=' + spectrum_data_dict['num'] + '"'
                       + "\n")
            mgf_file.write("RTINSECONDS=" + str(spectrum_data_dict['retentionTime']) + "\n")
            mgf_file.write("PEPMASS=" + str(spectrum_data_dict['precursorMz'][0]['precursorMz']) + " "
                           + str(spectrum_data_dict['precursorMz'][0]['precursorIntensity']) + "\n")
            for mz, intensity in zip(spectrum_data_dict['m/z array'], spectrum_data_dict['intensity array']):
                mgf_file.write(str(mz) + ' ' + str(intensity) + "\n")
            mgf_file.write("END IONS" + "\n")
    with open(output_dir + datatype + '_data_full.mgf', 'a') as mgf_file:
        mgf_file.write("BEGIN IONS" + "\n")
        mgf_file.write("TITLE=" + output_dir.split("\\")[-1].split("_blanka_")[0] + "." + spectrum_data_dict['num']
                   + "." + spectrum_data_dict['num'] + '. File:"' + output_dir.split("\\")[-1].split("_blanka_")[0]
                   + '.RAW", NativeID:"controllerType=0 controllerNumber=1 scan=' + spectrum_data_dict['num'] + '"'                       + "\n")
        mgf_file.write("RTINSECONDS=" + str(spectrum_data_dict['retentionTime'] * 60) + "\n")
        if spectrum_data_dict['msLevel'] >= 2:
            mgf_file.write("PEPMASS=" + str(spectrum_data_dict['precursorMz'][0]['precursorMz']) + " "
                           + str(spectrum_data_dict['precursorMz'][0]['precursorIntensity']) + "\n")
        for mz, intensity in zip(spectrum_data_dict['m/z array'], spectrum_data_dict['intensity array']):
            mgf_file.write(str(mz) + ' ' + str(intensity) + "\n")
        mgf_file.write("END IONS" + "\n")

def old_mgf_writer(spectrum_data_dict, output_dir, datatype):
    # deprecated function
    # write spectrum data to .mgf file
    processed_sample_params = {key: spectrum_data_dict[key] for key in spectrum_data_dict.keys()
                               if key != 'm/z array' and key != 'intensity array'}
    processed_sample_spectra = [{'m/z array': spectrum_data_dict['m/z array'],
                                 'intensity array': spectrum_data_dict['intensity array'],
                                 'params': processed_sample_params}]
    pytmgf.write(spectra=processed_sample_spectra, output=output_dir + datatype + "_data.mgf")
