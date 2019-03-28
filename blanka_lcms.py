import argparse, subprocess, os, copy, sys, pandas, numpy, timeit
import pyteomics.mzxml as pytmzxml
import pyteomics.mgf as pytmgf
from multiprocessing import Pool, cpu_count
from functools import partial

def mzxml_data_detection(directory):
    # scan directory for .mzXML files
    return [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(directory)
           for files in filenames if files.endswith('.mzXML')]

def raw_data_detection(args, file_directory):
    # scan directory for raw data
    if args['instrument'] == 'lcq':
        # looks for data generated from Sanchez Lab Thermo LCQ
        directory_list = [os.path.join(dirpath, directory) for dirpath, dirnames, filenames in os.walk(file_directory)
                          for directory in dirnames]
        if directory_list == []:
            directory_list = [file_directory]
        file_list = [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(file_directory)
                     for files in filenames if files.endswith('.RAW')]
        return [(files, directory) for files in file_list for directory in directory_list]
    elif args['instrument'] == 'qtof':
        # looks for data generated from COP QTOF
        file_list = [os.path.join(dirpath, directory) for dirpath, dirnames, filenames in os.walk(file_directory)
                     for directory in dirnames if directory.endswith('.d')]
        #msconvert_list = [(files, files.split('\\')[-1]) for files in file_list]
        msconvert_list = [(files, files[:files.rfind('\\')]) for files in file_list]
        return msconvert_list

def msconvert(args, msconvert_list):
    # convert raw data files detected into .mzXML format using MSConvert and default Sanchez Lab settings
    with open('config.ini', 'r') as config_file:
        msconvert_path = config_file.read().split('=')[1] + ' '
    for files in msconvert_list:
        if args['output'] == '':
            args['output'] = files[1]
        msconvert = msconvert_path + files[0] + " -o " + args['output'] + ' --mzXML --32 --mz32 --inten32\
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
    return pandas.DataFrame({'m/z': mz_array, 'intensity': intensity_array})

def load_control_data(args):
    # load control dataset
    if args['control'] == '':
        print "No control file specified. Exiting BLANKA."
        sys.exit(1)
    if args['control'].endswith('.mzXML'):
        return list(pytmzxml.read(args['control']))
    elif args['control'].endswith('d') or args['control'].endswith('.RAW'):
        args['control'] = args['control'].split('.')[0] + '.mzXML'
        return list(pytmzxml.read(args['control']))
    else:
        control_files = [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(args['control'])
                         for files in filenames if files.endswith('.mzXML')]
        if control_files == []:
            control_raw_list = raw_data_detection(args, args['control'])
            msconvert(args, control_raw_list)
            control_files = [os.path.join(dirpath, files) for dirpath, dirnames, filenames in os.walk(args['control'])
                             for files in filenames if files.endswith('.mzXML')]
        control_data = []
        for control in control_files:
            control_data += list(pytmzxml.read(control))
        print len(control_data)
        return control_data

def select_control_spectrum(args, ms_mode, ret_time, precursor_mz, control_spectra_data):
    # select spectrum to use as control from control dataset based on given parameters
    if ms_mode == 1:
        # search by retention time for MS1
        control_spectra_list = [i for i in control_spectra_data if i['msLevel'] == 1
                                if ((i['retentionTime'] * 60) + args['retention_time_tolerance']) >= ret_time >=
                                ((i['retentionTime'] * 60) - args['retention_time_tolerance'])]
    elif ms_mode >= 2:
        # search by retention time and precursor mz for MSn
        control_spectra_list = [i for i in control_spectra_data if i['msLevel'] == 2
                                if ((i['retentionTime'] * 60) + args['retention_time_tolerance']) >= ret_time >=
                                ((i['retentionTime'] * 60) - args['retention_time_tolerance']) and
                                (i['precursorMz'][0]['precursorMz'] + args['precursor_mz_tolerance']) >=
                                precursor_mz >= (i['precursorMz'][0]['precursorMz'] -
                                                 args['precursor_mz_tolerance'])]
    if len(control_spectra_list) != 0:
    # return None if no matching spectra found
        if len(control_spectra_list) > 1:
            # return best match if more than one spectra found meeting criteria
            if ms_mode == 1:
                ret_time_indexes = [abs(ret_time - (i['retentionTime'] * 60)) for i in control_spectra_list]
                return control_spectra_list[ret_time_indexes.index(min(ret_time_indexes))]
            elif ms_mode >= 2:
                # scores each spectrum based on how close retention time and precursor mz are compared to each spectra
                ret_time_indexes = [abs(ret_time - (i['retentionTime'] * 60)) for i in control_spectra_list]
                precursor_mz_indexes = [abs(precursor_mz - i['precursorMz'][0]['precursorMz'])
                                        for i in control_spectra_list]
                if len(ret_time_indexes) == len(precursor_mz_indexes):
                    score = [0] * (len(ret_time_indexes))
                    rettime_sorted_indexes = numpy.argsort(ret_time_indexes)
                    precursor_sorted_indexes = numpy.argsort(precursor_mz_indexes)
                    for count, i in enumerate(zip(rettime_sorted_indexes, precursor_sorted_indexes), 1):
                        score[i[0]] += count
                        score[i[1]] += count
                    return control_spectra_list[score.index(max(score))]
        elif len(control_spectra_list) == 1:
            return control_spectra_list[0]

def noise_removal(signal_noise_ratio, spectrum):
    # remove noise using average of 5% lowest intensity peaks as noise level and user defined signal/noise ratio
    if len(spectrum['m/z array']) == 0:
        spectrum_dataframe = spectrum_to_dataframe(spectrum)
        spectrum_dataframe = spectrum_dataframe.drop(spectrum_dataframe[spectrum_dataframe['intensity'] == 0].index)
        if not spectrum_dataframe.empty:
            sorted_intensities = sorted(spectrum_dataframe['intensity'].values.tolist())
            lowest_sorted_intensities = sorted_intensities[:int(round((len(sorted_intensities)*0.05)))]
            if len(lowest_sorted_intensities) == 0:
                spectrum_noise = sum(lowest_sorted_intensities)/len(lowest_sorted_intensities)
                print spectrum_noise
                spectrum_dataframe = spectrum_dataframe.drop(spectrum_dataframe[spectrum_dataframe['intensity'] <= (signal_noise_ratio * spectrum_noise)].index)
                spectrum['m/z array'] = spectrum_dataframe['m/z'].values
                spectrum['intensity array'] = spectrum_dataframe['intensity'].values
    return spectrum

def blank_removal(sample_spectrum, control_spectrum, peak_mz_tolerance):
    # remove control spectrum peaks from sample spectrum if m/z within specified tolerance
    # returns processed sample spectrum and dictionary with removed peaks
    sample_spectrum_df = spectrum_to_dataframe(sample_spectrum)
    control_spectrum_df = spectrum_to_dataframe(control_spectrum)
    if not sample_spectrum_df.empty and not control_spectrum_df.empty:
        blank_subtracted_df = pandas.merge_asof(sample_spectrum_df, control_spectrum_df, on='m/z',
                                                tolerance=peak_mz_tolerance, direction='nearest').fillna(0)
        changed_peaks = blank_subtracted_df[blank_subtracted_df['intensity_y'] != 0.0]
        changed_dict = copy.deepcopy(sample_spectrum)
        changed_dict['m/z array'] = changed_peaks['m/z'].values
        changed_dict['intensity array'] = changed_peaks['intensity_x'].values
        blank_subtracted_df = blank_subtracted_df.drop(blank_subtracted_df[blank_subtracted_df['intensity_y']
                                                                           != 0.0].index)
        sample_spectrum['m/z array'] = blank_subtracted_df['m/z'].values
        sample_spectrum['intensity array'] = blank_subtracted_df['intensity_x'].values
        return [sample_spectrum, changed_dict]
    else:
        return [sample_spectrum, None]

def spectra_compare(args, control_dataset, sample_spectrum):
    # select control spectrum and remove blanks
    ms_mode = sample_spectrum['msLevel']
    ret_time = sample_spectrum['retentionTime'] * 60
    if ms_mode > 1:
        precursor_mz = sample_spectrum['precursorMz'][0]['precursorMz']
    else:
        precursor_mz = None
    control_spectrum = select_control_spectrum(args, ms_mode, ret_time, precursor_mz, control_dataset)
    if control_spectrum == None:
        return [sample_spectrum, None]
    elif ms_mode != 2:
        processed_spectrum = blank_removal(sample_spectrum, control_spectrum, args['peak_mz_tolerance'])
        return processed_spectrum
        # processed_spectrum = [sample_spectrum, changed_dict]

def mgf_writer(spectrum_data_dict, output_dir, datatype):
    with open(output_dir + datatype + '_data_ms2.mgf', 'a') as mgf_file:
        if spectrum_data_dict['msLevel'] >= 2:
            mgf_file.write("BEGIN IONS" + "\n")
            mgf_file.write("TITLE=" + output_dir.split("\\")[-1].split("_blanka_")[0] + "." + spectrum_data_dict['num']
                       + "." + spectrum_data_dict['num'] + '. File:"' + output_dir.split("\\")[-1].split("_blanka_")[0]
                       + '.RAW", NativeID:"controllerType=0 controllerNumber=1 scan=' + spectrum_data_dict['num'] + '"'
                       + "\n")
            mgf_file.write("RTINSECONDS=" + str(spectrum_data_dict['retentionTime'] * 60) + "\n")
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
