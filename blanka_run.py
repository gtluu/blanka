import argparse, subprocess, os, copy, sys, pandas, numpy, timeit
import pyteomics.mzxml as pytmzxml
import pyteomics.mgf as pytmgf
from multiprocessing import Pool, cpu_count
from functools import partial
import blanka_lcms as lcms
import blanka_maldi_dd as dd

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cpu', help="number of threads to use - default = max-1", default=cpu_count() - 1, type=int)
    parser.add_argument('--instrument', help="instrument/experiment: choose 'lcq', 'qtof', or 'dd'",
                        default='', type=str)
    parser.add_argument('--sample', help="sample input directory with single or multiple files", type=str, default='')
    parser.add_argument('--control', help="control input file path with '.mzXML' file extension (lcq/qtof) or name of \
                                           control sample spot (dd)", type=str, default='')
    parser.add_argument('--dd_template', help='dried droplet excel sheet with sample names', default='', type=str)
    parser.add_argument('--output', help="output directory for all generated files; default=source folder",
                        type=str, default='')
    parser.add_argument('--signal_noise_ratio', help="integer signal to noise ratio - default = 4", default=4, type=int)
    parser.add_argument('-r', '--retention_time_tolerance', help="retention time error in seconds - default = 0.1 s",
                        default=0.1, type=float)
    parser.add_argument('-m', '--precursor_mz_tolerance', help="absolute precursor m/z error in Da - default = 0.02 Da",
                        default=0.02, type=float)
    parser.add_argument('-p', '--peak_mz_tolerance', help="absolute peak m/z error in Da - default = 0.02 Da",
                        default=0.02, type=float)
    parser.add_argument('--noise_removal_only', help='only perform noise removal; no blank removal',
                        default=False, type=bool)
    parser.add_argument('--blank_removal_only', help='only perform blank removal; no noise removal',
                        default=False, type=bool)
    arguments = parser.parse_args()
    return vars(arguments)

def run_lcms(args):
    if not args['sample'].endswith('.mzXML'):
        sample_file_list = lcms.mzxml_data_detection(args['sample'])
        # find .mzXML files in sample and control directory
        if sample_file_list == []:
            sample_raw_list = lcms.raw_data_detection(args, args['sample'])
            lcms.msconvert(args, sample_raw_list)
            # detect raw data if no .mzXML files found and convert to .mzXML
            if args['output'] == '':
                files_dir = args['sample']
            else:
                files_dir = args['output']
            sample_file_list = lcms.mzxml_data_detection(files_dir)
            # list of .mzXML files
    else:
        sample_file_list = [args['sample']]
        # single .mzXML file
    control_data = lcms.load_control_data(args)
    # list of control spectra from .mzXML file
    control_noise_args = partial(lcms.noise_removal, args['signal_noise_ratio'])
    control_noiseless_data = pool.map(control_noise_args, control_data, chunksize=200)
    # list of noise removed control spectra

    if args['noise_removal_only'] == False and args['blank_removal_only'] == False:
        for dataset in sample_file_list:
            if args['control'] != dataset:
                if args['output'] == '':
                    blanka_output = dataset.split('.')[0] + '_blanka_'
                    # ex: D:\folder\filename
                else:
                    blanka_output = args['output'] + dataset.split('\\')[-1].split('.')[0] + '_blanka_'
                # prep output directory/filenames
                print "Processing " + dataset.split("\\")[-1]
                sample_data = list(pytmzxml.read(dataset))
                print "Removing Noise"
                sample_noise_args = partial(lcms.noise_removal, args['signal_noise_ratio'])
                sample_noiseless_data = pool.map(sample_noise_args, sample_data, chunksize=200)
                for spectrum in sample_noiseless_data:
                    lcms.mgf_writer(spectrum, blanka_output, 'noise_removed')
                # remove noise and write to .mgf
                print "Removing Blank"
                spectra_compare_args = partial(lcms.spectra_compare, args, control_noiseless_data)
                processed_data = filter(None, pool.map(spectra_compare_args, sample_noiseless_data, chunksize=200))
                for processed_spectrum, changed_spectrum_data in processed_data:
                    lcms.mgf_writer(processed_spectrum, blanka_output, 'processed')
                    if changed_spectrum_data != None:
                        lcms.mgf_writer(changed_spectrum_data, blanka_output, 'removed_peaks')
                # remove blank and write to .mgf
    elif args['noise_removal_only'] == True:
        for dataset in sample_file_list:
            if args['control'] != dataset:
                if args['output'] == '':
                    blanka_output = dataset.split('.')[0] + '_blanka_'
                    # ex: D:\folder\filename
                else:
                    blanka_output = args['output'] + dataset.split('\\')[-1].split('.')[0] + '_blanka_'
                # prep output directory/filenames
                print "Processing " + dataset.split("\\")[-1]
                sample_data = list(pytmzxml.read(dataset))
                print "Removing Noise"
                sample_noise_args = partial(lcms.noise_removal, args['signal_noise_ratio'])
                sample_noiseless_data = pool.map(sample_noise_args, sample_data, chunksize=200)
                for spectrum in sample_noiseless_data:
                    lcms.mgf_writer(spectrum, blanka_output, 'noise_removed')
                # remove noise and write to .mgf
    elif args['blank_removal_only'] == True:
        for dataset in sample_file_list:
            if args['control'] != dataset:
                if args['output'] == '':
                    blanka_output = dataset.split('.')[0] + '_blanka_'
                    # ex: D:\folder\filename
                else:
                    blanka_output = args['output'] + dataset.split('\\')[-1].split('.')[0] + '_blanka_'
                # prep output directory/filenames
                print "Processing " + dataset.split("\\")[-1]
                sample_data = list(pytmzxml.read(dataset))
                print "Removing Blank"
                spectra_compare_args = partial(lcms.spectra_compare, args, control_data)
                processed_data = filter(None, pool.map(spectra_compare_args, sample_data, chunksize=200))
                for processed_spectrum, changed_spectrum_data in processed_data:
                    lcms.mgf_writer(processed_spectrum, blanka_output, 'processed')
                    if changed_spectrum_data != None:
                        lcms.mgf_writer(changed_spectrum_data, blanka_output, 'removed_peaks')
                # remove blank and write to .mgf

def run_maldi_dd(args):
    if not args['sample'].endswith('.mzXML'):
        file_list = dd.mzxml_data_detection(args['sample'])
        # find .mzXML files in sample directory
        if file_list == []:
            raw_file_list = dd.raw_data_detection(args)
            #raw_file_list = dd.parse_maldi_template(args, raw_file_list)
            raw_file_list = [(i, i.split('\\')[-5]) for i in raw_file_list]
            dd.msconvert(args, raw_file_list)
            # detect raw data if no .mzXML files found and convert to .mzXML
            if args['output'] == '':
                files_dir = args['sample']
            else:
                files_dir = args['output']
            file_list = dd.mzxml_data_detection(files_dir)
            # list of .mzXML files
    else:
        file_list = [args['sample']]
        # single .mzXML file

    control_data = dd.load_control_data(args)
    control_list = control_data[0][1]
    control_data = [i[0] for i in control_data]
    # list of control spectra (dict form)
    control_noise_args = partial(dd.noise_removal, args['signal_noise_ratio'])
    control_noiseless_data = pool.map(control_noise_args, control_data)
    # list of noise removed control spectra
    control_noiseless_data = dd.combine_control_spectra(control_noiseless_data)
    # single control dict

    if args['noise_removal_only'] == False and args['blank_removal_only'] == False:
        file_list = [i for i in file_list if not i.startswith(args['control']) and i not in control_list]
        sample_data = [[list(pytmzxml.read(mzxml))[0], mzxml] for mzxml in file_list]
        print "Removing noise from samples."
        sample_noise_args = partial(dd.noise_removal, args['signal_noise_ratio'])
        sample_noiseless_data = pool.map(sample_noise_args, sample_data)
        for spectrum, filename in sample_noiseless_data:
            if args['output'] == '':
                blanka_output = filename.split('.')[0] + '_blanka_'
                # ex: D:\folder\filename
            else:
                blanka_output = args['output'] + filename.split('\\')[-1].split('.')[0] + '_blanka_'
            dd.mgf_writer(spectrum, blanka_output, 'noise_removed')
        print "Removing blank from samples."
        sample_blank_args = partial(dd.blank_removal, args['peak_mz_tolerance'], control_noiseless_data)
        sample_blankless_data = pool.map(sample_blank_args, sample_noiseless_data)
        for spectrum, changed_spectrum_data in sample_blankless_data:
            if args['output'] == '':
                blanka_output = spectrum[1].split('.')[0] + '_blanka_'
                # ex: D:\folder\filename
            else:
                blanka_output = args['output'] + spectrum[1].split('\\')[-1].split('.')[0] + '_blanka_'
            dd.mgf_writer(spectrum[0], blanka_output, 'processed')
            if changed_spectrum_data != None:
                dd.mgf_writer(changed_spectrum_data, blanka_output, 'removed_peaks')
    elif args['noise_removal_only'] == True:
        file_list = [i for i in file_list if not i.startswith(args['control']) and i not in control_list]
        sample_data = [[list(pytmzxml.read(mzxml))[0], mzxml] for mzxml in file_list]
        print "Removing noise from samples."
        sample_noise_args = partial(dd.noise_removal, args['signal_noise_ratio'])
        sample_noiseless_data = pool.map(sample_noise_args, sample_data)
        for spectrum, filename in sample_noiseless_data:
            if args['output'] == '':
                blanka_output = filename.split('.')[0] + '_blanka_'
                # ex: D:\folder\filename
            else:
                blanka_output = args['output'] + filename.split('\\')[-1].split('.')[0] + '_blanka_'
            dd.mgf_writer(spectrum, blanka_output, 'noise_removed')
    elif args['blank_removal_only'] == True:
        file_list = [i for i in file_list if not i.startswith(args['control']) and i not in control_list]
        sample_data = [[list(pytmzxml.read(mzxml))[0], mzxml] for mzxml in file_list]
        print "Removing blank from samples."
        sample_blank_args = partial(dd.blank_removal, args['peak_mz_tolerance'], control_noiseless_data)
        sample_blankless_data = pool.map(sample_blank_args, sample_data)
        for spectrum, changed_spectrum_data in sample_blankless_data:
            if args['output'] == '':
                blanka_output = spectrum[1].split('.')[0] + '_blanka_'
                # ex: D:\folder\filename
            else:
                blanka_output = args['output'] + spectrum[1].split('\\')[-1].split('.')[0] + '_blanka_'
            dd.mgf_writer(spectrum[0], blanka_output, 'processed')
            if changed_spectrum_data != None:
                dd.mgf_writer(changed_spectrum_data, blanka_output, 'removed_peaks')

if __name__ == "__main__":

    arguments = get_args()
    pool = Pool(processes=arguments['cpu'])

    if arguments['instrument'] == 'lcq' or arguments['instrument'] == 'qtof':
        run_lcms(arguments)
    elif arguments['instrument'] == 'dd':
        run_maldi_dd(arguments)

    pool.close()
    pool.join()
