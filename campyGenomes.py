#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
campyGenomes.py - Download and run INNUca in Campylobacter HTS public available data
<https://github.com/B-UMMI/campyGenomes/>

Copyright (C) 2016 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: August 25, 2016

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import utils
import sys
import random
import multiprocessing
import time
import shutil

version = '0.2'

general_threads_to_use = [4, 8, 16, 32, 64, 128]


def requiredPrograms():
	programs_version_dictionary = {}

	programs_version_dictionary['getSeqENA.py'] = ['--version', '>=', '0.4']
	programs_version_dictionary['ascp'] = ['--version', '>=', '3.6.1']

	programs_version_dictionary['INNUca.py'] = ['--version', '>=', '1.6']
	programs_version_dictionary['gunzip'] = ['--version', '>=', '1.6']
	programs_version_dictionary['java'] = ['-version', '>=', '1.8']
	programs_version_dictionary['mlst'] = ['--version', '>=', '2.4']
	missingPrograms = utils.checkPrograms(programs_version_dictionary)
	if len(missingPrograms) > 0:
		sys.exit('\n' + 'Errors:' + '\n' + '\n'.join(missingPrograms))


def removeFiles(directory, file_extention):
	files = [f for f in os.listdir(directory) if not f.startswith('.') and (os.path.isfile(os.path.join(directory, f)) or os.path.islink(os.path.join(directory, f)))]
	for file_found in files:
		if file_found.endswith(file_extention):
			os.remove(os.path.join(directory, file_found))


def downloadAndINNUca(outdir, run_ID, asperaKey, threads):
	start_time = time.time()
	temp_file = os.path.join(outdir, run_ID + '.temp.runID_fileList.txt')
	with open(temp_file, 'wt') as writer:
		writer.write(run_ID + '\n')

	command = ['getSeqENA.py', '-l', temp_file, '-o', outdir, '-a', asperaKey, '--downloadLibrariesType', 'PE']
	getSeqENA_run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

	os.remove(temp_file)

	sample_directory = os.path.join(outdir, run_ID, '')

	innuca_run_successfully = False
	if getSeqENA_run_successfully:
		command = ['INNUca.py', '-i', sample_directory, '-s', '"Campylobacter jejuni"', '-g', '1.6', '-o', sample_directory, '-j', str(threads), '--spadesSaveReport', '--pilonKeepSPAdesAssembly']
		innuca_run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None)

		innuca_dir = os.path.join(sample_directory, run_ID, '')
		files = [f for f in os.listdir(innuca_dir) if not f.startswith('.') and os.path.isfile(os.path.join(innuca_dir, f))]
		for file_innuca in files:
			shutil.move(os.path.join(innuca_dir, file_innuca), os.path.join(sample_directory, file_innuca))
		utils.removeDirectory(innuca_dir)

	removeFiles(sample_directory, '.gz')
	removeFiles(sample_directory, '.log')
	removeFiles(sample_directory, '.cpu.txt')

	if innuca_run_successfully:
		time_taken = utils.runTime(start_time)
		utils.saveVariableToPickle(time_taken, sample_directory, run_ID + '_downloadAndINNUca_time')

	utils.saveVariableToPickle(innuca_run_successfully, sample_directory, run_ID + '_run_successfully')


def determineNumberProcess(threads_to_use):
	number_process = {}
	threads_to_use = sorted(threads_to_use, reverse=True)
	for threads in threads_to_use:
		number_process[threads] = threads_to_use[0] // threads
	return number_process


def determineBatchSamples(samples, threads_to_use):
	base_batch = len(samples) // len(threads_to_use)
	remaining_samples = len(samples) % len(threads_to_use)

	samples_each_threads = {}

	base_position = 0
	for i in range(0, len(threads_to_use)):
		number_samples = base_batch
		if i < remaining_samples:
				number_samples += 1
		samples_each_threads[threads_to_use[i]] = samples[base_position: base_position + number_samples]
		base_position += number_samples

	return samples_each_threads


def runCampyGenomes(args):
	start_time = time.time()

	listRunIDs = utils.getListIDs(os.path.abspath(args.listRunIDs.name))
	outdir = os.path.abspath(args.outdir)
	utils.check_create_directory(outdir)
	asperaKey = args.asperaKey.name
	threads_to_use = [j for j in general_threads_to_use if j <= args.threads]

	# Start logger
	logfile, time_str = utils.start_logger(outdir)

	# Get general information
	utils.general_information(logfile, version, outdir, time_str)

	# Check programms
	requiredPrograms()

	# Randomize the list with Run IDs
	random.shuffle(listRunIDs)

	number_process = determineNumberProcess(threads_to_use)

	samples_each_threads = determineBatchSamples(listRunIDs, threads_to_use)

	run_successfully = 0
	with open(os.path.join(outdir, 'samples_with_problems.' + time_str + '.tab'), 'wt') as writer_success:
		with open(os.path.join(outdir, 'running_times.' + time_str + '.tab'), 'wt') as writer_times:

			for threads in samples_each_threads:
				print '\n' + 'Running for ' + str(threads) + ' threads' + '\n'
				threads_dir = os.path.join(outdir, str(threads) + '_threads', '')
				utils.check_create_directory(threads_dir)

				pool = multiprocessing.Pool(processes=number_process[threads])
				for sample in samples_each_threads[threads]:
					pool.apply_async(downloadAndINNUca, args=(threads_dir, sample, asperaKey, threads,))
				pool.close()
				pool.join()

				removeFiles(threads_dir, '.log')
				removeFiles(threads_dir, 'getSeqENA.samples_with_problems.txt')
				removeFiles(threads_dir, '.cpu.txt')

				samples_directories = [d for d in os.listdir(threads_dir) if not d.startswith('.') and os.path.isdir(os.path.join(threads_dir, d, ''))]
				for sample_dir in samples_directories:
					sample_dir_path = os.path.join(threads_dir, sample_dir, '')

					files = [f for f in os.listdir(sample_dir_path) if not f.startswith('.') and os.path.isfile(os.path.join(sample_dir_path, f))]
					for file_found in files:
						file_path = os.path.join(sample_dir_path, file_found)
						if file_found == sample_dir + '_run_successfully.pkl':
							sample_run_successfully = utils.extractVariableFromPickle(file_path)
							if not sample_run_successfully:
								writer_success.write(sample_dir + '\t' + threads_dir + '\n')
							else:
								run_successfully += 1
							os.remove(file_path)
						elif file_found == sample_dir + '_downloadAndINNUca_time.pkl':
							time_taken = utils.extractVariableFromPickle(file_path)
							writer_times.write(sample_dir + '\t' + threads_dir + '\t' + str(time_taken) + '\n')
							os.remove(file_path)

	time_taken = utils.runTime(start_time)
	del time_taken

	if run_successfully == 0:
		sys.exit('No RunIDs were successfully run!')
	else:
		print str(run_successfully) + ' samples out of ' + str(len(listRunIDs)) + ' run successfully'


def main():

	parser = argparse.ArgumentParser(prog='campyGenomes.py', description=" Download and run INNUca in Campylobacter HTS public available data", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

	parser_required = parser.add_argument_group('Required options')
	parser_required.add_argument('-l', '--listRunIDs', type=argparse.FileType('r'), metavar='/path/to/list/RunIDs.txt', help='Path to list containing the RunIDs to be downloaded', required=True)
	parser_required.add_argument('-a', '--asperaKey', type=argparse.FileType('r'), metavar='/path/to/asperaweb_id_dsa.openssh', help='Tells analyseSero38.py to download fastq files from ENA using Aspera Connect. With this option, the path to Private-key file asperaweb_id_dsa.openssh is provided (normaly found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).', required=True)

	parser_optional = parser.add_argument_group('Facultative options')
	parser_optional.add_argument('-o', '--outdir', type=str, metavar='/output/directory/', help='Path for output directory', required=False, default='.')
	parser_optional.add_argument('-j', '--threads', type=int, metavar='N', help='Maximum number of threads', required=False, default=4, choices=general_threads_to_use)

	parser.set_defaults(func=runCampyGenomes)

	args = parser.parse_args()

	args.func(args)


if __name__ == "__main__":
	main()
