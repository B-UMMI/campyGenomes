import os
import sys
import shlex
import subprocess
import os.path
import time
import shutil
from threading import Timer
import pickle


def start_logger(workdir):
	time_str = time.strftime("%Y%m%d-%H%M%S")
	sys.stdout = Logger(workdir, time_str)
	logfile = sys.stdout.getLogFile()
	return logfile, time_str


class Logger(object):
	def __init__(self, out_directory, time_str):
		self.logfile = os.path.join(out_directory, str('run.' + time_str + '.log'))
		self.terminal = sys.stdout
		self.log = open(self.logfile, "w")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)
		self.log.flush()

	def flush(self):
		pass

	def getLogFile(self):
		return self.logfile


def general_information(logfile, version, outdir, time_str):
	# Check if output directory exists

	print '\n' + '==========> campyGenomes <=========='
	print '\n' + 'Program start: ' + time.ctime()

	# Tells where the logfile will be stored
	print '\n' + 'LOGFILE:'
	print logfile

	# Print command
	print '\n' + 'COMMAND:'
	script_path = os.path.abspath(sys.argv[0])
	print sys.executable + ' ' + script_path + ' ' + ' '.join(sys.argv[1:])

	# Print directory where programme was lunch
	print '\n' + 'PRESENT DIRECTORY :'
	present_directory = os.path.abspath(os.getcwd())
	print present_directory

	# Print program version
	print '\n' + 'VERSION:'
	scriptVersionGit(version, present_directory, script_path)

	# Print PATH variable
	print '\n' + 'PATH variable:'
	print os.environ['PATH']

	# Save CPU information
	with open(os.path.join(outdir, 'cpu_information.' + time_str + '.cpu.txt'), 'wt') as reader:
		command = ['cat', '/proc/cpuinfo']
		run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
		reader.write(stdout)


def scriptVersionGit(version, directory, script_path):
	print 'Version ' + version
	os.chdir(os.path.dirname(script_path))
	command = ['git', 'log', '-1', '--date=local', '--pretty=format:"%h (%H) - Commit by %cn, %cd) : %s"']
	run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
	print stdout
	command = ['git', 'remote', 'show', 'origin']
	run_successfully, stdout, stderr = runCommandPopenCommunicate(command, False, None)
	print stdout
	os.chdir(directory)


def runTime(start_time):
	end_time = time.time()
	time_taken = end_time - start_time
	hours, rest = divmod(time_taken, 3600)
	minutes, seconds = divmod(rest, 60)
	print 'Runtime :' + str(hours) + 'h:' + str(minutes) + 'm:' + str(round(seconds, 2)) + 's'
	return time_taken


# USADO
def check_create_directory(directory):
	if not os.path.isdir(directory):
		os.makedirs(directory)


# USADO
def runCommandPopenCommunicate(command, shell_True, timeout_sec_None):
	run_successfully = False
	if isinstance(command, basestring):
		command = shlex.split(command)
	else:
		command = shlex.split(' '.join(command))

	print 'Running: ' + ' '.join(command)
	if shell_True:
		command = ' '.join(command)
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	else:
		proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if timeout_sec_None is None:
		stdout, stderr = proc.communicate()
	else:
		timer = Timer(timeout_sec_None, proc.kill)
		timer.start()
		stdout, stderr = proc.communicate()
		timer.cancel()

	if proc.returncode == 0:
		run_successfully = True
	else:
		print 'STDOUT'
		print stdout.decode("utf-8")
		print 'STDERR'
		print stderr.decode("utf-8")
	return run_successfully, stdout, stderr


# Remove directory
def removeDirectory(directory):
	if os.path.isdir(directory):
		shutil.rmtree(directory)


# USADO
def getListIDs(fileListIDs):
	list_ids = []

	with open(fileListIDs, 'rtU') as lines:
		for line in lines:
			line = line.splitlines()[0]
			if len(line) > 0:
				list_ids.append(line)

	if len(list_ids) == 0:
		sys.exit('No runIDs were found in ' + fileListIDs)

	return list_ids


# Check programs versions
def checkPrograms(programs_version_dictionary):
	print '\n' + 'Checking dependencies...'
	programs = programs_version_dictionary
	which_program = ['which', '']
	listMissings = []
	for program in programs:
		which_program[1] = program
		run_successfully, stdout, stderr = runCommandPopenCommunicate(which_program, False, None)
		if not run_successfully:
			listMissings.append(program + ' not found in PATH.')
		else:
			if programs[program][0] is None:
				print program + ' (impossible to determine programme version) found at: ' + stdout.splitlines()[0]
			else:
				check_version = [stdout.splitlines()[0], programs[program][0]]
				run_successfully, stdout, stderr = runCommandPopenCommunicate(check_version, False, None)
				if stdout == '':
					stdout = stderr
				if program == 'bunzip2':
					version_line = stdout.splitlines()[0].rsplit(',', 1)[0].split(' ')[-1]
				else:
					version_line = stdout.splitlines()[0].split(' ')[-1]
				replace_characters = ['"', 'v', 'V', '+']
				for i in replace_characters:
					version_line = version_line.replace(i, '')
				print program + ' (' + version_line + ') found'
				if programs[program][1] == '>=':
					program_found_version = version_line.split('.')
					program_version_required = programs[program][2].split('.')
					if float('.'.join(program_found_version[0:2])) < float('.'.join(program_version_required[0:2])):
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
					elif float('.'.join(program_found_version[0:2])) == float('.'.join(program_version_required[0:2])):
						if len(program_version_required) == 3:
							if len(program_found_version) == 2:
								program_found_version.append(0)
							if program_found_version[2].split('_')[0] < program_version_required[2]:
								listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
				else:
					if version_line != programs[program][2]:
						listMissings.append('It is required ' + program + ' with version ' + programs[program][1] + ' ' + programs[program][2])
	return listMissings


def saveVariableToPickle(variableToStore, outdir, prefix):
	pickleFile = os.path.join(outdir, str(prefix + '.pkl'))
	with open(pickleFile, 'wb') as writer:
		pickle.dump(variableToStore, writer)


def extractVariableFromPickle(pickleFile):
	with open(pickleFile, 'rb') as reader:
		variable = pickle.load(reader)
	return variable
