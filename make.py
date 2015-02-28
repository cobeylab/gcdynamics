#!/usr/bin/env python

import os
import sys
import subprocess

EXEC_NAME = 'gcdynamics'
COMPILERS = [
	('icc', 'icpc', '-std=c++11'),
	('clang', 'clang++', '-std=c++11 -stdlib=libc++'),
	('gcc', 'g++', '-std=c++11')
]

BOOST_VERSION = ('1', '57', '0')
BOOST_URL = 'http://sourceforge.net/projects/boost/files/boost/{0}/boost_{1}.tar.bz2/download'.format(
	'.'.join(BOOST_VERSION),
	'_'.join(BOOST_VERSION)
)

def execPresent(execName):
	return os.system('which {0}'.format(execName)) == 0

def run_cmd(cmd_str):
	print(cmd_str)
	proc = subprocess.Popen(cmd_str, shell=True)
	proc.wait()

def download_boost():
	print('Checking for boost...')
	
	if os.path.exists('boost'):
		print('boost appears to be present.\n')
		assert os.path.isdir('boost')
		return
	
	boost_version_underscore = '_'.join(BOOST_VERSION)
	
	boost_basename = 'boost_{0}'.format('_'.join(BOOST_VERSION))
	boost_filename = '{0}.tar.bz2'.format(boost_basename)
	
	print('Downloading boost...')
	if os.path.exists(boost_filename):
		run_cmd('rm {0}'.format(boost_filename))
	run_cmd('curl -L {0} -o boost_{1}.tar.bz2'.format(BOOST_URL, boost_version_underscore))
	
	print('Unpacking boost...')
	if os.path.exists(boost_basename):
		run_cmd('rm -r {0}'.format(boost_basename))
	run_cmd('tar xzf {0}'.format(boost_filename))
	run_cmd('rm {0}'.format(boost_filename))
	run_cmd('mv {0} boost'.format(boost_basename))
	
	print('')

def git_submodules():
	print('Updating git submodules...')
	run_cmd('git submodule init')
	run_cmd('git submodule update')
	print('')

if __name__ == '__main__':
	os.chdir(os.path.dirname(__file__))
	
	download_boost()
	git_submodules()
	
	print('Looking for compilers...')
	for compiler in COMPILERS:
		if execPresent(compiler[0]) and execPresent(compiler[1]):
			cCompiler, cppCompiler, flags = compiler
			break
	print('')
	
	includeDirs = [
		'boost/',
		'zppjson/src',
		'zppdb/src',
		'zppsim/src',
		'libjson'
	]
	
	srcDirs = [
		'src',
		'zppjson/src',
		'zppdb/src',
		'zppsim/src',
	]
	
	print('Compiling...')
	
	run_cmd('mkdir -p bin')
	run_cmd('{0} -O3 -c libjson/json.c -o bin/libjson.o'.format(cCompiler))
	run_cmd('{0} -O3 -lsqlite3 {1} {2} bin/libjson.o {3} -o bin/{4}'.format(
		cppCompiler,
		flags,
		' '.join(['-I{0}'.format(os.path.expanduser(x)) for x in includeDirs]),
		' '.join(['{0}/*.cpp'.format(x) for x in srcDirs]),
		EXEC_NAME
	))
	os.system('rm bin/libjson.o')
	
	print('Build complete.')
