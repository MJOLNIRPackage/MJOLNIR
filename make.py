#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 11:54:33 2021

@author: Jakob Lass

Tool to replace makefile in order to fully port development to windows
"""

import argparse
from ntpath import join

import os,glob

from MJOLNIR.Data.DataSet import OxfordList

possibleArguments = ['test', 'tutorials', 'wheel', 'html', 'upload', 'testVersion', 'version']

parser = argparse.ArgumentParser(description="Make tool replacing linux Makefile.")
parser.add_argument("task", nargs='?', default='help', type=str, help="Type of task to be performed. Possible tasks are: {}. Run without argument to see help menu.".format(OxfordList(possibleArguments)))
parser.add_argument("version", nargs='?',default=None, type=str, help="Version number to be created. Default is previous plus 1, i.e. x.y.z -> x.y.z+1")

def callHelp(parser=parser):
    parser.print_help()

def test():
    os.system("python -m pytest -vv test")

def cleanTutorialFolders():
    folders = [f.path for f in os.scandir(os.path.join('docs','Tutorials')) if f.is_dir()]
    print('Cleaning tutorial folders: '+OxfordList(folders))
    for folder in folders:
        os.system('python clean.py '+folder)

def generateTutorials():
    os.system('python '+str(os.path.join('Tutorials','tutorials.py')))
    

def makeHTML():
    print('Creating HTML through Sphinx')
    os.system("sphinx-build docs build")

def makeWheel():
    os.system("python setup.py sdist")

def getLatestBuild():
    list_of_files = glob.glob(os.path.join('dist','*')) # * means all if need specific format then *.csv
    return max(list_of_files, key=os.path.getctime)

def getCurrentVersion():
    string  = getLatestBuild()
    versionList = string.split('-')[1].split('.')[:-2]
    if 'post' in versionList[-1]:
        versionList[-1] = versionList[-1].replace('post','')
        adder = 'post'
    else:
        adder = ''
    version = [int(x) for x in versionList]
    
    if adder != '':
        version[-1]=adder+str(version[-1]+1)
    else:
        version[-1]+=1
    return '.'.join([str(x) for x in version])

def uploadPyPi(server='testpypi'):
    build = getLatestBuild()
    print('Uploading ' + str(build) + ' to '+server)
    os.system('twine upload {} -r {}'.format(build,server))


def upload():
    uploadPyPi(server='testpypi')
    uploadPyPi(server='pypi')

def update(version):
    os.system('python Update.py '+version)

def makeTutorials():
    cleanTutorialFolders()
    generateTutorials()
    makeHTML()

args = parser.parse_args()


if args.task == '' or args.task.lower() == 'help':
    callHelp()
elif args.task.lower() == 'test':
    test()

elif args.task.lower() == 'tutorials':
    makeTutorials()

elif args.task.lower() == 'wheel':
    makeWheel()

elif args.task.lower() == 'upload':
    upload()

elif args.task.lower() == 'testversion':
    if args.version is None:
        version = getCurrentVersion()
    else:
        version = args.version
    print('Creating test version '+version)
    update(version=version)
    makeWheel()
    uploadPyPi(server='testpypi')

elif args.task.lower() == 'html':
    makeHTML()

elif args.task.lower() == 'version':
    if args.version is None:
        version = getCurrentVersion()
    else:
        version = args.version
    print('Creating version '+version)
    update(version=version)
    makeTutorials()
    addFiles = ['setup.py',
                os.path.join('docs','conf.py'),
                os.path.join('docs','Tutorials','*'),
                os.path.join('docs','index.rst'),
                os.path.join('MJOLNIR','__init__.py'),
                ]
    os.system("git add {}".format(' '.join(addFiles)))
    os.system('git commit -m "Update version"')
    os.system('git tag -a {:} -m "{:}"'.format(version))
    makeWheel()
    os.system("git push")
    os.system("git push --tags")

else:
    print('Provided argument not understood. Recieved',args.task,'\n\n')
    callHelp()