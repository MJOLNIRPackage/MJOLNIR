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
import re, numpy as np
from MJOLNIR.Data.DataSet import OxfordList

DEBUG = False

possibleArguments = ['test', 'tutorials', 'wheel', 'html', 'upload', 'testVersion', 'version']

parser = argparse.ArgumentParser(description="Make tool replacing linux Makefile.")
parser.add_argument("task", nargs='?', default='help', type=str, help="Type of task to be performed. Possible tasks are: {}. Run without argument to see help menu.".format(OxfordList(possibleArguments)))
parser.add_argument("version", nargs='?',default=None, type=str, help="Version number to be created. Default is previous plus 1, i.e. x.y.z -> x.y.z+1")

def callHelp(parser=parser):
    parser.print_help()

def test():
    os.system("python -m pytest -vv test")

# def cleanTutorialFolders():
#     folders = [f.path for f in os.scandir(os.path.join('docs','Tutorials')) if f.is_dir()]
#     print('Cleaning tutorial folders: '+OxfordList(folders))
#     for folder in folders:
#         os.system('python clean.py '+folder)

def generateTutorials():
    if not DEBUG:
        os.system('python '+str(os.path.join('docs','tutorials.py')))
    else:
        print('python '+str(os.path.join('docs','tutorials.py')))

def makeHTML():
    print('Creating HTML through Sphinx')
    if not DEBUG:
        os.system("sphinx-build docs build")

def makeWheel():
    if not DEBUG:
        os.system("python -m build")
    else:
        print("python -m build")

def getLatestVersions():
    list_of_files = glob.glob(os.path.join('dist','*')) # * means all if need specific format then *.csv
    list_of_files = sorted(list_of_files, key=os.path.getctime)
    
    versions = {}
    pattern = r"-([0-9]*)\.([0-9]*)\.([0-9]*)\.?(post)?([0-9]*)?"

    for b in list_of_files:
        #b = 'dist\\MJOLNIR-1.1.19.tar.gz'
        a = re.findall(pattern,b.lower())[0]
        if not a[0] == '1':
            print(b,a)
        else:
            a = np.asarray(a)
            a = a[np.asarray([len(x)>0 for x in a])]
            a = '.'.join(a)
            a = a.replace('post.','post')
            if not a in versions:
                versions[a] = [os.path.getctime(os.path.join(r'C:\Users\lass_j\Documents\Software\MJOLNIR',b)),[b]]
            else:
                versions[a][1].append(b)
    sortedVersions = dict(sorted(versions.items(), key=lambda item: -item[1][0]))
    return sortedVersions

def getLatestBuilds():
    sortedVersions = getLatestVersions()
    key  = list(sortedVersions.keys())[0]
    return sortedVersions[key][1]

def getCurrentVersion():
    sortedVersions = getLatestVersions()
    string  = list(sortedVersions.keys())[0]
    #versionList = string.split('-')[1].split('.')[:-2]
    # if 'post' in versionList[-1]:
    #     versionList[-1] = versionList[-1].replace('post','')
    #     adder = 'post'
    # else:
    #     adder = ''
    # version = [int(x) for x in versionList]
    
    # if adder != '':
    #     version[-1]=adder+str(version[-1]+1)
    # else:
    #     version[-1]+=1
    return string#'.'.join([str(x) for x in version])

def uploadPyPi(server='testpypi'):
    build = getLatestBuilds()
    for b in build:
        print('Uploading ' + str(b) + ' to '+server)
        if not DEBUG:
            os.system('twine upload {} -r {}'.format(b,server))


def upload():
    uploadPyPi(server='testpypi')
    uploadPyPi(server='pypi')

def update(version):
    if not DEBUG:
        os.system('python Update.py '+version)
    else:
        print('python Update.py '+version)

def makeTutorials():
    # cleanTutorialFolders()
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
    if not DEBUG:
        os.system("git add {}".format(' '.join(addFiles)))
        os.system('git commit -m "Update version"')
        os.system('git tag -a {0} -m "{0}"'.format(version))
    else:
        print("git add {}".format(' '.join(addFiles)))
        print('git commit -m "Update version"')
        print('git tag -a {0} -m "{0}"'.format(version))
    makeWheel()
    if not DEBUG:
        os.system("git push")
        os.system("git push --tags")
    else:
        print("git push")
        print("git push --tags")

else:
    print('Provided argument not understood. Received',args.task,'\n\n')
    callHelp()