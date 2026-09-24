import os
import sys
import json
#import jsonschema
import pathlib
import shutil
from datetime import datetime
import re
import multiprocessing as mp
import xds_input
import subprocess as sub

class Reprocess(object):

      def __init__(self, jData):
          self._ioDict = dict()
          self._ioDict['inData'] = json.dumps(jData, default=str)
          self._ioDict['outData'] = json.dumps(dict(), default=str)
          self._ioDict['success'] = True
          self._ioDict['WorkingDirectory'] = None
          self.old_folders = []
          return

      def get_inData(self):
          return json.loads(self._ioDict['inData'])

      def set_inData(self, jData):
          self._ioDict['inData'] = json.dumps(jData, default=str)

      jshandle = property(get_inData, set_inData)

      def get_outData(self):
          return json.loads(self._ioDict['outData'])

      def set_outData(self, results):
          self._ioDict['outData'] = json.dumps(results, default=str)

      results = property(get_outData, set_outData)

      def writeInputData(self, inData):
          # Write input data
          if self._ioDict['WorkingDirectory'] is not None:
              jsonName = "inData_" + self.__class__.__name__ + ".json"
              with open(str(self.getOutputDirectory() / jsonName), 'w') as f:
                  f.write(json.dumps(inData, default=str, indent=4))
          return

      def writeOutputData(self, results):
          self.set_outData(results)
          if self._ioDict['WorkingDirectory'] is not None:
              jsonName = "outData_" + self.__class__.__name__ + ".json"
              with open(str(self.getOutputDirectory() / jsonName), 'w') as f:
                  f.write(json.dumps(results, default=str, indent=4))
          return

      def setFailure(self):
          self._ioDict['success'] = False

      def is_success(self):
          return self._ioDict['success']

      def setOutputDirectory(self, path=None):
          if path is None:
              directory = self.jshandle.get('processing_directory', os.getcwd())
              self._ioDict['WorkingDirectory'] = pathlib.Path(directory)
          else:
              self._ioDict['WorkingDirectory'] = pathlib.Path(path)

      def getOutputDirectory(self):
          return self._ioDict['WorkingDirectory']

      @property
      def datadir_search(self):
          datadir = pathlib.Path(self.jshandle['raw_data'])
          listofImageFolders = []
          if datadir.exists():
             imagelist = os.listdir(str(datadir))
             if len(imagelist) > 0:
                for imagefolder in imagelist:
                    listofImageFolders.append(str(datadir / imagefolder))
          else:
             self.setFailure()
          return listofImageFolders

      def xdsfolder_search(self):
          procdir = pathlib.Path(self.jshandle['process_data'])
          old_folders = []
          if procdir.exists():
             old_folders = list(procdir.glob(self.jshandle['prefix'] + self.jshandle['suffix']))
             if len(old_folders) == 0:
                print('No data or processing folder found! plesae check the paths')
             self.setFailure()
          return old_folders

      def makeOutputDirectory(self):
          self.setOutputDirectory()
          new_output = self.getOutputDirectory() / datetime.now().strftime('proc_%Y%m%d')
          try:
              new_output.mkdir(parents=True)
          except FileExistsError:
              pass
          os.chdir(str(new_output))
          self.setOutputDirectory(path=str(new_output))

          old_folders = self.xdsfolder_search()
          listofxdsdirs = []
          listofLinks = []

          if len(old_folders) > 0:
             for ii in range(len(old_folders)):
                 old_folder_path = pathlib.Path(old_folders[ii])
                 folder_structure = old_folder_path.parents
                 xtal_name = 'xtal_%d' %ii
                 new_xds_dir = self.getOutputDirectory() / folder_structure[3].name / folder_structure[2].name / folder_structure[1].name / xtal_name
                 new_xds_dir.mkdir(parents=True)
                               
                 old_xds = old_folders[ii] / 'XDS.INP'
                 if old_xds.exists() and folder_structure[0].exists():
                    new_xds_inp = new_xds_dir / 'XDS.INP_old'
                    shutil.copy(str(old_xds), str(new_xds_inp))
                    # print("%s old_xds --> new path %s" % (old_xds, xtal_name))                                 
                    listofxdsdirs.append(new_xds_dir)
                    folder_raw = str(folder_structure[0]).replace("PROCESSED_DATA", "RAW_DATA")
                    # listofLinks.append(folder_structure[1] / 'links')
                    listofLinks.append(folder_raw)
                 else:
                    pass

          else:
              self.setFailure()
          return listofxdsdirs, listofLinks

      def DictionaryOldXDS(self, XDSINP):
          fname = pathlib.Path(XDSINP)
          xds_inData = dict()
          if not fname.exists():
             print('%s file does not exist; cannot create XDS.INP file' %fname)
             return xds_inData
          if fname.stat().st_size == 0:
             print('%s file exists but empty' %fname)
             
             return xds_inData
          rex_dict = dict(template=re.compile(r'NAME_TEMPLATE_OF_DATA_FRAMES=(?P<template>.*)\n'),
               detZ=re.compile(r'DETECTOR_DISTANCE=\s(?P<detZ>([0-9.]+))\n'),
               data_range=re.compile(r'DATA_RANGE=\s(?P<data_range>.*)\n'),
               beam=re.compile(r'ORGX=\s(?P<beam>.*)\n'),
               osc_range=re.compile(r'OSCILLATION_RANGE=\s(?P<osc_range>([0-9.]+))\n'))
          with open(str(fname), 'r') as fh:
               for line in fh:
                   for k, pat in rex_dict.items():
                       match = pat.search(line)
                       if match:
                          if k == 'template':
                             template_old = match.group('template')
                             xds_inData['template_name'] = pathlib.Path(template_old).name
                             print(f"new xds.inp file: {xds_inData['template_name']}")
                          if k == 'detZ':
                             xds_inData['detZ'] = match.group('detZ')
                          if k == 'data_range':
                             xds_inData['data_range'] = match.group('data_range')
                          if k == 'beam':
                             r = match.group('beam')
                             r = r.split()
                             xds_inData['beamX'] = r[0]
                             xds_inData['beamY'] = r[2]
                          if k == 'osc_range':
                             xds_inData['osc_range'] = match.group('osc_range')

          
          xds_inData['SG'] = self.jshandle.get('space_group', '0')
          xds_inData['unit_cell'] = self.jshandle.get('unit_cell', '79 79 38 90 90 90')
          xds_inData['res_cut'] = self.jshandle.get('resolution_cutoff', '0.0')
          xds_inData['friedel'] = self.jshandle.get('friedels_law', 'TRUE')
          xds_inData['ref'] = self.jshandle.get('reference_data', '')
          xds_inData['njobs'] = 4
          xds_inData['nproc'] = 10
          # xds_inData['nodes'] = 'mxhpc2-1704.esrf.fr mxhpc2-1705.esrf.fr'

          return xds_inData

      def xds_index_inp(self, XDSINP_old, linkname):
          xds_inData = self.DictionaryOldXDS(XDSINP_old)
          if not xds_inData:
             pass
             return
          xds_inData['res_cut'] = self.jshandle.get('index_res','5.0')
          xds_inData['jobs'] = 'XYCORR INIT COLSPOT IDXREF'
          try:
             xds_inData['template'] = os.path.join(linkname, xds_inData['template_name'])
             xds_string = xds_input.INP[self.jshandle.get('beamline', 'ID23-2')]
             new_xds_path = pathlib.Path(os.path.join(os.getcwd(), 'XDS.INP'))
             if not new_xds_path.exists():
                fh = open(str(new_xds_path), 'w')
             
                fh.write(xds_string.format(**xds_inData))
                fh.close()
             else:
                pass
          except Exception as err:
             print(err)
             pass
          return

      def xds_integrate_inp(self, XDSINP_old, linkname):
          xds_inData = self.DictionaryOldXDS(XDSINP_old)
          if not xds_inData:
             pass
             return
          xds_inData['res_cut'] = self.jshandle.get('resolution_cutoff', '0.0')
          xds_inData['jobs'] = 'DEFPIX INTEGRATE CORRECT'
          xds_inData['template'] = os.path.join(linkname, xds_inData['template_name'])
          xds_string = xds_input.INP[self.jshandle.get('beamline', 'ID23-2')]
          
          new_xds_path = pathlib.Path(os.path.join(os.getcwd(), 'XDS.INP'))
          if new_xds_path.exists():
             shutil.copy('XDS.INP', 'Indexing.INP')
             try:
               fh = open(str(new_xds_path), 'w')
               fh.write(xds_string.format(**xds_inData))
               fh.close()
             except Exception as err:
               print(err)
               pass
          else:
             pass
          return

      def runxds(self, xdsdir, linkname):
          '''
          listofxdsdirs = self.makeOutputDirectory()
          if len(listofxdsdirs) == 0:
             print("old processing path was not found")
             return
          '''
          try:
              old_xds = xdsdir / 'XDS.INP_old'
              os.chdir(str(xdsdir))
              self.xds_index_inp(old_xds, linkname)
              sub.call(['xds_par > /dev/null'], shell=True)
              
              self.xds_integrate_inp(old_xds, linkname)
              sub.call(['xds_par > /dev/null'], shell=True)
              os.chdir(str(self.getOutputDirectory()))
              HKLfile = pathlib.Path(xdsdir / 'XDS_ASCII.HKL')
              if not HKLfile.exists():
                 print("Failed to process %s" %xdsdir)
              else:
                 print('processed')
          except OSError:
             self.setFailure()
             raise

          return
    
      def runpilatus(self):
          xdsjoblist, listofImageLinks = self.makeOutputDirectory()
          if len(xdsjoblist) == 0:
             print("old processing path was not found")
             return
          else:
             print("No. of xtals found in total: %d" %len(xdsjoblist))

          try:
             for j in range(len(xdsjoblist)):

                proc = []; 
                for i in range (0,5):
                    try:
                       jobid = mp.Process(target=self.runxds(xdsjoblist[(j*5)+i], listofImageLinks[(j*5)+i]))
                       proc.append(jobid)

                    except IndexError:
                       pass
                for p in proc:
                    p.start()

                for p in proc:
                    p.join()
                print("%d crystals have been attempted\n" %(j+1))

          except Exception:
             self.setFailure()
             raise
          return

if __name__ == '__main__':
   jdata = open(sys.argv[1], 'r')
   rep = Reprocess(json.load(jdata))
   rep.runpilatus()
   

