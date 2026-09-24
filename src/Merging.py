# -*- coding: utf-8 -*-

__author__ = ["S. Basu"]
__license__ = "M.I.T"
__date__ = "26/06/2017"
__refactordate__ = "06/05/2021"



import errno
import shutil
import time
import subprocess as sub
import pathlib
import jsonschema

from src.cellprobe import Cell
import src.index_check as index_check
from src.space_group_lib import *
from src.xscale_output import OutputParser
from src.run_command import *
from src.scale_utl import ScaleUtils
from src.ascii import ASCII
from src.dendro2highcharts import dendro2highcharts
import src.correlation as corc

from src.abstract import Abstract


logger = logging.getLogger('sxdm')


class Merging(Abstract):
    _command = "xscale_par > /dev/null"

    @staticmethod
    def getInDataScheme():
        return {
            "type": "object",
            "properties": {
                "xtallist": {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                "dirlist": {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                "experiment": {"type": "string"},
                "isXtal": {"type": "boolean"},
                "suffix": {"type": "string"},
                "reference": {"type": "string"},
                "resolution": {"type": "string"},
                "friedels_law": {"type": "string"},
                "unit_cell": {"type": "string"},
                "space_group": {"type": "string"},
                "nproc": {"type": "integer"},
                "isa_cutoff": {"type": "number"}
            }
        }

    @staticmethod
    def getOutDataScheme():
        return {
            "type": "object",
            "properties": {
                "hklpaths_found": {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                "isa_selected": {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                "hkl_cc_sorted": {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                "xtals_expected": {"type": "integer"},
                "xtals_found": {"type": "integer"},
                "xtals_after_idx_check": {"type": "integer"},
                "xtals_after_isa": {"type": "integer"},
                "xtals_after_cell": {"type": "integer"},
                "xtals_after_pCC": {"type": "integer"},
                "xtals_rejected": {"type": "integer"},
                "space_group": {"type": "string"},
                "unit-cell": {"type": "string"},
                "reference": {"type": "string"},
                "friedel": {"type": "string"},
                "anisotropicity": {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                'multiplicity': {
                    "type": "array",
                    "items": {
                        "type": "string"
                    }
                },
                'cell_array': {
                    "type": "array",
                    "items": {
                        "type": "number"
                    }
                },
                'cell_n_clusters': {"type": "integer"},
                'cc_n_clusters': {"type": "integer"},
                'clusters': {"type": "integer"},
                'hclust_matrix': {
                    "type": "array",
                    "items": {
                        "type": "number"
                    }
                },
                'cell_dendrogram': {}, #stores dendrogram dictionary, use 'ivl' key-values for labels when make plot
                'cc_dendrogram':{},
                'histogram': {} #cell_array converted to dictionary for highcharts

            }
        }

    def outdir_make(self):
        try:
            outname = 'adm_' + self.jshandle['experiment']
            self.setOutputDirectory()
            subadm = 'adm_%d' % self.results['xtals_expected']
            suffix = self.jshandle.get('suffix', None)
            if suffix is not None:
                subadm = subadm + "_" + suffix
            else:
                pass
            subadm = self.getOutputDirectory() / outname / subadm
            subadm.mkdir(parents=True, exist_ok=True)
            os.chdir(subadm)
            self.setOutputDirectory(path=subadm)
        except KeyError as e:
            logger.error(e)
            self.setFailure()
            return

    @staticmethod
    def sort_xtal(lists):
        sortnum = []
        sortlist = []
        if lists:
            for val in lists:
                num_sep = val.split('_')
                num = num_sep[1].split('.')
                sortnum.append(int(num[0]))

            for v in sorted(sortnum):
                element = 'xtal_'+str(v)+'.HKL'
                sortlist.append(element)
        return sortlist

    def find_HKLs(self):
        hklpaths = []
        xtallist = self.jshandle.get('xtallist', [])
        isXtal = self.jshandle.get('isXtal', False)
        if not isXtal:
            for folder in self.jshandle['dirlist']:
                dirs = pathlib.Path(folder)
                posix = list(dirs.glob('*/*'))
                xtallist += list(map(lambda x: str(x), posix))
        else:
            msg = "List of xtal folders directly provided, not searching\n"
            logger.info(msg)
            pass
        try:
            self.results['xtals_expected'] = len(xtallist)

            msg = "# of xtals expected: %d\n" %self.results['xtals_expected']
            logger.info('MSG: {}'.format(msg))

            for path in xtallist:
                filepath = os.path.join(path, "XDS_ASCII.HKL")
                if os.path.isfile(filepath):
                    hklpaths.append(filepath)

                elif os.path.isfile(path) and path.endswith('.HKL'):
                    hklpaths.append(path)

                else:
                    pass
                    #err = "couldn't find XDS_ASCII.HKLs %s" % filepath
                    #logger.info('Error:{}'.format(err))
        except KeyError as e:
            logger.error(e)
            self.setFailure()
        self.results['xtals_found'] = len(hklpaths)
        self.results['hklpaths_found'] = hklpaths
        return


    def create_file_links(self):

        self.results['filelinks'] = []
        try:
            for ii in range(len(self.results['hklpaths_found'])):
                try:
                    linkname = "xtal_"+str(ii)+".HKL"

                    os.symlink(self.results['hklpaths_found'][ii], linkname)
                    self.results['filelinks'].append(linkname)
                except OSError as e:
                    if e.errno == errno.EEXIST:
                        os.remove(linkname)
                        os.symlink(self.results['hklpaths_found'][ii], linkname)
                        self.results['filelinks'].append(linkname)

                    else:
                         logger.error(e)
                         self.setFailure()
        except (KeyError, IndexError) as err:
            logger.error(err)
            self.setFailure()
        return

    def indexing_(self):
        self.results['reference'] = self.jshandle.get('reference', self.results['hklpaths_found'][0])
        try:
            for ii in range(1, len(self.results['hklpaths_found'])):
                if not index_check.similar_symmetry(self.results['reference'], self.results['hklpaths_found'][ii]):
                    self.results['hklpaths_found'].pop(ii)
                    logger.info("wrong indexing\n")
                else:
                    pass

        except (IndexError, ValueError) as err:
            pass

        logger.info('MSG: # of cleaned xtals %s' %len(self.results['hklpaths_found']))
        self.results['xtals_after_idx_check'] = len(self.results['hklpaths_found'])
        return

    def create_inp(self, filelist, inData):

        friedel = self.jshandle.get("friedels_law", "TRUE")
        reso_cut = self.jshandle.get("resolution", "1.0")
        reference = inData.get('reference', self.results['hklpaths_found'][0])

        fh = open("XSCALE.INP",'w')
        fh.write("OUTPUT_FILE=XSCALE.HKL\n")
        fh.write("MAXIMUM_NUMBER_OF_PROCESSORS=%d\n" %self.jshandle.get('nproc', 12))
        fh.write("SAVE_CORRECTION_IMAGES=FALSE\n")
        fh.write("FRIEDEL'S_LAW=%s\n\n" %friedel)
        fh.write('REFERENCE_DATA_SET= %s\n' %reference)
        try:
            fh.write("SPACE_GROUP_NUMBER=%s\n" %inData['space_group'])
            fh.write("UNIT_CELL_CONSTANTS=%s\n" %inData['unit_cell'])
        except KeyError as err:
            logger.info("didnot add spg/unit_cell")

        for f in filelist:
            fh.write("INPUT_FILE=%s\n" %f)
            fh.write("SNRC=0.0\n")
            fh.write("INCLUDE_RESOLUTION_RANGE= 50 %f\n" %float(reso_cut))

        fh.close()
        return

    def Isa_select(self, config):
        # method to perform ISa based selection of 'good' HKL files for next round of XSCALEing
        isa_threshold = self.jshandle.get('isa_cutoff', 3.0)
        self.results['isa_selected'] = []

        if os.path.isfile("XSCALE.LP"):
            fh = open("XSCALE.LP", 'r')
            all_lines = fh.readlines()
            fh.close()
            try:
                start = all_lines.index('     a        b          ISa    ISa0   INPUT DATA SET\n')
                end = all_lines.index(' (ASSUMING A PROTEIN WITH 50% SOLVENT)\n')
                # select the table with ISa values from the file..
                isa_list = all_lines[start+1:end-3]
            # except (ValueError, IndexError) as err:
            #     error = "check if XSCALE ran properly or XSCALE.INP screwed up \n"
            #     logger.error(error)


                for lines in isa_list:
                    line = lines.split()
                    try:
                        if float(line[2]) > float(isa_threshold):
                            self.results['isa_selected'].append(line[4])
                    except IndexError:
                        pass

                self.results['isa_selected'] = Merging.sort_xtal(self.results['isa_selected'])
                self.create_inp(self.results['isa_selected'], config)

                msg = "%d xtals selected by ISa-cut out of %d xtals\n" %(len(self.results['isa_selected']),
                                                                         len(self.results['hklpaths_found']))
                self.results['xtals_after_isa'] = len(self.results['isa_selected'])

                logger.info('MSG:{}'.format(msg))
                try:
                    # user = self.jshandle.get('user', " ")
                    run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], Merging._command, self.getLogFileName())
                except KeyError as e:
                    sub.call(Merging._command, shell=True)

            except Exception as e:
                logger.error(e)
                self.setFailure()
        else:
            err = "XSCALE.LP file does not exist"
            logger.error(err)
            self.setFailure()
        return

    def xscale_for_sad(self):
        self.find_HKLs()

        self.outdir_make()
        config = dict()
        if len(self.results['hklpaths_found']) > 0:
            if len(self.results['hklpaths_found']) < 1000:
                self.indexing_()
            else:
                logger.info('Too many hkls, so skipping the indexing check')
                pass
            self.create_file_links()
            indict = {"listofHKLfiles": self.results['hklpaths_found'],
                      "fom":'bfac'}
            sc = ScaleUtils(indict)
            sc.ref_choice(indict)

            logging.info('lowest-Bfactor file: %s' %sc.results['reference'])
            indata_ascii = {"xds_ascii": sc.results['reference']}
            ref_for_cell_sg = ASCII(indata_ascii)
            ref_for_cell_sg.get_data(indata_ascii)
            self.results['space_group'] = space_group[ref_for_cell_sg.results['spg']][0]
            self.results['unit-cell'] = ref_for_cell_sg.results['unit_cell']
            self.results['friedel'] = ref_for_cell_sg.results['anom']
            config['reference'] = sc.results['reference']
            self.create_file_links()
            self.create_inp(self.results['filelinks'], config)
        else:
            # self.create_inp(self.results['hklpaths_found'], config)
            logger.error("No hkl file found")

        msg = "xscale-ing of native SAD data\n"
        logger.info('MSG:{}'.format(msg))
        try:
            run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], Merging._command, self.getLogFileName())
        except KeyError:
            sub.call(Merging._command, shell=True)

        try:
            indict = {"LPfile": "XSCALE.LP"}
            xscale_parse = OutputParser(indict)
            xscale_parse.parse_xscale_output(indict)
            logger.info('stat_dict:{}'.format(xscale_parse.results))
            self.results['nSAD_xscale_stats'] = xscale_parse.results
            shutil.copyfile('XSCALE.HKL', 'nSAD.HKL')
            shutil.copyfile('XSCALE.LP', 'nSAD.LP')

        except Exception as err:
            logger.error(err)
            self.setFailure()

        self.Isa_select(config)
        shutil.copyfile('XSCALE.LP', 'ISa_Select.LP')
        shutil.copyfile('XSCALE.HKL', 'ISa_Select.HKL')
        indict = {"LPfile": "ISa_Select.LP"}
        xscale_parse = OutputParser(indict)
        xscale_parse.parse_xscale_output(indict)
        self.results['ISa_selection'] = xscale_parse.results
        return

    def xscale_for_sx(self):
        self.find_HKLs()
        config = dict()
        try:
            self.outdir_make()

        except OSError:
            err = "adm folder either not there or not accessible\n"
            logger.error(err)
            self.setFailure()
            return

        try:
            if len(self.results['hklpaths_found']) < 1000:
                self.indexing_()
            else:
                msg = "Too many hkls, so skipping indexing check"
                logger.info('MSG:{}'.format(msg))
                pass

            indict = {"listofHKLfiles": self.results['hklpaths_found'],
                      "fom":'bfac'}
            sc = ScaleUtils(indict)
            sc.ref_choice(indict)

            logging.info('lowest-Rmeas file: %s' %sc.results['reference'])
            indata_ascii = {"xds_ascii": sc.results['reference']}
            ref_for_cell_sg = ASCII(indata_ascii)
            ref_for_cell_sg.get_data(indata_ascii)

            self.results['space_group'] = space_group[ref_for_cell_sg.results['spg']][0]
            self.results['unit-cell'] = ref_for_cell_sg.results['unit_cell']
            self.results['friedel'] = ref_for_cell_sg.results['anom']
            config['reference'] = sc.results['reference']
            self.create_file_links()
            self.create_inp(self.results['filelinks'], config)

            msg = "Running 1st round of xscale-ing with Rmeas based ranking\n"
            logger.info('MSG:{}'.format(msg))
            try:
                run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], Merging._command, self.getOutputDirectory())
            except KeyError:
                sub.call(Merging._command, shell=True)
            try:
                indict = {"LPfile": "XSCALE.LP"}
                xscale_parse = OutputParser(indict)

                xscale_parse.parse_xscale_output(indict)
                logger.info('stat_dict:{}'.format(xscale_parse.results))
                self.results['no_selection'] = xscale_parse.results
                shutil.copyfile("XSCALE.INP", "noSelect.INP")
                shutil.copyfile("XSCALE.LP", "noSelect.LP")
                shutil.copyfile("XSCALE.HKL", "noSelect.HKL")

            except OSError:
                err = "xscaling may not have run\n"
                logger.error(err)
                self.setFailure()
                return


            msg = "running xscale after ISa selection\n"
            logger.info(msg)
            self.Isa_select(config)

            try:
                indict = {"LPfile": "XSCALE.LP"}
                xscale_parse = OutputParser(indict)
                xscale_parse.parse_xscale_output(indict)

                logger.info('stat_dict:{}'.format(xscale_parse.results))
                self.results['ISa_selection'] = xscale_parse.results
                shutil.copyfile("XSCALE.INP", "ISa_Select.INP")
                shutil.copyfile("XSCALE.LP", "ISa_Select.LP")
                shutil.copyfile("XSCALE.HKL", "ISa_Select.HKL")

            except OSError:
                err = "xscaling after ISa selection may not work\n"
                logger.info('OSError:{}'.format(err))
                self.setFailure()

                return
            celldict = {"listofHKLfiles": self.results['filelinks']}
            celler = Cell(celldict)
            celler.setter(celldict)
            if len(celler.results['hklList']) > 200:
                celler.cell_analysis(celldict)
                self.results['cell_array'] = celler.results['cell_ar'].tolist()
                self.results['histogram'] = celler.dict_for_histogram()
            else:
                celler.clustering(celldict)
                self.results['cell_array'] = celler.results['cell_ar_best_cluster'].tolist()
                self.results['dendro_labels'] = celler.results['data_points']



            mode_cells = str(celler.results['a_mode']) + "  " + str(celler.results['b_mode']) + '  ' + str(celler.results['c_mode']) + \
            "  " + str(celler.results['al_mode']) + "  " + str(celler.results['be_mode']) + "   " + str(celler.results['ga_mode'])

            # self.results['unit-cell'] = mode_cells
            # convert dendro dictionary for easy plottable dictionary in adp tracker
            try:
                self.results['cell_dendrogram'] = dendro2highcharts(celler.results['dendo'])
                # self.scale_results['hclust_matrix'] = celler.hclust_matrix
                self.results['cell_n_clusters'] = celler.results['n_clusters_cell']
            except Exception as e:
                logger.info('skipping-dendro-cell:{}'.format(e))

            try:
                config['unit_cell'] = mode_cells
                config['space_group'] = ref_for_cell_sg.results['spg']

                self.create_inp(celler.results['cell_select'], config)
                msg = "xscaling after cell-analysis and rejecting outliers\n"
                msg += "%d xtals got selected after cell analysis out of %d xtals" \
                      %(len(celler.results['cell_select']), len(self.results['hklpaths_found']))
                logger.info('MSG:{}'.format(msg))

                try:
                    run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], Merging._command, self.getLogFileName())
                except (OSError, TypeError, KeyError) as e:
                    sub.call(Merging._command, shell=True)

                self.results['xtals_after_cell'] = len(celler.results['cell_select'])

            except Exception as e:
                logger.info('Error:{}'.format(e))
                self.setFailure()

            try:
                indict = {"LPfile": "XSCALE.LP"}
                xscale_parse = OutputParser(indict)
                xscale_parse.parse_xscale_output(indict)

                logger.info('stat_dict:{}'.format(xscale_parse.results))
                self.results['cell_selection'] = xscale_parse.results
                shutil.copyfile("XSCALE.INP", "Cell_Select.INP")
                shutil.copyfile("XSCALE.LP", "Cell_Select.LP")
                shutil.copyfile("XSCALE.HKL", "Cell_Select.HKL")

            except OSError:
                err = "xscaling after Cell selection may not work\n"
                logger.info('OSError:{}'.format(err))

            ccDict = {"xds_ascii": 'noSelect.HKL'}
            CC = corc.CCestimator(ccDict)
            logger.info('ASCII loaded')
            CC.cc_select(fom='pcc')
            try:
                self.results['cc_dendrogram'] = dendro2highcharts(CC.results['cc_dendo'])
                self.results['cc_n_clusters'] = CC.results['n_clusters_cc']
                self.results['hkl_cc_sorted'] = CC.results['cc_cluster_list']
                self.results['xtals_after_pCC'] = len(self.results['hkl_cc_sorted'])
                msg = "pair-correlation sorting over Isa_select: %d hkls" %self.results['xtals_after_pCC']
                logger.info('MSG:{}'.format(msg))
                self.create_inp(self.results['hkl_cc_sorted'], config)
            except Exception as err:
                logger.info('cc-dendro-empty:{}'.format(err))
            try:
                run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], Merging._command, self.getLogFileName())
            except (OSError, TypeError, Exception) as e :
                sub.call(Merging._command, shell=True)



            try:
                indict = {"LPfile": "XSCALE.LP"}
                xscale_parse = OutputParser(indict)
                xscale_parse.parse_xscale_output(indict)

                logger.info('stat_dict:{}'.format(xscale_parse.results))
                self.results['pCC_selection'] = xscale_parse.results
                shutil.copyfile("XSCALE.INP", "pCC_Select.INP")
                shutil.copyfile("XSCALE.LP", "pCC_Select.LP")
                shutil.copyfile("XSCALE.HKL", "pCC_Select.HKL")

            except OSError:
                err = "xscaling after pair-correlation selection may not work\n"
                logger.info('OSError:{}'.format(err))
                self.setFailure()

        except Exception as e:
            logger.info('Error: {}'.format(e))
            self.setFailure()


        return

    def isocluster(self, xscalefile):

        if os.path.isfile(xscalefile):
            comd = "xscale_isocluster -dim 3 -dmax %s %s" %(self.jshandle['resolution'], xscalefile)
            msg = "iso-clustering based on Correlation\n"
            logger.info('MSG:{}'.format(msg))
            try:
                run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], comd, 'isocluster.log')
            except KeyError:
                comd1 = "xscale_isocluster %s > isocluster.log" %xscalefile
                sub.call(comd1, shell=True)

        else:
            self.setFailure()
            return

        if os.path.isfile(self.getOutputDirectory() / "isocluster.log"):
            fkey = "best number of clusters"

            fh = open("isocluster.log", "r")
            _all = fh.readlines(); fh.close()
            for lines in _all:
                if "Error" in lines:
                    self.setFailure()
                elif fkey in lines:
                    line = lines.split(':')
                    val = line[-1].strip('\n')
                    self.results["clusters"] = val.strip(' ')
                else:
                    pass
        if os.path.isfile(str(self.getOutputDirectory(),"XSCALE.1.INP")):
            shutil.copyfile("XSCALE.1.INP", "XSCALE.INP")
            msg = "xscale-ing over 1st cluster only..\n"
            logger.info('MSG:{}'.format(msg))
            try:
                run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], Merging._command, self.getLogFileName())
            except KeyError:
                sub.call(Merging._command, shell=True)
            try:
                indict = {"LPfile": "XSCALE.LP"}
                xscale_parse = OutputParser(indict)
                xscale_parse.parse_xscale_output(indict)

                logger.info('stat_dict:{}'.format(xscale_parse.results))
                shutil.copyfile("XSCALE.INP", "cluster1.INP")
                shutil.copyfile("XSCALE.LP", "cluster1.LP")
                shutil.copyfile("XSCALE.HKL", "cluster1.HKL")
                self.results['iso-cluster'] = xscale_parse.results

            except OSError:
                err = "xscaling after iso-clustering may not work\n"
                logger.info('OSError:{}'.format(err))
                self.setFailure()
                return
        return

    def aniso_check(self):
        point_cmd = "pointless -xdsin noSelect.HKL hklout noSelect_point.mtz"
        try:
            run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], point_cmd, 'pointless.log')
        except KeyError:
            point_cmd1 = "pointless -xdsin noSelect.HKL hklout noSelect_point.mtz > pointless.log"
            sub.call(point_cmd1, shell=True)
        time.sleep(2)

        if os.path.isfile(str(self.getOutputDirectory() / "noSelect_point.mtz")):

            fh = open("onlymerge.inp", 'w')
            fh.write("ONLYMERGE\n")
            try:
                from iotbx import reflection_file_reader
                mtzfile = reflection_file_reader.any_reflection_file("noSelect_point.mtz")
                mtz_content = mtzfile.file_content()
                col_labels = mtz_content.column_labels()
                if "BATCH" in col_labels:
                    n_batches = mtz_content.n_batches()
                    fh.write("run 1 batch %d to %d\n" %(1, n_batches))
                else:
                    pass
            except (ImportError, RuntimeError) as err:
                logger.info('iotbx-import-error:{}'.format(err))
                fh.write("run 1 batch %d to %d\n" %(1, 50))
            fh.close()
            aim_cmd = "aimless hklin noSelect_point.mtz hklout noSelect_aim.mtz < onlymerge.inp"
            try:
                run_command("sxdm", self.getOutputDirectory(), self.jshandle['user'], aim_cmd, 'aimless.log')
            except (OSError, TypeError, Exception) as e:
                aim_cmd1 = aim_cmd + '| tee aimless.log'
                sub.call(aim_cmd1, shell=True)
        else:
            err = "Could be pointless did not run\n"
            logger.info('aimless-error:{}'.format(err))
            self.setFailure()
            return
        try:
            fh = open("aimless.log", 'r')
            _all = fh.readlines()
            fh.close()
        except (OSError, IOError):
            err = "aimless.log file doesn't exist"
            logger.info('OSError:{}'.format(err))
            self.setFailure()
            return

        keyword = "Estimated maximum resolution limits"

        for lines in _all:

            if "Error" in lines:
                self.setFailure()
            elif keyword in lines:
                line = lines.split(',')
                try:
                    a = line[1].split(':')[1]
                    b = line[2].split(':')[1]
                    c = line[3].split(':')[1].strip('\n')
                    aniso_string = "a* =%s, b*=%s, c*=%s" %(a,b,c)
                    self.results['anisotropicity'] = aniso_string
                    logger.info('Anisotropy:{}'.format(aniso_string))

                except IndexError:
                    hk_plane = line[1].split(':')[1]
                    l_axis = line[2].split(':')[1].strip('\n')
                    aniso_string = "hk-plane = %s, l-axis = %s" %(hk_plane, l_axis)
                    self.results['anisotropicity'] = aniso_string

            else:
                pass
        return

    def create_mtzs(self):
        lst_of_hkls = ['noSelect.HKL', 'ISa_Select.HKL', 'Cell_Select.HKL', 'pCC_Select.HKL']
        indata_ascii = {"mtz_format": "CCP4_I+F",
                        "resolution": float(self.jshandle['resolution']),
                        }
        for hklfile in lst_of_hkls:
            if os.path.isfile(hklfile):
                indata_ascii['xds_ascii'] = hklfile
                hkl = ASCII(indata_ascii)
                hkl.get_data(indata_ascii)
                hkl.xdsconv(indata_ascii)
            else:
                logger.info('mtz-create-info:{}'.format('could not create mtz; file may not exist'))
                pass
        return

    def run_(self):

        try:
            # jsonschema.validate(instance=Merging.jshandle, schema=Merging.getInDataScheme())
            if self.jshandle['experiment'] == 'native-sad':
                self.xscale_for_sad()
            elif self.jshandle['experiment'] == 'serial-xtal':
                self.xscale_for_sx()
                indict = {'xds_ascii': 'noSelect.HKL'}
                xscale = ASCII(indict)
                xscale.get_data(indict)
                xscale.multiplicity_check()
                self.results['multiplicity'] = xscale.results['multiplicity']

                self.aniso_check()
                self.create_mtzs()

            elif self.jshandle['experiment'] == 'inverse-beam' or self.jshandle['experiment'] == 'interleave-and-inverse-first':
                pass
            else:
                logger.info('Error: experiment type not supported %s' %self.jshandle['experiment'])
                self.setFailure()
        except Exception as e:
            logger.info('Error:{}'.format(e))
            self.setFailure()
        return


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
    datefmt='%y-%m-%d %H:%M',
    filename='merge.log',
    filemode='a')

    import json
    fh = open(sys.argv[1], 'r')
    inData = json.load(fh)
    fh.close()

    mm = Merging(inData)
    # val = jsonschema.Draft3Validator(Merging.getInDataScheme()).is_valid(inData)
    # print(val)

    mm.run_()
    if mm.is_success():
        mm.writeOutputData(mm.results)
