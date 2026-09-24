
__author__ = ["S. Basu"]
__license__ = "M.I.T"
__date__ = "28/07/2017"
__refactordate__ = "10/05/2021"

import os, sys,glob
import logging
from src.abstract import Abstract
from src.xscale_output import OutputParser

logger = logging.getLogger('sxdm')

class ScaleUtils(Abstract):

    def find_corrects(self, inData):
        self.results['listofCORRECTfiles'] = []
        try:
            for fname in inData['listofHKLfiles']:
                folder = os.path.dirname(fname)
                path = os.path.join(folder, 'CORRECT.LP')

                if os.path.isfile(path):
                    self.results['listofCORRECTfiles'].append(path)

                else:
                    logger.info('CORRECT.LP could not be found in %s' %folder)
                    self.setFailure()
        except KeyError:
            self.setFailure()
        return

    def check_bfactor(self, inData):
        self.find_corrects(inData)
        bfac_dicts = {}

        if len(self.results['listofCORRECTfiles']) == 0:
            err = 'ValueError: no CORRECT.LP found'
            logger.info('ValueError: {}'.format(err))

        else:
            for fname, cor_name in zip(inData['listofHKLfiles'], self.results['listofCORRECTfiles']):
                fh = open(cor_name, 'r')
                _all = fh.readlines()
                fh.close()
                xasci = fname
                for lines in _all:
                    if "WILSON LINE" in lines:
                        line = lines.split()
                        try:
                            bfac_dicts[xasci] = float(line[9])
                        except Exception:
                            logger.info('B-factor might be negative, not considered')
                    else:
                        pass

            self.results['bfac_sorted_hkls'] = sorted(bfac_dicts.items(), key=lambda x : x[1])
        return

    def rank_rmeas(self, inData):
        self.find_corrects(inData)
        rmeas_dict = {}

        if len(self.results['listofCORRECTfiles']) == 0:
            err = 'ValueError: no CORRECT.LP found'
            logger.info('ValueError: {}'.format(err))
        else:
            for fname, cor_name in zip(inData['listofHKLfiles'], self.results['listofCORRECTfiles']):
                indict = {'CORRECT_file': cor_name}
                correct_parse = OutputParser(indict)
                correct_parse.parse_xds_stats(indict)
                mean_rmeas = correct_parse.mean_rmeas_calc(correct_parse.results['xds_stat'])
                rmeas_dict[fname] = mean_rmeas
            self.results['rmeas_sorted_hkls'] = sorted(rmeas_dict.items(), key=lambda x:x[1])
        return


    def ref_choice(self, inData):
        reference = None
        if inData['fom'] == 'bfac':
            self.check_bfactor(inData)
            try:
                reference = self.results['bfac_sorted_hkls'][0][0]

            except (IndexError, ValueError):
                err = 'bfactor selection may not work'
                logger.error(err)
                self.setFailure()


        elif inData['fom'] == 'rmeas':
            self.rank_rmeas(inData)
            try:
                reference = self.results['rmeas_sorted_hkls'][0][0]
            except (IndexError, ValueError):
                err = 'Rmeas based referenceing may not have worked'
                logger.error(err)
                self.setFailure()
        else:
            pass
        self.results['reference'] = reference
        print(self.results)
        return

    def Bfact_sorter(self, inData):
        bfac_sorted_hkls = []
        self.check_bfactor(inData)
        if len(self.results['bfac_sorted_hkls']) > 0:
            for i in range(len(self.results['bfac_sorted_hkls'])):
                bfac_sorted_hkls.append(self.results['bfac_sorted_hkls'][i][0])
        else:
            err = "Rmeas based sorting did not work, check"
            logger.error(err)
            self.setFailure()

        self.results['bfact_sorted_hkls'] = bfac_sorted_hkls
        return

    def rmeas_sorter(self, inData):
        rmeas_sorted_hkls = []
        self.rank_rmeas(inData)
        if len(self.results['rmeas_sorted_hkls']) > 0:
            for i in range(len(self.results['rmeas_sorted_hkls'])):
                rmeas_sorted_hkls.append(self.results['rmeas_sorted_hkls'][i][0])
        else:
            err = "Rmeas based sorting did not work, check"
            logger.error(err)
            self.setFailure()
        self.results['rmeas_sorted_hkls'] = rmeas_sorted_hkls
        return

def main():
    hklpaths = glob.glob(os.path.join(sys.argv[1], 'XDS_ASCII.HKL'))
    inData = dict()
    inData['listofHKLfiles'] = hklpaths
    sc = ScaleUtils(inData)
    sc.rmeas_sorter(inData)
    print(sc.results['rmeas_sorted_hkls'])

if __name__ == '__main__':
    main()
