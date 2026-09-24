
__author__ = ["S. Basu"]
__license__ = "M.I.T"
__date__ = "28/06/2017"
__refactordate__ = "07/05/2021"

import os, sys
import logging
from src.abstract import Abstract

logger = logging.getLogger('sxdm')

class OutputParser(Abstract):

    @property
    def getInDataScheme(self):
        return {
            "type": "object",
            "properties": {
                "LPfile": {"type": "string"},
                "CORRECT_file": {"type": "string"}
            }
        }

    def parse_xds_stats(self, inData):
        self.results['xds_stat'] = []
        count = 1
        skip_string = ['RESOLUTION', 'LIMIT']
        search_string = ' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n'
        try:
            filename = inData['CORRECT_file']
            if not os.path.isfile(filename):
                e = IOError("%s file does not exist\n" %filename)
                logger.error(e)
                self.setFailure()
                return
        except KeyError:
            self.setFailure()
            logger.error("XSCALE LP file is not provided")
            return

        fh = open(filename, 'r')
        _all = fh.readlines()
        fh.close()
        try:
           start_idx = [idx for idx, search in enumerate(_all) if search == search_string]
           useful = _all[start_idx[-1]+1:start_idx[-1]+13] # get only last chunk of statistics tables

           for lines in useful:
              if any(k in lines for k in skip_string):
                 pass
              else:
                 line = lines.split()
                 if len(line) > 0:
                    each_row = dict()
                    each_row['order'] = count
                    each_row['resolution_limit'] = float(line[0])
                    each_row['observed_reflections'] = int(line[1])
                    each_row['unique_reflections'] = int(line[2])
                    each_row['possible_reflections'] = int(line[3])
                    each_row['completeness'] = float(line[4].strip('%'))
                    each_row['r_factor_observed'] = float(line[5].strip('%'))
                    each_row['r_factor_expected'] = float(line[6].strip('%'))
                    each_row['compared_reflections'] = int(line[7])
                    each_row['Isigma'] = float(line[8])
                    each_row['rmeas'] = float(line[9].strip('%'))
                    each_row['sigAno'] = float(line[12])
                    each_row['anom_compared'] = float(line[13])
                    if line[10][-1] == '*':
                        each_row['cc_half'] = float(line[10].strip('*'))
                    else:
                        each_row['cc_half'] = float(line[10])
                    if line[11][-1] == '*':
                        each_row['anomalous_correlation'] = float(line[11].strip('*'))
                    else:
                        each_row['anomalous_correlation'] = float(line[11])
                    self.results['xds_stat'].append(each_row)
                    count += 1
                 else:
                     pass
        except Exception as e:
             print(e)
             pass


    def parse_xscale_output(self, inData):
        self.results['stat'] = []
        count = 1
        try:
            filename = inData['LPfile']
            if not os.path.isfile(filename):
                e = IOError("%s file does not exist\n" %filename)
                logger.error(e)
                self.setFailure()
                return
        except KeyError:
            self.setFailure()
            logger.error("XSCALE LP file is not provided")
            return

        skip_line = ['RESOLUTION', 'LIMIT', 'total']


        fh = open(filename, 'r')
        _all = fh.readlines()
        fh.close()
        try:
            start = _all.index(' SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION\n')
            end = _all.index(' ========== STATISTICS OF INPUT DATA SET ==========\n')
            useful = _all[start+1:end-1]

        except (TypeError, ValueError, IndexError) as err:
            logger.error('there was no XDS_ASCII files, so nothing to XSCALE. Check input/output directory naming')
            logger.error(err)
            self.setFailure()
            return

        for lines in useful:
            if any(k in lines for k in skip_line):
                pass
            else:
                logger.info(lines[:-1])
                line = lines.split()
                if len(line) == 14:
                    each_row = dict()
                    each_row['order'] = count
                    each_row['resolution_limit'] = float(line[0])
                    each_row['observed_reflections'] = int(line[1])
                    each_row['unique_reflections'] = int(line[2])
                    each_row['possible_reflections'] = int(line[3])
                    each_row['completeness'] = float(line[4].strip('%'))
                    each_row['r_factor_observed'] = float(line[5].strip('%'))
                    each_row['r_factor_expected'] = float(line[6].strip('%'))
                    each_row['compared_reflections'] = int(line[7])
                    each_row['Isigma'] = float(line[8])
                    each_row['rmeas'] = float(line[9].strip('%'))
                    each_row['sigAno'] = float(line[12])
                    each_row['anom_compared'] = float(line[13])
                    if line[10][-1] == '*':
                        each_row['cc_half'] = float(line[10].strip('*'))
                    else:
                        each_row['cc_half'] = float(line[10])
                    if line[11][-1] == '*':
                        each_row['anomalous_correlation'] = float(line[11].strip('*'))
                    else:
                        each_row['anomalous_correlation'] = float(line[11])
                    self.results['stat'].append(each_row)
                    count += 1
                elif 0 < len(line) < 14:
                    msg = "float formatting error in xscale table"
                    logger.info('xscale-output-error:{}'.format(msg))
                    each_row = dict()
                    each_row['order'] = count
                    each_row['resolution_limit'] = float(line[0])
                    each_row['observed_reflections'] = int(line[1])
                    each_row['unique_reflections'] = int(line[2])
                    each_row['possible_reflections'] = int(line[3])
                    each_row['completeness'] = float(line[4].strip('%'))
                    each_row['r_factor_observed'] = 'NaN'
                    each_row['r_factor_expected'] = 'NaN'
                    each_row['compared_reflections'] = 'NaN'
                    each_row['Isigma'] = 'NaN'
                    each_row['rmeas'] = 'NaN'
                    each_row['sigAno'] = 'NaN'
                    each_row['anom_compared'] = 'NaN'
                    each_row['cc_half'] = 'NaN'
                    each_row['anomalous_correlation'] = 'NaN'
                    self.results['stat'].append(each_row)
                    count += 1
                else:
                    pass

        return

    @staticmethod
    def mean_rmeas_calc(stat):
        rmeas_lst = []
        mean = 0.0
        if len(stat) > 0:
            for each_row in stat:
                rmeas_lst.append(each_row['rmeas']/100.0)
            mean = (sum(rmeas_lst)/len(rmeas_lst))*100.0

        else:
            err = "CORRECT.LP file doesn't have stats"
            logger.error(err)
        return mean

    @staticmethod
    def display_table(inData):

        """
        :param inData: basically pass self.results or output dictionary of parse_xscale_output function
        :return: void
        """

        ofh = open('table.txt', 'w')
        ofh.write('Resolution   observed    unique   possible  completeness    Robs      Rexpect   compared  I/sig      Rmeas     CC1/2   CCano  anom_comp  sig_ano\n')
        try:
            for row in inData['stat']:
                ofh.write("%6.2f %12d %10d %10d %12.1f %11.1f %11.1f %9d %8.1f %10.1f %9.1f %7.1f %8d %9.3f\n"\
                    %(row['resolution_limit'], row['observed_reflections'], \
                        row['unique_reflections'], row['possible_reflections'], \
                        row['completeness'], row['r_factor_observed'], \
                        row['r_factor_expected'], row['compared_reflections'], \
                        row['Isigma'], row['rmeas'], row['cc_half'], \
                        row['anomalous_correlation'], row['anom_compared'], row['sigAno']))
            ofh.close()
        except KeyError:
            logger.error("run parse_xscale_output to generate self.results; then call this function with correct argument \n")

    @staticmethod
    def print_table(stat):
        if stat:
            print('Resolution   observed    unique   possible  completeness      Robs      Rexpect   compared  I/sig      Rmeas     CC1/2   CCano  anom_comp  sig_ano')
            for row in stat:
                print("%6.2f %12d %10d %10d %12.1f %11.1f %11.1f %9d %8.1f %10.1f %9.1f %7.1f %8d %9.3f"\
                    %(row['resolution_limit'], row['observed_reflections'], \
                    row['unique_reflections'], row['possible_reflections'], \
                    row['completeness'], row['r_factor_observed'], \
                    row['r_factor_expected'], row['compared_reflections'], \
                    row['Isigma'], row['rmeas'], row['cc_half'], \
                    row['anomalous_correlation'], row['anom_compared'], row['sigAno']))
        else:
            print("no stats to display, xscale did not run\n")

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%a, %d %b %Y %H:%M:%S',
    filename='stats.log',
    filemode='w')

    indict = {"CORRECT_file": '/nfs/ssx/shbasu/MEmmery/processed/proc_20201029/AD029A/atx_01/xtal_0/CORRECT.LP'}
    xscale = OutputParser(indict)
    xscale.parse_xds_stats(indict)
    print(xscale.results['xds_stat'])
    xscale.print_table(xscale.results['xds_stat'])

