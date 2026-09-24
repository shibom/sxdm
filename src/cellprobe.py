

__author__ = ["S. Basu"]
__license__ = "M.I.T"
__date__ = "26/06/2017"
__refactordate__ = "25/05/2021"

import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import scipy.cluster.hierarchy as sch
import logging
from src.abstract import Abstract

logger = logging.getLogger('sxdm')

class Cell(Abstract):

    @property
    def getInDataScheme(self):
        return {
            "type": "object",
            "properties": {
                "listofHKLfiles": {
                    "type": "array",
                    "items": {"type": "string"}
                },
                "XSCALEfile": {"type": "string"}
            }
        }

    def setter(self, inData):

        try:
            self.results['hklList'] = inData["listofHKLfiles"]
        except KeyError:
            if inData['XSCALEfile'].endswith('.LP'):
                self.results['hklList'] = Cell.get_filenames(inData['XSCALEfile'])
            else:
                err = ValueError("Reprocess with a lower ISa_cutoff, no XDS_ASCII.HKL file has enough ISa\n")
                logger.error(err)

        try:
            self.results['cell_ar'] = np.empty((len(self.results['hklList']), 6))
            self.results['cell_vector'] = np.empty((len(self.results['hklList']), 3))
            self.results['cell_select'] = []
            self.results['dendo'] = {}
            self.results['hclust_matrix'] = np.empty((len(self.results['hklList'])-1, 4))
        except Exception as err:
            logger.error(err)
            self.setFailure()
        return

    @staticmethod
    def get_filenames(filename):
        fLst1 = []; fLst2 = []
        if os.path.isfile(filename):
            fh = open(filename, 'r')
            all_lines = fh.readlines()
            fh.close()

            for lines in all_lines:
                if "INPUT_FILE" in lines:
                    line = lines.split('=')
                    fLst1.append(line[1])
                #else:
                #	pass
            for name in fLst1:
                tmp = name.strip('\n')
                fLst2.append(tmp)
        else:
            err = OSError("%s does not exist" %filename)
            logger.err(err)

        return fLst2

    @staticmethod
    def get_cells(filename):
        unit_cell = dict()
        try:
            fh = open(filename, 'r')
            all_lines = fh.readlines()
            fh.close()
            for lines in all_lines:
                if "!UNIT_CELL_CONSTANTS" in lines:
                    line = lines.split()
                    unit_cell['a'] = float(line[1])
                    unit_cell['b'] = float(line[2])
                    unit_cell['c'] = float(line[3])
                    unit_cell['al'] = float(line[4])
                    unit_cell['be'] = float(line[5])
                    unit_cell['ga'] = float(line[6])
                elif "!SPACE_GROUP_NUMBER" in lines:
                    line = lines.split()
                    unit_cell['spg'] = int(line[1])
        except Exception as e:
            logger.error(e)
        return unit_cell

    @staticmethod
    def calc_vector(cell_dict):
        Vab = np.sqrt(cell_dict['a']**2 + cell_dict['b']**2 - 2*cell_dict['a']*cell_dict['b']*np.cos(np.pi-np.deg2rad(cell_dict['ga'])))
        Vbc = np.sqrt(cell_dict['b']**2 + cell_dict['c']**2 - 2*cell_dict['b']*cell_dict['c']*np.cos(np.pi-np.deg2rad(cell_dict['al'])))
        Vca = np.sqrt(cell_dict['c']**2 + cell_dict['a']**2 - 2*cell_dict['c']*cell_dict['a']*np.cos(np.pi-np.deg2rad(cell_dict['be'])))

        return Vab, Vbc, Vca

    @staticmethod
    def reject_outlier(array):
        ohne_outlier = array[abs(array - np.mean(array)) <= 1.5*np.std(array)]
        ohne_indx = np.where(abs(array - np.mean(array)) <= 1.5*np.std(array))
        ohne_indx = ohne_indx[0].tolist()

        return ohne_outlier, ohne_indx

    def lcv_(self):

        for i in range(len(self.results['hklList'])):
            unit_cell = Cell.get_cells(self.results['hklList'][i])
            Vab, Vbc, Vca = Cell.calc_vector(unit_cell)

            self.results['cell_vector'][i,0] = Vab
            self.results['cell_vector'][i,1] = Vbc
            self.results['cell_vector'][i,2] = Vca
        return

    @staticmethod
    def cluster_indices(cluster_assignments):
        n = cluster_assignments.max()
        indices = []
        for cluster_number in range(1, n + 1):
            indices.append(np.where(cluster_assignments == cluster_number)[0])

        return indices

    def clustering(self, inData):
        self.setter(inData)
        self.lcv_()
        self.results['data_points'] = []
        for item in self.results['hklList']:
            tmp = os.path.basename(item)
            self.results['data_points'].append(tmp.strip('.HKL'))


        #hierarchical clustering using resultant vectors from unit cells
        Y = sch.linkage(self.results['cell_vector'], method='ward')
        self.results['hclust_matrix'] = Y
        assign = sch.fcluster(Y, 3, 'distance')
        MSG = "# clusters: %d" %assign.max()
        self.results['n_clusters_cell'] = int(assign.max())
        logger.info('MSG: {}'.format(MSG))
        idx = Cell.cluster_indices(assign)
        max_size = []
        for k, i in enumerate(idx):
            max_size.append(len(i))
            logger.info("cluster: %d; size: %d" %((k+1),len(i)))
        msg = "most populated cluster: %d" %(max_size.index(max(max_size))+1)
        logger.info('MSG: {}'.format(msg))
        best_cluster_name = max_size.index(max(max_size))
        self.results['cell_ar_best_cluster'] = np.empty((len(idx[best_cluster_name]), 6))
        for i in range(len(idx[best_cluster_name])):
            item = idx[best_cluster_name][i]
            cell_selected_file = self.results['hklList'][item]
            self.results['cell_select'].append(cell_selected_file)
            unit_cell = Cell.get_cells(cell_selected_file)
            self.results['cell_ar_best_cluster'][i,0] = unit_cell['a']
            self.results['cell_ar_best_cluster'][i,1] = unit_cell['b']
            self.results['cell_ar_best_cluster'][i,2] = unit_cell['c']
            self.results['cell_ar_best_cluster'][i,3] = unit_cell['al']
            self.results['cell_ar_best_cluster'][i,4] = unit_cell['be']
            self.results['cell_ar_best_cluster'][i,5] = unit_cell['ga']
        #calculate median values for unit cell parameters from most populated cluster..
        self.results['a_mode'] = stats.mode(self.results['cell_ar_best_cluster'][:,0])[0][0]
        self.results['b_mode'] = stats.mode(self.results['cell_ar_best_cluster'][:,1])[0][0]
        self.results['c_mode'] = stats.mode(self.results['cell_ar_best_cluster'][:,2])[0][0]
        self.results['al_mode'] = stats.mode(self.results['cell_ar_best_cluster'][:,3])[0][0]
        self.results['be_mode'] = stats.mode(self.results['cell_ar_best_cluster'][:,4])[0][0]
        self.results['ga_mode'] = stats.mode(self.results['cell_ar_best_cluster'][:,5])[0][0]

        logger.info('Cell selection # HKLs: %d' %len(self.results['cell_select']))

        self.results['dendo'] = sch.dendrogram(Y, p=10, labels=self.results['data_points'], no_plot=True)
        #logger.info('dendrogram switched off')

        return

    def cell_analysis(self, inData):
        self.setter(inData)
        tmp_fileList = [] #Local variable - a list of lists

        for i in range(len(self.results['hkllist'])):
            unit_cell = self.get_cells(self.results['hkllist'][i])


            self.results['cell_ar'][i,0] = unit_cell['a']
            self.results['cell_ar'][i,1] = unit_cell['b']
            self.results['cell_ar'][i,2] = unit_cell['c']
            self.results['cell_ar'][i,3] = unit_cell['al']
            self.results['cell_ar'][i,4] = unit_cell['be']
            self.results['cell_ar'][i,5] = unit_cell['ga']

        a_ohne, a_indx = self.reject_outlier(self.results['cell_ar'][:,0])
        b_ohne, b_indx = self.reject_outlier(self.results['cell_ar'][:,1])
        c_ohne, c_indx = self.reject_outlier(self.results['cell_ar'][:,2])
        al_ohne, al_indx = self.reject_outlier(self.results['cell_ar'][:,3])
        be_ohne, be_indx = self.reject_outlier(self.results['cell_ar'][:,4])
        ga_ohne, ga_indx = self.reject_outlier(self.results['cell_ar'][:,5])

        self.results['cell_ar'] = np.column_stack(a_ohne, b_ohne, c_ohne, al_ohne, be_ohne, ga_ohne)

        # calculate median values for unit cell parameters..
        self.results['a_mode'] = stats.mode(a_ohne)[0][0]
        self.results['b_mode'] = stats.mode(b_ohne)[0][0]
        self.results['c_mode'] = stats.mode(c_ohne)[0][0]
        self.results['al_mode'] = stats.mode(al_ohne)[0][0]
        self.results['be_mode'] = stats.mode(be_ohne)[0][0]
        self.results['ga_mode'] = stats.mode(ga_ohne)[0][0]

        for index in (a_indx, b_indx, c_indx, al_indx, be_indx, ga_indx):
            if index: # no outlier may lead to an empty array, avoiding it..
                tmp_fileList.append([self.results['hkllist'][i] for i in index])

        cell_select = set(tmp_fileList[0]) #use set method to find common HKL files after outlier rejection
        for s in tmp_fileList[1:2]:
            cell_select.intersection_update(s)
            self.results['cell_select'] = list(cell_select)

        return

    def dict_for_histogram(self):
        self.results['hist_dict'] = dict()
        if self.results['cell_ar'].size == 0:
            msg = "Cell selection did not work, so cell_array is empty"
            logger.info('MSG:{}'.format(msg))

        self.results['hist_dict']['a'] = list(self.results['cell_ar'][:,0])
        self.results['hist_dict']['b'] = list(self.results['cell_ar'][:,1])
        self.results['hist_dict']['c'] = list(self.results['cell_ar'][:,2])
        self.results['hist_dict']['al'] = list(self.results['cell_ar'][:,3])
        self.results['hist_dict']['be'] = list(self.results['cell_ar'][:,4])
        self.results['hist_dict']['ga'] = list(self.results['cell_ar'][:,5])

        return


    def cell_histogram(self, inData):
        self.cell_analysis(inData)

        if len(self.results['hklList']) == 0:
            err = ValueError("xscaling with ISa-selection did not work, reduce isa-cutoff")
            logger.error(err)
            self.setFailure()

        fig, axs = plt.subplots(2,3, squeeze=True, facecolor='w', edgecolor='k')
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace = 0.4, wspace=0.3)
        axs = axs.flatten()


        axs[0].hist(self.results['cell_ar'][:,0], 50, histtype='step')
        axs[0].set_xlabel('a-axis(Ang)', fontweight='bold')
        axs[0].set_ylabel('frequency', fontweight='bold')
        axs[0].set_title("a=%5.2f A" %self.results['a_mode'])


        axs[1].hist(self.results['cell_ar'][:,1], 50, histtype='step')
        axs[1].set_xlabel('b-axis(Ang)', fontweight='bold')
        axs[1].set_ylabel('frequency', fontweight='bold')
        axs[1].set_title("b=%5.2f A" %self.results['b_mode'])


        axs[2].hist(self.results['cell_ar'][:,2], 50, histtype='step')
        axs[2].set_xlabel('c-axis(Ang)', fontweight='bold')
        axs[2].set_ylabel('frequency', fontweight='bold')
        axs[2].set_title("c=%5.2f A" %self.results['c_mode'])


        axs[3].hist(self.results['cell_ar'][:,3],50, histtype='step')
        axs[3].set_xlabel('alpha(deg)', fontweight='bold')
        axs[3].set_ylabel('frequency', fontweight='bold')
        axs[3].set_title("alpha=%5.2f A" %self.results['al_mode'])

        axs[4].hist(self.results['cell_ar'][:,4], 50, histtype='step')
        axs[4].set_xlabel('beta(deg)', fontweight='bold')
        axs[4].set_ylabel('frequency', fontweight='bold')
        axs[4].set_title("beta=%5.2f A" %self.results['be_mode'])

        axs[5].hist(self.results['cell_ar'][:,5], 50, histtype='step')
        axs[5].set_xlabel('gamma[deg]', fontweight='bold')
        axs[5].set_ylabel('frequency', fontweight='bold')
        axs[5].set_title("gamma=%5.2f A" %self.results['ga_mode'])

        plt.savefig("cell-histogram.png")
        '''
        # run_command will not work because this is a python command and not bash shell command.
        try:
            run_command("Scale&Merge", os.getcwd(), os.environ['USER'], cmd, 'merge.log')
        except Exception:
            sub.call(cmd, shell=True)
        return self.status
        '''


def finder():
    root = sys.argv[1]
    import glob
    paths = glob.glob(os.path.join(root, 'minisets_*/xtal_*/XDS_ASCII.HKL'))
    return paths
'''
import glob
paths = glob.glob('/Users/shibom/work/CY_IMISX/all_pepT_sad_v2/xtal_*.HKL')
cell = Cell(sorted(paths[0:100]))
cell.clustering()
'''
