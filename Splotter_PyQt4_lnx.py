#!/usr/bin/env python

import matplotlib
#from PyQt4.QtCore import QStringList
matplotlib.use('Qt4Agg')
#import sys
import os
import touchstone as ts
from PyQt4 import QtGui #, QtCore
import pylab
pylab.ion()
from numpy import size, array, log10, abs, angle, pi, real, imag, unwrap
import matplotlib.pyplot as plt
from Jsonify import *
from Json_Plotter import *


def dB(x):
    return 20.*log10(abs(x))

def magLinear(x):
    return abs(x)

def phaseDeg(x):
    return unwrap(angle(x)) * 180./pi

def phaseRad(x):
    return unwrap(angle(x))

#########################################################################################

class sParamPlotter_(QtGui.QMainWindow):

    def __init__(self, session):

        super(sParamPlotter_, self).__init__()
        self.sessionFileName = session
        self.Json_Session = JsonifySession(session)
        self.models = []
        self.traces = {}
        self.lastPath = None
        self.initialize()

    def initialize(self):
        # create a menu bar for the window
        menubar = self.menuBar()
        # Enable user to load a .s4p file
        file_action = QtGui.QAction(QtGui.QIcon('file.png'), '&Load .s*p File', self)
        file_action.setShortcut('Ctrl+O')
        file_action.setStatusTip('Import data from a file')
        file_action.triggered.connect(self.loadFile)
        file_menu = menubar.addMenu('&Load File')
        file_menu.addAction(file_action)
        # Enable user to plot Loaded traces
        plot_action = QtGui.QAction(QtGui.QIcon('plot.png'), '&Plot', self)
        plot_action.setShortcut('Ctrl+P')
        plot_action.setStatusTip('Plot a loaded trace')
        plot_action.triggered.connect(self.plotSelected)
        plot_menu = menubar.addMenu('&Plot')
        plot_menu.addAction(plot_action)
        # Enable to load and plot json files
        json_action = QtGui.QAction(QtGui.QIcon('file.png'), '&Plot .json File', self)
        json_action.setShortcut('Ctrl+J')
        json_action.setStatusTip('Display plot from .json file')
        json_action.triggered.connect(self.loadJson)
        plot_menu.addAction(json_action)
        # Let user clear the tree
        clear_action = QtGui.QAction(QtGui.QIcon('clear.png'), '&Clear', self)
        clear_action.setShortcut('Ctrl+C')
        clear_action.setStatusTip('Clear all imported traces')
        clear_action.triggered.connect(self.clear_traces)
        clear_menu = menubar.addMenu('&Clear')
        clear_menu.addAction(clear_action)
        # Add widgets and display window
        self.w = Label_Window(self)
        self.setCentralWidget(self.w)
        self.setGeometry(300,300,500,500)
        self.setWindowTitle('Splotter')
        self.show()

    def loadJson(self):
        """
        loads pre-plotted json plotSelected
        """
        if self.lastPath is None:
            fullpath = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File', '.', '*.json'))
            #fullpath = str(QtGui.QFileDialog.getOpenFileName(self, caption='Open File', filter='*.json'))
        else:
            fullpath = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File', self.lastPath, '*.json'))
        if fullpath == '':
            print('No *.json file selected')
        else:
            text = load_json_file(fullpath)
            plot_json(text)

    def clear_traces(self):
        """
        Removes all entries in the tree menu
        """
        self.traces.clear()
        self.refreshList()

    def refreshList(self):
        """
        Refreshes the tree menu
        """
        self.w.refresh_tree(self.traces)
        self.w.tree_menu.sortItems(0,0) # Sorts the Tree menu in Acending order based on label

    def get_x_fn(self):
        return str(self.w.x_func_edit.text())

    def get_y_fn(self):
        return str(self.w.y_func_edit.text())

    def get_x_label(self):
        return str(self.w.x_label_edit.text())

    def get_y_label(self):
        return str(self.w.y_label_edit.text())

    def get_title(self):
        return str(self.w.title_edit.text())

    def get_selected_items(self):
        """
        Returns the all selected leaves in the tree menu
        """
        plot_list = []
        root = self.w.tree_menu.invisibleRootItem()
        self.w.tree_menu.get_all_selected_leaves(root, plot_list, root)
        return plot_list

    def loadFile(self):
        """
        Load a touchstone file. Put each parameter in the Tree Menu
        """

        if self.lastPath is None:
            #fullpath = str(QtGui.QFileDialog.getOpenFileName(self, caption='Open File', filter='*.s*p'))
            #fullpath = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File', '.', 'Touchstone (*.s*p)'))
            pathlist = QtGui.QFileDialog.getOpenFileNames(self, 'Open File', '.', 'Touchstone (*.s*p)')            
            
        else:
            #fullpath = str(QtGui.QFileDialog.getOpenFileName(self, caption='Open File', directory=self.lastPath, filter='*.s*p'))
            pathlist = QtGui.QFileDialog.getOpenFileNames(self, 'Open File', self.lastPath, 'Touchstone (*.s*p)')

        
        for pi in pathlist: # Iterates over the QStringList provided from the file dialog.
            fullpath = str(pi)
            dirname, filename = os.path.split(fullpath)
            self.lastPath = dirname
            labelName, useBalSparms, useOEOrdering = getEntryPlusCheckbuttonValues_(self, 'Label For File', EntryDefault=filename, MixedModeDefault=0, OEOrderingDefault=1)
    
            SModel = None
            model_index = len(self.models)
            if useBalSparms == 1:
                SModel = ts.readmodel(fullpath, oe_ordering=bool(useOEOrdering))
                SModel.toMixedMode()
            else:
                SModel = ts.readmodel(fullpath, oe_ordering=bool(useOEOrdering))
    
            N = SModel.getPorts()
            self.models.append(SModel)
    
            portLetterList = ['d',]*int(N/2) + ['c',]*int(N/2)
            portNumberList = list(range(int(N/2))) + list(range(int(N/2)))
            
            for ii in range(0,N):
                for jj in range(0,N):
                    if useBalSparms == 0:
                        paramLabel = 'S%d,%d' % (ii+1,jj+1)
                        mode = ''
                    else:
                        paramLabel = 'S%s%s%s,%s' % (portLetterList[ii],portLetterList[jj], portNumberList[ii]+1,portNumberList[jj]+1)
                        mode = 'S%s%s' % (portLetterList[ii],portLetterList[jj])
                    trace_label = labelName + ' - ' + paramLabel
                    self.traces[trace_label] = {'label': trace_label,
                                                    'filename': filename,
                                                    'mixedmode': useBalSparms,
                                                    'mode': mode,
                                                    'oe_ordering': useOEOrdering,
                                                    'model_index': model_index,
                                                    'i': ii,
                                                    'j': jj,
                                'paramLabel': paramLabel,
                                'label name': labelName}
        self.refreshList()

    def plotSelected(self):

        """
        Get the selected traces from the traceList.  Then, plot the
        data in matplotlib window with the pylab interface.
        """
        plot_dict = PlotDict(self.get_title(), self.get_y_label(), self.get_x_label())

        xfcnStr = self.get_x_fn()
        yfcnStr = self.get_y_fn()

        xfcn = eval('lambda x, y: %s' % xfcnStr)
        yfcn = eval('lambda x, y: %s' % yfcnStr)

        plotObj = plt.figure()
        axes = plotObj.add_subplot(1,1,1)
        selected_plots = self.get_selected_items()
        for item in selected_plots:
            item = str(item)
            Data = self.traces[item]
            M = self.models[Data['model_index']]
            x = M.farray # Data['frequency']
            y = M.Sarray[:,Data['i'],Data['j']] # Data['data']
            axes.plot(xfcn(x, y), yfcn(x, y), label=Data['label'], picker=5)
            plot_dict.add_line(xfcn(x,y).tolist(), yfcn(x,y).tolist(), Data['label'])

        axes.grid(True)
        axes.legend(loc='best')
        axes.set_xlabel(self.get_x_label())
        axes.set_ylabel(self.get_y_label())
        axes.set_title(self.get_title())

        plot_dict.set_limits(axes.get_xlim(), axes.get_ylim())

        plt.show()

        self.Json_Session.add_plot(plot_dict)
        self.Json_Session.save()

########################################################################################

class Label_Window(QtGui.QWidget):

    """
    Creates a widget where labels, single line text entries, and the tree menu will live
    """

    def __init__(self, parent):

        super(Label_Window, self).__init__(parent)
        self.tree_menu = None
        self.init_widgets()

    def init_widgets(self):

        title = QtGui.QLabel('Title')
        x_label = QtGui.QLabel('X Label')
        y_label = QtGui.QLabel('Y Label')
        X_func = QtGui.QLabel('X Function')
        Y_func = QtGui.QLabel('Y Function')

        self.title_edit = QtGui.QLineEdit()
        self.title_edit.setText('')
        self.x_label_edit = QtGui.QLineEdit()
        self.x_label_edit.setText('Frequency (GHz)')
        self.y_label_edit = QtGui.QLineEdit()
        self.y_label_edit.setText('dB')
        self.x_func_edit = QtGui.QLineEdit()
        self.x_func_edit.setText('x*1.e-9')
        self.y_func_edit = QtGui.QLineEdit()
        self.y_func_edit.setText('20.*log10(abs(y))')

        grid = QtGui.QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(title, 1, 0)
        grid.addWidget(self.title_edit, 1, 1)

        grid.addWidget(x_label, 2, 0)
        grid.addWidget(self.x_label_edit, 2, 1)

        grid.addWidget(y_label, 3, 0)
        grid.addWidget(self.y_label_edit, 3, 1)

        grid.addWidget(X_func, 4, 0)
        grid.addWidget(self.x_func_edit, 4, 1)

        grid.addWidget(Y_func, 5, 0)
        grid.addWidget(self.y_func_edit, 5, 1)

        self.tree_menu = TreeMenu(self)
        grid.addWidget(self.tree_menu, 6, 0, 5, 5)

        self.setLayout(grid)

    def refresh_tree(self, traces):
        #remove all entries from the current tree menu
        self.tree_menu.remove_all()
        for item in traces.keys():
            data = traces[item]
            self.tree_menu.add_item(data)

#############################################################################

class TreeMenu(QtGui.QTreeWidget):

    """
    A collapsable tree menu. Will Always have 3 levels
    """

    def __init__(self, parent):
        super(TreeMenu,self).__init__(parent)
        self.setColumnCount(1)
        self.header().close()  # Removes Unneccessary header
        # Below enables Multiple List Items to be selected at once
        self.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)

    def remove_all(self):
        # root = self.invisibleRootItem()
        # while (root.childCount() > 0):
        #       root.removeChild(root.child(0))
        self.clear()

    def add_item(self, item_dict):
        label = item_dict['label name']
        paramLabel = item_dict['paramLabel']
        mixedmode = item_dict['mixedmode']

        root = self.invisibleRootItem()
        item = Node()
        item.setText(0, label)

        if mixedmode:
            mode = item_dict['mode']
            child1 = Node()
            child1.setText(0, mode)
        else:
            mode = 'S%d,*' % (item_dict['i']+1,)
            child1 = Node()
            child1.setText(0, mode)

        child2 = QtGui.QTreeWidgetItem()
        child2.setText(0, paramLabel)

        lvl0_flag, item = self.is_child(root, item)
        lvl1_flag, child1 = self.is_child(item, child1)

        if not lvl1_flag:
            child1 = Node(parent=item)
            child1.setText(0, mode)

        lvl2_flag, child2 = self.is_child(child1, child2)
        if not lvl2_flag:
            child2 = Node(parent=child1, full_label=item_dict['label'], is_leaf=True)
            child2.setText(0, paramLabel)

        if not lvl0_flag:
            self.addTopLevelItem(item)

        child1.add_leaf(item_dict['label'])
        item.add_leaf(item_dict['label'])

    def is_child(self, parent, item):
        """
        Checks if QTreeWidgetItem item is a child to parent
        """
        child_count = parent.childCount()
        for i in range(child_count):
            child = parent.child(i)
            if child.text(0) == item.text(0):
                return True, child
        return False, item

    def get_all_selected_leaves(self, root, item_list, top_item):
        """
        returns all leaves that are selected
        """
        child_count = root.childCount()
        for i in range(child_count):
            child = root.child(i)
            if child.isSelected():
                if child.is_leaf:
                    if child.full_label not in item_list:
                        item_list.append(child.full_label)
                else:
                    for leaf in child.leaf_list:
                        if leaf not in item_list:
                            item_list.append(leaf)
            self.get_all_selected_leaves(child, item_list, top_item)

######################################################################################

class Node(QtGui.QTreeWidgetItem):
    """
    An extension of the QTreeWidgetItem Class that allows for additional info to be stored in the nodes
    """
    def __init__(self, full_label='', parent=None, is_leaf=False):
        super(Node,self).__init__(parent)
        self.is_leaf = is_leaf
        self.full_label = full_label
        self.leaf_list = []

    def add_leaf(self, leaf):
        if not leaf in self.leaf_list:
            self.leaf_list.append(leaf)

##########################################################################################

class EntryPlusCheckButtonDialog_(QtGui.QDialog):

    def __init__(self, parent=None, title=None, EntryDefault = '', MixedModeDefault = 0, OEOrderingDefault = 0):
        super(EntryPlusCheckButtonDialog_, self).__init__(parent)
        self.title = title
        self.EntryDefault = EntryDefault
        self.MixedModeDefault = MixedModeDefault
        self.OEOrderingDefault = OEOrderingDefault
        self.body()

    def body(self):

        grid = QtGui.QGridLayout()
        grid.setSpacing(10)

        label = QtGui.QLabel(self.title)
        self.labelValue = QtGui.QLineEdit()
        self.labelValue.setText(self.EntryDefault)
        grid.addWidget(label, 1, 0)
        grid.addWidget(self.labelValue, 1, 1)

        self.balsparmValue = QtGui.QCheckBox('Balanced S-parameters', self)
        grid.addWidget(self.balsparmValue, 2, 0)

        self.oeorderingValue = QtGui.QCheckBox('Odd-Even Ordering', self)
        self.oeorderingValue.setCheckState(2)  # Starts out as checked
        grid.addWidget(self.oeorderingValue, 3,0)

        accept_button = QtGui.QPushButton('Okay')
        accept_button.clicked.connect(self.accept)
        grid.addWidget(accept_button, 4, 0)

        cancel_button = QtGui.QPushButton('Cancel')
        cancel_button.clicked.connect(self.reject)
        grid.addWidget(cancel_button, 4, 1)

        self.setLayout(grid)
        self.show()


def getEntryPlusCheckbuttonValues_(parent=None, title = None, EntryDefault = '', MixedModeDefault = 0, OEOrderingDefault = 0):
    labelDialog =  EntryPlusCheckButtonDialog_(parent=parent, title=title,
            EntryDefault=EntryDefault,
            MixedModeDefault=MixedModeDefault,
            OEOrderingDefault=OEOrderingDefault)

    if labelDialog.exec_():
        return str(labelDialog.labelValue.text()), int(labelDialog.balsparmValue.isChecked()), int(labelDialog.oeorderingValue.isChecked())
    else:
        return str(EntryDefault), int(MixedModeDefault), int(OEOrderingDefault)

#####################################################################################

codeImportStr = """import cramsens.touchstone as ts
import pylab
from numpy import array, log10, abs, angle, pi, real, imag, unwrap

def dB(x):
    return 20.*log10(abs(x))

def magLinear(x):
    return abs(x)

def phaseDeg(x):
    return angle(x) * 180./pi

def phaseRad(x):
    return angle(x)
"""


def main():

    app = QtGui.QApplication(sys.argv)
    print (os.popen('pwd').read())

    argv = sys.argv[1:]
    if len(argv) > 0:
        session=argv[0]
    else:
        session='Splotter.Session.json'

    e = sParamPlotter_(session)
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
