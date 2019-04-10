# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import numpy as np
import pickle

from pyworkflow.em.viewers import ObjectView, DataViewer, views
from pyworkflow.em.viewers.showj import MODE, MODE_MD, ORDER, VISIBLE, RENDER, SORT_BY
from pyworkflow.protocol.params import IntParam, LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer

# import xmippLib
from xmipp3.protocols.protocol_consensus_classes3D import XmippProtConsensusClasses3D
from .plotter import XmippPlotter

class XmippConsensus3DClassesViewer(ProtocolViewer):
    """         
        Viewer for the 'Xmipp - consensus classes 3D'.\n
        Select those particles/cooridantes with high (close to 1.0) 'zScoreDeepLearning1' value and save them.
        The Histogram may help you to decide a threshold.
    """
    _label = 'viewer consensus 3D classes'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtConsensusClasses3D]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('noteViz', LabelParam, label="\n")
        form.addParam('visualizeDendogram', LabelParam, important=True,
                      label="Vizualoze dendogram",
                      help="Visualize dendogram. It helps to decide how many "
                           "independent classes are in the data set.")
        form.addParam('visualizeAllClasses', LabelParam,
                      label="Visualize all node classes",
                      help="Visualize all node classes. The associated volume is "
                           "only for illustrative purposes and it should be "
                           "recalculated.")
        form.addParam('visualizeNClasses', IntParam, default=3,
                      label="Visualize N classes",
                      help="Visualize N classes via joining nearest nodes."
                           "The associated volume is only for illustrative "
                           "purposes and it should be recalculated.")

    def _getVisualizeDict(self):
        return {'visualizeDendogram': self._visualizeDendogram,
                'visualizeAllClasses': self._visualizeClasses}

    def _visualizeClasses(self, e=None):
        _views = []

        labels = 'id enabled _size _representative._filename _xmipp_nodeLabel'

        obj = self.protocol.nodeIntersectionClasses
        _views.append(views.Classes3DView(self._project, obj.strId(),
                                          obj.getFileName(),
                                          viewParams={VISIBLE: labels}))

        return _views

    def _visualizeDendogram(self, e=None):
        views = []

        classDistFn = 'classDistances.pkl'
        with open(self.protocol._getExtraPath(classDistFn), 'rb') as pklFile:
            classDistDict = pickle.load(pklFile)

        classDistances = classDistDict['classDistances']
        nodesLabels = classDistDict['nodesLabels']
        nInputs = len(self.protocol.inputMultiClasses)

        # ijk = ['i', 'j', 'k', 'l', 'm', 'n']
        intersectionStr = u' \u2229 '.join(["Par%s"%(i+1) for i in
                                          range(0, nInputs)])
        Y = classDistances.values()
        Z = hierarchy.linkage(np.asarray(Y), 'single')
        xplotter = XmippPlotter(windowTitle="Dendogram")
        xplotter.createSubPlot(u"Nodes: %s" % intersectionStr,
                               "",
                               "Distance between nodes")

        labels = []
        for x in nodesLabels:
            xList = x.split(':')
            label = ''
            for i, y in enumerate(xList):
                label += 'Par%d:cl '%(i+1) + y + '\n'  # + u' \u2229\n'
            labels.append(label)

        nplabels = np.asarray(labels)
        dn = hierarchy.dendrogram(Z, labels=nplabels)
        # for k, v in dn.iteritems():
        #     print("%s (%s)" % (k, type(k)))
        #     if k == 'icoord' or k == 'dcoord':
        #         for vv in v:
        #             print(vv)
        #     else:
        #         print("%s: %s" % (k, v))
        #     print('-----')
        print
        print([nodesLabels[id] for id in dn['leaves']])
        print

        plt.style.use("seaborn-whitegrid")

        views.append(xplotter)

        return views
