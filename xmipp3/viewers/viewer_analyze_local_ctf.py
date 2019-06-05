# **************************************************************************
# *
# * Authors:     Carlos Oscar Sanchez Sorzano
# *              Estrella Fernandez Gimenez
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

import numpy as np
import xmippLib
from os.path import exists

from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
from pyworkflow.protocol.params import IntParam, LabelParam
from pyworkflow.em.viewers import ObjectView
import pyworkflow.em.viewers.showj as showj


from xmipp3.protocols.protocol_analyze_local_ctf import XmippProtAnalyzeLocalCTF
from .plotter import XmippPlotter


class XmippAnalyzeLocalCTFViewer(ProtocolViewer):
    """ Visualization of the output of analyze_local_ctf protocol
    """
    _label = 'viewer analyze local defocus'
    _targets = [XmippProtAnalyzeLocalCTF]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        self._data = None

    def getData(self):
        if self._data is None:
            self._data = self.loadData()
        return self._data

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayR2', LabelParam, default=False, label='Display micrograph R2')
        form.addParam('displayLocalDefocus', IntParam, default=1,
                      label='Display local defocus of micrograph:',
                      help="""Type the ID of the micrograph to see particle local defocus of the selected micrograph. 
                      It is possible that not all the micrographs are available""")

    def _getVisualizeDict(self):
        return {'displayR2' : self._viewR2,
                'displayLocalDefocus': self._viewLocalDefocus,
                }

    def _viewLocalDefocus(self, paramName=None):
        views=[]
        fnDefoci="%s"%(self.protocol._getExtraPath("micrographDefoci.xmd"))
        if exists(fnDefoci):
            try:
                mdPoints = xmippLib.MetaData("mic_%d@%s"%(self.displayLocalDefocus.get(),fnDefoci))

                x = mdPoints.getColumnValues(xmippLib.MDL_XCOOR)
                y = mdPoints.getColumnValues(xmippLib.MDL_YCOOR)
                defocusA = mdPoints.getColumnValues(xmippLib.MDL_CTF_DEFOCUSA)
                residuals = mdPoints.getColumnValues(xmippLib.MDL_CTF_DEFOCUS_RESIDUAL)

                title="Micrograph %d defocus"%self.displayLocalDefocus.get()
                xplotter = XmippPlotter(windowTitle=title)
                a = xplotter.createSubPlot(title, 'x', 'y', projection='3d')
                a.set_zlabel('Defocus')
                a.scatter(x, y, defocusA, c='r', marker='o')
                a.scatter(x, y, np.asarray(defocusA)-np.asarray(residuals), c='b', marker='^')
                legends = ['Avg. defocus','Adjusted defocus']
                xplotter.showLegend(legends, loc=1)
                views.append(xplotter)
            except:
                pass
        return views

    def _viewR2(self, paramName=None):
        views = []
        if hasattr(self.protocol, "outputMicrographs"):
            obj = self.protocol.outputMicrographs
            fn = obj.getFileName()
            labels = 'id _filename _xmipp_ctfDefocusR2'
            views.append(ObjectView(self._project, obj.strId(), fn,
                                          viewParams={showj.ORDER: labels,
                                                      showj.VISIBLE: labels,
                                                      showj.MODE: showj.MODE_MD,
                                                      showj.RENDER:'_filename'}))
        return views



