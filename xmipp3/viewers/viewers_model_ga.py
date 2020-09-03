# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * BCU, Centro Nacional de Biotecnologia, CSIC
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

import os, glob
import numpy as np
from matplotlib import cm

import pyworkflow.viewer as pwviewer

import pwem.viewers.views as vi
from pwem.viewers.viewer_chimera import ChimeraView

from pwem.emlib.image import ImageHandler

from ..protocols.protocol_model_ga import XmippProtModelGA

class XmippProtModelGAViewer(pwviewer.Viewer):
    """ Wrapper to visualize the output of GA modeling
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [XmippProtModelGA]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, XmippProtModelGA):
            num_regions = len(self.protocol.inputSeqs)
            regions_id = [x+1 for x in range(num_regions)]
            jet = cm.get_cmap('viridis')(np.linspace(0, 1, num_regions+1))
            jet[0] = [0, 0, 0, 1.]
            jet = tuple(tuple(x) for x in jet)
            levels = tuple(tuple((x-0.01, 1)) for x in regions_id)
            levels = ((0, 0),) + levels
            self.new_ap = {
                'regions': {
                    'image_levels': levels,
                    'image_colors': jet,
                    'dim_transparent_voxels': True,
                    'transparency_depth': 0.5,
                    'image_brightness_factor': 1,
                    }
                }
            filePath = self.chimeraViewFile()
            views.append(ChimeraView(filePath))

        return views

    def chimeraViewFile(self):
        filePath = self.protocol._getTmpPath('viewChimera.py')
        f = open(filePath, "w")
        f.write('import ast\n')
        f.write('from chimerax.core.commands import run\n')
        f.write('from chimerax.map.colortables import colortables\n')
        f.write('colortables._appearances.update(%s)\n' % self.new_ap)
        f.write('run(session, "open %s")\n' % self.protocol.outputMasks.getFirstItem().getFileName())
        f.write('run(session, "volume style image appearance regions step 1")\n')
        f.close()
        return filePath
