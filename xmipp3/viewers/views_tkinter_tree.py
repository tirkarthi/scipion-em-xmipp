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
import os
import tables
from matplotlib import cm
import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.viewers.viewer_chimera import Chimera

import pyworkflow.utils as pwutils
from pyworkflow.gui import TreeProvider, FileTreeProvider, ObjectBrowser


class SeggerTreeProvider(FileTreeProvider):
    """ Populate Tree from SetOfVolumes. """

    def __init__(self, volumeList, path):
        TreeProvider.__init__(self)
        self.volumeList = volumeList
        self._path = path

    def getColumns(self):
        return [('Volume', 150), ('Score', 50)]

    def getObjectInfo(self, volume):
        volumeName = os.path.basename(volume.getFileName())
        volumeName = os.path.splitext(volumeName)[0]
        return {'key': volumeName, 'parent': None,
                'text': volumeName, 'values': (round(volume.score.get(), 6))}

    def getObjectActions(self, obj):
        return []

    def getObjectPreview(self, obj):
        filename = obj.getFileName()
        fileHandler = self.getFileHandler(obj)
        return fileHandler._getImagePreview(filename), fileHandler._getImageString(filename)

    def _getObjectList(self):
        """Retrieve the object list"""
        return self.volumeList

    def getObjects(self):
        objList = self._getObjectList()
        return objList


class SeggerBrowser(ObjectBrowser):
    """
    This class extend from ListDialog to allow calling
    an Segger subprocess from a list of Volumes.
    """
    def __init__(self, window, provider, path):
        self.path = path
        ObjectBrowser.__init__(self, window.root, provider)
        self.tree.itemDoubleClick = self._itemDoubleClick

    def _itemDoubleClick(self, obj):
        self.volume = obj
        self.lanchSeggerForVolume(self.volume)

    def lanchSeggerForVolume(self, volume):
        segFile = self.idMask2Segger(volume)
        filePath = self.chimeraViewFile(segFile)
        program = Chimera.getProgram()
        pwutils.runJob(None, program, filePath, env=Chimera.getEnviron())

    def idMask2Segger(self, volume):
        ih = ImageHandler()
        filePath = volume.getFileName()
        image = ih.read(filePath)
        mask = np.squeeze(image.getData())
        mask = mask.astype(np.uint32)
        outPath = os.path.join(self.path, 'mask.seg')
        h5file = tables.open_file(outPath, mode='w')
        root = h5file.root
        a = root._v_attrs
        a.format = 'segger'
        a.format_version = 2
        a.name = os.path.basename(filePath)
        a.map_size = np.array(mask.shape, np.int32)
        a.map_path = filePath
        a.map_level = 0.01
        a.ijk_to_xyz_trasform = np.array([[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.]], np.float32)
        atom = tables.Atom.from_dtype(mask.dtype)
        filters = tables.Filters(complevel=5, complib='zlib')
        ma = h5file.create_carray(root, 'mask', atom, mask.shape, filters=filters)
        ma[:] = mask
        num_regions = int(np.amax(mask))
        rids = np.array([r + 1.0 for r in range(num_regions)], np.int32)
        h5file.create_array(root, 'region_ids', rids)
        rcolors = cm.get_cmap('viridis')(np.linspace(0, 1, num_regions))
        h5file.create_array(root, 'region_colors', rcolors)
        slev = np.array([0] * num_regions, np.float32)
        h5file.create_array(root, 'smoothing_levels', slev)
        pids = np.array([0] * num_regions, np.int32)
        h5file.create_array(root, 'parent_ids', pids)
        refpts = np.array([np.asarray(
            [np.where(mask == (r + 1))[0][0], np.where(mask == (r + 1))[1][0], np.where(mask == (r + 1))[2][0]]) for
            r in range(num_regions)], np.float32)
        h5file.create_array(root, 'ref_points', refpts)
        h5file.close()
        return outPath

    def chimeraViewFile(self, segFile):
        filePath = os.path.join(self.path, 'viewChimera.py')
        f = open(filePath, "w")
        f.write('from chimerax.core.commands import run\n')
        f.write('run(session, "open %s")\n' % os.path.abspath(segFile))
        f.write('run(session, "vol all hide")\n')
        f.close()
        return filePath
