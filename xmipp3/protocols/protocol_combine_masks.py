# *****************************************************************************
# *
# * Authors:     David Herreros Calero         dherreros@cnb.csic.es
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
# *****************************************************************************

import numpy as np

from pwem.objects import Volume
from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler

import pyworkflow.protocol.params as params


class XmippProtCombineMasks(ProtAnalysis3D):
    """Combines several segmented regions in a map in a single mask"""
    _label = 'combine masks'

    # --------------------------- DEFINE param functions ---------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMasks', params.PointerParam, pointerClass='SetOfVolumes', label='Input Masks',
                      important=True, help='Input maksk stored in a SetOfVolumes to be combined')
        form.addParam('binary', params.BooleanParam, label='Binarize mask', default=True,
                      help='If False, a mask containing identifiers for each region will be returned')

    # --------------------------- INSERT steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('combineMasks')
        self._insertFunctionStep('createOutputStep')

    def combineMasks(self):
        dim = self.inputMasks.get().getDim()
        binary = self.binary.get()
        ih = ImageHandler()
        outMask = ih.createImage()
        mask_combined = np.zeros(dim, float)
        for idp, mask_img in enumerate(self.inputMasks.get().iterItems()):
            mask = ih.read(mask_img).getData()
            if binary:
                mask_combined += mask
            else:
                idm = self.intersectionMask(np.copy(mask_combined), np.copy(mask))
                mask = mask * (idp + 1)
                mask_combined = mask_combined + mask
                mask_combined[idm] = (idp + 1)
        outMask.setData(mask_combined)
        ih.write(outMask, self._getExtraPath('Combined_Mask.mrc'))


    def createOutputStep(self):
        volume = Volume()
        volume.setSamplingRate(self.inputMasks.get().getSamplingRate())
        volume.setLocation(self._getExtraPath('Combined_Mask.mrc'))
        self._defineOutputs(outputMask=volume)
        self._defineSourceRelation(self.inputMasks, volume)

    # --------------------------- DEFINE utils functions ----------------------
    def intersectionMask(self, minMask, maxMask):
        minMask[minMask != 0.0] = 1
        maxMask[maxMask != 0.0] = 1
        sumMask = minMask + maxMask
        idm = sumMask == 2
        return idm

    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            msg = ("Maks combined to %s\n" % self.outputMask.getFileName())
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Mask not ready yet")
        return methodsMsgs

    def _summary(self):
        summary = []
        summary.append("Number of input masks provided: %d\n"
                       % len(self.inputMasks.get()))
        if self.getOutputsSize() >= 1:
            msg = ("Maks combined to %s\n" % self.outputMask.getFileName())
            summary.append(msg)
        else:
            summary.append("Maks not ready yet.")
        return summary