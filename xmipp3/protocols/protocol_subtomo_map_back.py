# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

from pyworkflow.em import ALIGN_3D
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, FloatParam
from tomo.objects import Tomogram
from tomo.protocols import ProtTomoBase


from xmipp3.convert import alignmentToRow
import xmippLib


class XmippProtSubtomoMapBack(EMProtocol, ProtTomoBase):
    """ This protocol takes a tomogram, a reference subtomogram and a metadata with geometrical parameters
   (x,y,z) and places the reference subtomogram on the tomogram at the designated locations (map back).
   It has different representation options."""

    _label = 'map back subtomos'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputClasses', PointerParam, pointerClass="SetOfClassesSubTomograms", label='Class',
                      help="Subtomogram class from which the coordinates of the subtomograms and the reference will be "
                           "used. It should be a SetOfClassesSubTomograms with just 1 item.")
        form.addParam('inputTomograms', PointerParam, pointerClass="SetOfTomograms",
                      label='Original tomograms', help="Original tomograms from which the subtomograms were extracted")
        form.addParam('invertContrast', BooleanParam, default=False, label='Invert reference contrast',
                      help="Invert the contrast if the reference is black over a white background.  Xmipp, Spider, "
                           "Relion and Eman require white particles over a black background. ")
        form.addParam('paintingType', EnumParam,
                      choices=['Copy', 'Average', 'Highlight', 'Binarize'],
                      default=0, important=True,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Painting mode',
                      help='The program has several painting options:\n*Copy*: Copying the reference onto the tomogram.'
                           '\n*Average*: Setting the region occupied by the reference in the tomogram to the average '
                           'value of that region.\n*Highlight*: Add the reference multiplied by a constant to the '
                           'location specified.\n*Binarize*: Copy a binarized version of the reference onto the '
                           'tomogram.')
        form.addParam('removeBackground', BooleanParam, default=False, label='Remove background',
                      help="Set tomogram to 0", condition="paintingType == 0 or paintingType == 3")
        form.addParam('threshold', FloatParam, default=0.5, label='Threshold',
                      help="threshold applied to tomogram", condition="paintingType == 1 or paintingType == 3")
        form.addParam('constant', FloatParam, default=2, label='Multiplier',
                      help="constant to multiply the reference",
                      condition="paintingType == 2")

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInput')
        for subtomoClass in self.inputClasses.get():
            self._insertFunctionStep('runMapBack', subtomoClass.getObjId())
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions -------------------------------
    def convertInput(self):
        for tomo in self.inputTomograms.get().iterItems():
            img = ImageHandler()
            fnTomo = self._getExtraPath('tomogram_%d.mrc' % tomo.getObjId())
            img.convert(tomo, fnTomo)
            for classSubt in self.inputClasses.get().iterItems():
                cId = classSubt.getFirstItem().getClassId()
                fnRef = self._getExtraPath('reference%d.mrc' % cId)
                img.convert(classSubt.getRepresentative(), fnRef)
                if self.invertContrast.get() == True:
                    self.runJob("xmipp_image_operate", " -i %s  --mult -1" % fnRef)
            if self.paintingType.get() == 0 or self.paintingType.get() == 3:
                if self.removeBackground.get() == True:
                    self.runJob("xmipp_image_operate", " -i %s  --mult 0" % fnTomo)

    def runMapBack(self, classId):
        for tomo in self.inputTomograms.get().iterItems():
            TsSubtomo = self.inputClasses.get().getSamplingRate()
            TsTomo = tomo.getSamplingRate()
            scaleFactor = TsSubtomo/TsTomo
            mdGeometry = xmippLib.MetaData()
            for subtomo in self.inputClasses.get().getFirstItem().iterItems():
                if subtomo.getCoordinate3D().getVolId() == tomo.getObjId():
                    nRow = md.Row()
                    nRow.setValue(xmippLib.MDL_ITEM_ID, long(subtomo.getObjId()))
                    nRow.setValue(xmippLib.MDL_XCOOR, int(subtomo.getCoordinate3D().getX()*scaleFactor))
                    nRow.setValue(xmippLib.MDL_YCOOR, int(subtomo.getCoordinate3D().getY()*scaleFactor))
                    nRow.setValue(xmippLib.MDL_ZCOOR, int(subtomo.getCoordinate3D().getZ()*scaleFactor))
                    # Convert transform matrix to Euler Angles (rot, tilt, psi)
                    alignmentToRow(subtomo.getTransform(), nRow, ALIGN_3D)
                    nRow.addToMd(mdGeometry)
            fnGeometry = self._getExtraPath("geometry%d.xmd" % classId)
            mdGeometry.write(fnGeometry)

            if self.inputClasses.get().getSamplingRate() != tomo.getSamplingRate():
                factor = self.inputClasses.get().getSamplingRate()/tomo.getSamplingRate()
                args = "-i %s -o %s --scale %d" % (self._getExtraPath("reference%d.mrc" % classId),
                                                   self._getExtraPath("reference%d.mrc" % classId), factor)
                self.runJob('xmipp_transform_geometry', args)

            if self.paintingType.get() == 0:
                painting = 'copy'
            elif self.paintingType.get() == 1:
                painting = 'avg %d' % self.threshold.get()
            elif self.paintingType.get() == 2:
                painting = 'highlight %d' % self.constant.get()
            elif self.paintingType.get() == 3:
                painting = 'copy_binary %d' % self.threshold.get()

            tomogram = self._getExtraPath("tomogram_%d.mrc" % tomo.getObjId())
            args = " -i %s -o %s --geom %s --ref %s --method %s" % (tomogram, tomogram,
                                                                    self._getExtraPath("geometry%d.xmd" % classId),
                                                                    self._getExtraPath("reference%d.mrc" % classId),
                                                                    painting)
            self.runJob("xmipp_tomo_map_back", args)

    def createOutput(self):
        inputTomos = self.inputTomograms.get()
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(inputTomos)
        for j, _ in enumerate(inputTomos):
            tomo = Tomogram()
            tomo.setLocation(self._getExtraPath("tomogram_%d.mrc" % int(j+1)))
            outputTomos.append(tomo)
        self._defineOutputs(outputTomograms=outputTomos)
        self._defineSourceRelation(self.inputTomograms, outputTomos)
        self._defineSourceRelation(self.inputClasses, outputTomos)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        for subtomo in self.inputClasses.get().getFirstItem().iterItems():
            if not subtomo.hasCoordinate3D():
                validateMsgs.append('Please provide a class which contains subtomograms with 3D coordinates.')
                break
        return validateMsgs

    def _summary(self):
        summary = []
        summary.append("%d subtomogram reference mapped back %d times to original tomograms" %
                       (len(self.inputClasses.get()), len(self.inputClasses.get().getFirstItem())))
        return summary

    def _methods(self):
        methods = []
        methods.append("References from %d subtomogram classes mapped back %d times to original tomograms %s" %
                       (len(self.inputClasses.get()), len(self.inputClasses.get().getFirstItem()),
                        self.getObjectTag('inputTomograms')))
        return methods