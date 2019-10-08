# **************************************************************************
# *
# * Authors:     Erney Ramirez Aportela (eramirez@cnb.csic.es)
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

from pyworkflow import VERSION_2_0
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils import getExt, copyFile
from pyworkflow.em.constants import ALIGN_PROJ
from xmipp3.protocols.protocol_preprocess.protocol_process import XmippProcessParticles
from xmipp3.utils import validateDLtoolkit
from xmipp3.convert import writeSetOfParticles, locationToXmipp
import xmippLib


class XmippProtUmbral(XmippProcessParticles):
    """ Given a set of particles, the protocol resize the particle from Mask center of mass.
    """
    _label = 'particle umbral'
    _lastUpdateVersion = VERSION_2_0
    
    # --------------------------- DEFINE param functions ----------------------
#     def _defineProcessParams(self, form):
#         form.addParam('inputVolume', PointerParam, pointerClass='Volume',
#                       label="Input Map", important=True,
#                       help='Select a volume to calculate the center of mass.')             

    # --------------------------- INSERT steps functions -----------------------
    def _insertProcessStep(self):
        # Convert input into xmipp Metadata format
        self._insertFunctionStep('umbralStep')

    def convertInputStep(self):
        # We do not need conversion since we are doing all inside python
        pass

        #XmippProcessParticles.convertInputStep(self)

        #writeSetOfParticles(self.inputParticles.get(),
        #                    self._getExtraPath('input_particles.xmd'))
        #writeSetOfParticles(self.inputParticles.get(), self.inputFn,
        #                    alignType=ALIGN_PROJ)

        #copyFile(self._getTmpPath("input_particles.xmd"),self._getExtraPath("input.xmd"))


    def umbralStep(self):

        particles = self.inputParticles.get()

        x = y = particles.getXDim()
        size = particles.getSize()
        xmippLib.createEmptyFile(self.outputStk, x, y, 1, size)

        for particle in self.inputParticles.get():

            index, filename = particle.getLocation()
            fn = locationToXmipp(index, filename)

        #md = xmippLib.MetaData(self._getExtraPath("input.xmd"))
        #mdout=md
        #InputIm = xmippLib.Image()       
        #inputI=md.getValue(xmippLib.MDL_IMAGE, xmippLib.MDL_ITEM_ID)

        #InputIm = xmipp3.Image(self.inputFn)
            inputIm=xmippLib.Image(fn)
            print("Filename %s " % fn)
            I = inputIm.getData()
            OutI = I*0
            Ydim, Xdim = I.shape

            for y in range(1, Ydim-1):
                for x in range(1, Xdim-1):
                    ct=0
                    if I[y,x]<0.0001:
                        ct=0
                        if I[y+1,x]<0.0001:
                            ct+=1
                        if I[y-1,x]<0.0001:
                            ct+=1
                        if I[y+1,x+1]<0.0001:
                            ct+=1
                        if I[y-1,x+1]<0.0001:
                            ct+=1
                        if I[y+1,x-1]<0.0001:
                            ct+=1
                        if I[y-1,x-1]<0.0001:
                            ct+=1
                        if I[y,x+1]<0.0001:
                            ct+=1
                        if I[y,x-1]<0.0001:
                            ct+=1
                        if ct>=6:
                            OutI[y,x]=0 
                        else:
                            OutI[y,x]=I[y,x] 
 
                             
                    else:
                        OutI[y,x]=I[y,x]                               


#                     if I[y,x]<0.000001:
#                         OutI[y,x]=0.0
#                     else:
#                         OutI[y,x]=I[y,x]

            inputIm.setData(OutI)
            inputIm.write((index, self.outputStk))
        print(self.outputStk)
        print("kakoraaaaaaaa2")
        #md.setValue(xmippLib.MDL_IMAGE, inputIm, xmippLib.MDL_ITEM_ID)

    def createOutputStep(self):
        inputSet = self.inputParticles.get()
        # outputSet could be SetOfParticles, SetOfAverages or any future sub-class of SetOfParticles

        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputSet)

        outputSet.copyItems(inputSet,
                            updateItemCallback=self._updateItem)

        self._defineOutputs(outputParticles=outputSet)
        self._defineTransformRelation(inputSet, outputSet)

    def _updateItem(self, item, row):
        """ Update just the filename
        """
        item.setFileName(self.outputStk)
    # --------------------------- INFO functions ------------------------------
    def _methods(self):
        messages = []
        if hasattr(self, 'boxsize'):
            messages.append('Estimated box size: %s pixels' % self.boxsize)
        return messages
    
    def _summary(self):
        messages = []
        if hasattr(self, 'boxsize'):
            messages.append('Estimated box size: %s pixels' % self.boxsize)
            messages.append("Open 'Analyze results' to see it in context. "
                            "(The displayed coordinate is just to show the size)")
        return messages

    def _citations(self):
        return ['']

