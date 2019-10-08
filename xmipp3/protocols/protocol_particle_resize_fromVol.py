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
from xmipp3.convert import writeSetOfParticles


class XmippProtParticleResizeFromVol(XmippProcessParticles):
    """ Given a set of particles, the protocol resize the particle from Mask center of mass.
    """
    _label = 'particle resize'
    _lastUpdateVersion = VERSION_2_0
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Map", important=True,
                      help='Select a volume to calculate the center of mass.')   
        form.addParam('size', IntParam, 
                      label="dim", important=True,
                      help='dimensions of window to apply to images. '
                           'It must be an integer value.')             

    # --------------------------- INSERT steps functions -----------------------
    def _insertProcessStep(self):
        # Convert input into xmipp Metadata format
        self._insertFunctionStep('resizeStep')

    def convertInputStep(self):

        #XmippProcessParticles.convertInputStep(self)

        #writeSetOfParticles(self.inputParticles.get(),
        #                    self._getExtraPath('input_particles.xmd'))
        writeSetOfParticles(self.inputParticles.get(), self.inputFn,
                            alignType=ALIGN_PROJ)

        copyFile(self._getTmpPath("input_particles.xmd"),self._getExtraPath("input.xmd"))


        self.volFn = self.inputVolume.get().getFileName()
        extVol = getExt(self.volFn)       
        if (extVol == '.mrc') or (extVol == '.map'):
            self.volFn = self.volFn + ':mrc'          
           

    def resizeStep(self):

        params  = ' -i %s' % self.inputFn
        params += ' --vol %s' % self.volFn
        params += ' --boxSize %s' % self.size
        params += ' -o %s' % self.outputStk

        self.runJob('xmipp_transform_window_from_volume', params)
        


    # def createOutputStep(self):
    #     """ The output is just an Integer. Other protocols can use it in those
    #         IntParam if it has set allowsPointer=True
    #     """
    #     particlesSet = self._createSetOfParticles()
    #     particlesSet.copyInfo(self.inputParticles.get())
    #     inputMd = self._getExtraPath('particle_resize.xmd')
    #     particlesSet.copyItems(self.inputParticles.get(),
    #                            updateItemCallback=self._updateParticle,
    #                            itemDataIterator=md.iterRows(inputMd))
    #     self._defineOutputs(outputParticles=particlesSet)
    #     self._defineSourceRelation(self.inputParticles.get(), particlesSet)
    #

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

    def _validate(self):
        return validateDLtoolkit(model=[('boxsize', 'weights.hdf5'),
                                        ('boxsize', 'feature_scaler.pkl')])
