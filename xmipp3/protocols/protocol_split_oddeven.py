# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import (PointerParam, EnumParam)
from pyworkflow.em import ProtCTFMicrographs
from pyworkflow.object import Float
from pyworkflow.em import ImageHandler
from pyworkflow.utils import getExt
from pyworkflow.em.data import Volume
import pyworkflow.em.metadata as md


class XmippProtSplitOddEven(ProtCTFMicrographs):
    """    
    Given a set of micrographs, they can be split in two sets like odd and even.
    """
    _label = 'split odd even'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtCTFMicrographs.__init__(self, **args)

       
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMicrographs', PointerParam,
                      pointerClass='SetOfMicrographs',
                      label="Micrographs",
                      help="Select a set of micrographs to be split.")
        
        form.addParam('typeSplit', EnumParam, label='split type', 
                      choices=['Chess', 'Rows', 'Columns'],
                      default=0, help="Type of splitting: The generated migragraph"
                      "will be the result of interpolating rows, columns or a chess"
                      "board pattern")
        
#         form.addParallelSection(threads = 1, mpi = 1)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
            # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('splitingStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        """ Read the input volume.
        """
        print('Starting...')
 
 
    def splitingStep(self):
        # Number of frequencies
        print('splitting step')
        numberofmic = self.inputMicrographs.get()
        print(numberofmic)
        
        if self.typeSplit.get() == 0:
            splittype = 'chess'
        if self.typeSplit.get() == 1:
            splittype = 'rows'
        if self.typeSplit.get() == 2:
            splittype = 'columns'
            
        for mic in self.inputMicrographs.get():
            # TODO: xmipp program  
            print(mic.getFileName()) 
            params = ' -i %s' % mic.getFileName()
            params += ' --odd %s' % self._getExtraPath('aa')
            params += ' --even %s' % self._getExtraPath('aa')
            params += ' --type %s' % splittype

            print(params)
#         self.runJob('xmipp_resolution_monogenic_signal', params)


    def createOutputStep(self):
        print('output')
#         volume=Volume()
#         volume.setFileName(self._getFileName(OUTPUT_RESOLUTION_FILE))
#         if (self.halfVolumes):
#             volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
#             self._defineOutputs(resolution_Volume=volume)
#             self._defineSourceRelation(self.inputVolume, volume)
#         else: 
#             volume.setSamplingRate(self.inputVolumes.get().getSamplingRate())
#             self._defineOutputs(resolution_Volume=volume)
#             self._defineSourceRelation(self.inputVolumes, volume)
            

    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'resolution_Volume'):
            messages.append(
                'Information about the method ')
        return messages
    
    def _summary(self):
        summary = []
        summary.append("Highest resolution %.2f Å,   "
                       "Lowest resolution %.2f Å. \n" % (1,5))
        return summary


