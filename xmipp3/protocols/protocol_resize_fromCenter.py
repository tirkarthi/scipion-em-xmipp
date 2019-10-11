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
from xmipp3.protocols.protocol_preprocess.protocol_process import XmippProcessParticles, XmippProcessVolumes
from xmipp3.utils import validateDLtoolkit
from xmipp3.convert import writeSetOfParticles

class XmippWindowHelper():
    """ Common features to change dimensions of either SetOfParticles
    or Volume objects.
    """
    
        #--------------------------- DEFINE param functions --------------------------------------------
    @classmethod
    def _defineProcessParams(cls, protocol, form):
        # Common parameters
        form.addParam('inputMask', PointerParam, pointerClass='VolumeMask',
                      label="Mask", important=True,
                      help='Select a mask to calculate the center of mass.')   
        form.addParam('size', IntParam, 
                      label="window resize", important=True,
                      help='dimensions of window to apply to images. '
                           'It must be an integer value.') 
        
    #--------------------------- INSERT steps functions ------------------------
    @classmethod      
    def _insertProcessStep(cls, protocol):
        cls._checkExtension(protocol)
        
        args = protocol._resizeArgs()        
        protocol._insertFunctionStep('resizeStep', args)
        
    @staticmethod           
    def _checkExtension(protocol):

        volFn = protocol.inputMask.get().getFileName()
        extVol = getExt(volFn)       
        if (extVol == '.mrc') or (extVol == '.map'):
            volFn = volFn + ':mrc'   
        protocol.volFn=volFn
            
    #--------------------------- STEP functions --------------------------------   
    @staticmethod
    def resizeStep(protocol, args):
        protocol.runJob("xmipp_transform_window_from_volume", args)
        
   

class XmippProtParticleResizeFromVol(XmippProcessParticles):
    """ Given a set of particles, the protocol resize the particle from Mask center of mass.
    """
    _label = 'cut particle from center'
    _lastUpdateVersion = VERSION_2_0
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):
        XmippWindowHelper._defineProcessParams(self, form)
          

    # --------------------------- INSERT steps functions -----------------------    
    def _insertProcessStep(self):
        XmippWindowHelper._insertProcessStep(self)       
        

    def convertInputStep(self):
        
        XmippWindowHelper._checkExtension(self)

        writeSetOfParticles(self.inputParticles.get(), self.inputFn,
                            alignType=ALIGN_PROJ)

        copyFile(self._getTmpPath("input_particles.xmd"),self._getExtraPath("input.xmd")) 
                      
           

    def resizeStep(self, args):
        XmippWindowHelper.resizeStep(self, self._resizeArgs())
       
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeArgs(self):
            return "-i %s -o %s --mask %s --boxSize %s" % (self.inputFn, self.outputStk, self.volFn, self.size)
  

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



class XmippProtVolumeResizeFromVol(XmippProcessVolumes):
    """ Given a set of particles, the protocol resize the particle from Mask center of mass.
    """
    _label = 'cut volume from center'
    _lastUpdateVersion = VERSION_2_0
    
    def __init__(self, **kwargs):
        XmippProcessVolumes.__init__(self, **kwargs)    
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):
        XmippWindowHelper._defineProcessParams(self, form)
          

    # --------------------------- INSERT steps functions -----------------------    
    def _insertProcessStep(self):
        XmippWindowHelper._insertProcessStep(self)       
        

    def convertInputStep(self):
        
        XmippWindowHelper._checkExtension(self)
                 
          

    def resizeStep(self, args):
        XmippWindowHelper.resizeStep(self, self._resizeArgs())
       
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _resizeArgs(self):
            return "-i %s -o %s --mask %s --boxSize %s" % (self.inputFn, self.outputStk, self.volFn, self.size)
  

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
