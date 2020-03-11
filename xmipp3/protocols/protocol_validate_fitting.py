# -*- coding: utf-8 -*-
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

import numpy as np

from pyworkflow import VERSION_2_0
from pyworkflow.protocol.params import (PointerParam,
                                        IntParam, FileParam)
from pyworkflow.em.protocol.protocol_3d import ProtAnalysis3D
from pyworkflow.utils import getExt
from shutil import copyfile
from pyworkflow.em.data import VolumeMask
from pyworkflow.em.convert import Ccp4Header

VALIDATE_METHOD_URL = 'http://github.com/I2PC/scipion-em-xmipp/wiki/XmippProtValFit'
OUTPUT_PDBVOL_FILE = 'pdbVol'
OUTPUT_PDBMRC_FILE = 'pdb_volume.mrc'
BLOCRES_AVG_FILE = 'blocresAvg'
BLOCRES_HALF_FILE = 'blocresHalf'
RESTA_FILE = 'diferencia.vol'
PDB_VALUE_FILE = 'pdb_diferencia.pdb'
MASK_FILE_MRC = 'mask.mrc'
MASK_FILE = 'mask.vol' 


class XmippProtValFit(ProtAnalysis3D):
    """    
    The protocol evaluates the quality of the fitting.
    """
    _label = 'validate fitting'
    _lastUpdateVersion = VERSION_2_0
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.stepsExecutionMode = 1
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      help='Select a volume.')
        
#         form.addParam('inputPDB', PointerParam, pointerClass='PdbFile',
#                       label="Refined PDB", important=True, )
        form.addParam('inputPDB', FileParam,
                      label="File path", important=True,
                      help='Specify a path to desired PDB structure.')
        
        form.addParam('pdbMap', PointerParam, pointerClass='Volume',
                      label="Volume from PDB", allowsNull=True,
                      help='Volume created from the PDB.'
                           ' The volume should be aligned with the reconstruction map.'
                           ' If the volume is not entered,' 
                           ' it is automatically created from the PDB.')        

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask', 
                      allowsNull=True,
                      label="Soft Mask", 
                      help='The mask determines which points are specimen'
                      ' and which are not. If the mask is not passed,' 
                      ' the method creates an automatic mask from the PDB.')
        
        form.addParam('box', IntParam, default=20,
                      label="window size",
                      help='Kernel size for determining'
                      ' local resolution (pixels/voxels).')
        form.addParallelSection(threads=8, mpi=1)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {
                 OUTPUT_PDBVOL_FILE: self._getTmpPath('pdb_volume'),
                 OUTPUT_PDBMRC_FILE: self._getExtraPath('pdb_volume.mrc'),
                 BLOCRES_AVG_FILE: self._getExtraPath('blocres_avg.map'),
                 BLOCRES_HALF_FILE: self._getExtraPath('blocres_half.map'),
                 RESTA_FILE: self._getExtraPath('diferencia.vol'),
                 PDB_VALUE_FILE:  self._getExtraPath('pdb_diferencia.pdb'), 
                 MASK_FILE_MRC : self._getExtraPath('mask.mrc'),  
                 MASK_FILE: self._getExtraPath('mask.vol')               
                 }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
 
        self._createFilenameTemplates() 
        input = self._insertFunctionStep('convertInputStep')
        Id = []
        for i in range(2):
            Id.append(self._insertFunctionStep('runBlocresStep', i, prerequisites=[input]))   
        input1 =self._insertFunctionStep('substractBlocresStep',prerequisites=Id)    

        self._insertFunctionStep('assignPdbStep', prerequisites=[input1])     


    def convertInputStep(self):
        """ Read the input volume."""
        volume = self.inputVolume.get()
        
        """Read the Origin."""
        self.shifts =volume.getOrigin(force=True).getShifts()
#         print("kakakakakakakaka")
#         print("origen = %f %f %f") %(shifts[0], shifts[1],shifts[2])
#         
#         fileName = volume.getFileName()
#         newFileName = self._getExtraPath('prueba.mrc'),
#         Ccp4Header.fixFile(fileName, newFileName, shifts,
#                            volume.getSamplingRate(), Ccp4Header.ORIGIN)
        
        self.vol = volume.getFileName()
        self.half1, self.half2 = volume.getHalfMaps().split(',')

        extVol = getExt(self.vol)        
        if (extVol == '.mrc') or (extVol == '.map'):
            self.vol_xmipp = self.vol + ':mrc'  
            
        """Create map from PDB """         
        if self.pdbMap.hasValue():   
            pdbvolume = self.pdbMap.get()
            self.pdbvol = pdbvolume.getFileName()
            ext = getExt(self.pdbvol) 
            if (ext == '.mrc') or (ext == '.map'):                    
                copyfile(pdbvolume.getFileName(),self._getFileName(OUTPUT_PDBMRC_FILE))
            else:
                params = ' -i %s' % self.pdbvol 
                params += ' -o %s' % self._getFileName(OUTPUT_PDBMRC_FILE)   
                params += ' -t vol'            
                self.runJob('xmipp_image_convert', params)
                
                params = ' -i %s' % self._getFileName(OUTPUT_PDBMRC_FILE) 
                params += ' -s %f' % self.inputVolume.get().getSamplingRate()           
                self.runJob('xmipp_image_header', params)
      
        else:         
            """ Convert PDB to Map """           
            params = ' --centerPDB '
            params += ' -v 0 '        
            params += ' --sampling %f' % self.inputVolume.get().getSamplingRate()        
            params += ' --size %d' % self.inputVolume.get().getXDim()
            params += ' -i %s' % self.inputPDB.get()        
            params += ' -o %s' % self._getFileName(OUTPUT_PDBVOL_FILE)
            self.runJob('xmipp_volume_from_pdb', params)            
    
            """ Align pdbMap to reconstruction Map """  
              
            params = ' --i1 %s' % self.vol_xmipp        
            params += ' --i2 %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol'
            params += ' --local --apply'  
            params += ' %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol'                    
            self.runJob('xmipp_volume_align', params)  
            
            """ convert align vol to mrc format """
            
            params = ' -i %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol' 
            params += ' -o %s' % self._getFileName(OUTPUT_PDBMRC_FILE)   
            params += ' -t vol'            
            self.runJob('xmipp_image_convert', params)
            
            params = ' -i %s' % self._getFileName(OUTPUT_PDBMRC_FILE) 
            params += ' -s %f' % self.inputVolume.get().getSamplingRate()           
            self.runJob('xmipp_image_header', params)  
        
        
        """ Create a mask"""               
        if self.Mask.hasValue():
            self.maskFn = self.Mask.get().getFileName()
            extM = getExt(self.maskFn)        
            if (extM == '.mrc') or (extM == '.map'):
                self.mask_xmipp = self.maskFn + ':mrc' 
        else:
            
            self.maskFn = self._getFileName(MASK_FILE_MRC)
            self.mask_xmipp = self._getFileName(MASK_FILE) 
            
            if (not self.pdbMap.hasValue()):
                            
                params = ' -i %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol'          
                params += ' -o %s' % self.mask_xmipp
                params += ' --select below 0.02 --substitute binarize'                    
                self.runJob('xmipp_transform_threshold', params) 
                 
                params = ' -i %s' % self.mask_xmipp        
                params += ' -o %s' % self.mask_xmipp
                params += ' --binaryOperation dilation --size 3'                    
                self.runJob('xmipp_transform_morphology', params) 
                 
                """ convert mask.vol to mrc format """
             
                params = ' -i %s' % self.mask_xmipp  
                params += ' -o %s' % self.maskFn   
                params += ' -t vol'            
                self.runJob('xmipp_image_convert', params)
            
                params = ' -i %s' % self.maskFn
                params += ' -s %f' % self.inputVolume.get().getSamplingRate()           
                self.runJob('xmipp_image_header', params) 
                
            else:
                """ Convert PDB to Map """           
                params = ' --centerPDB '
                params += ' -v 0 '        
                params += ' --sampling %f' % self.inputVolume.get().getSamplingRate()        
                params += ' --size %d' % self.inputVolume.get().getXDim()
                params += ' -i %s' % self.inputPDB.get()        
                params += ' -o %s' % self._getFileName(OUTPUT_PDBVOL_FILE)
                self.runJob('xmipp_volume_from_pdb', params)            
        
                """ Align pdbMap to reconstruction Map """  
                  
                params = ' --i1 %s' % self.vol_xmipp        
                params += ' --i2 %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol'
                params += ' --local --apply'  
                params += ' %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol'                    
                self.runJob('xmipp_volume_align', params)  
                
                params = ' -i %s' % self._getFileName(OUTPUT_PDBVOL_FILE)+'.vol'          
                params += ' -o %s' % self.mask_xmipp
                params += ' --select below 0.02 --substitute binarize'                    
                self.runJob('xmipp_transform_threshold', params) 
                 
                params = ' -i %s' % self.mask_xmipp        
                params += ' -o %s' % self.mask_xmipp
                params += ' --binaryOperation dilation --size 3'                    
                self.runJob('xmipp_transform_morphology', params) 
                 
                """ convert mask.vol to mrc format """
             
                params = ' -i %s' % self.mask_xmipp  
                params += ' -o %s' % self.maskFn   
                params += ' -t vol'            
                self.runJob('xmipp_image_convert', params)
            
                params = ' -i %s' % self.maskFn
                params += ' -s %f' % self.inputVolume.get().getSamplingRate()           
                self.runJob('xmipp_image_header', params) 
                
            
#             self.maskFn = VolumeMask()
#             self.maskFn.setSamplingRate(self.inputVolume.get().getSamplingRate())
 

                        
                              
    def runBlocresStep(self, i):
        # Local import to prevent discovery errors
        import bsoft

        if (i==0):

            """ Calculate FSC map-PDB """

            params = ' -criterio FSC -nofill -smooth -pad 1 '
            params += ' -cutoff 0.67'
            params += ' -maxresolution 2 '
            params += ' -step 1 '
            params += ' -box %d ' % self.box.get()
            params += ' -sampling %f,%f,%f' % (self.inputVolume.get().getSamplingRate(),
                                               self.inputVolume.get().getSamplingRate(),
                                               self.inputVolume.get().getSamplingRate())
#             params += ' -origin %f,%f,%f' % ((self.shifts[0], self.shifts[1], self.shifts[2]))          
            params += ' -Mask %s' % self.maskFn
            params += ' %s  %s' % (self.vol, self._getFileName(OUTPUT_PDBMRC_FILE))
            params += ' %s' % self._getFileName(BLOCRES_AVG_FILE)

            self.runJob(bsoft.Plugin.getProgram('blocres'), params,
                        env=bsoft.Plugin.getEnviron())
        else:

            """ Calculate FSC half1-half2 """

            params = ' -criterio FSC -nofill -smooth -pad 1 '
            params += ' -cutoff 0.5'
            params += ' -maxresolution 2 '
            params += ' -step 1 '
            params += ' -box %d ' % self.box.get()
            params += ' -sampling %f,%f,%f' % (self.inputVolume.get().getSamplingRate(),
                                               self.inputVolume.get().getSamplingRate(),
                                               self.inputVolume.get().getSamplingRate()) 
#             params += ' -origin %f,%f,%f' % ((self.shifts[0], self.shifts[1], self.shifts[2])) 
            params += ' -Mask %s' % self.maskFn
            params += ' %s  %s' % (self.half1, self.half2)
            params += ' %s' % self._getFileName(BLOCRES_HALF_FILE)

            self.runJob(bsoft.Plugin.getProgram('blocres'), params,
                        env=bsoft.Plugin.getEnviron())

    def substractBlocresStep(self):
        
        params = ' -i %s' % self._getFileName(BLOCRES_AVG_FILE)+':mrc'  
        params += ' --minus %s' % self._getFileName(BLOCRES_HALF_FILE)+':mrc'     
        params += ' -o %s ' % self._getFileName(RESTA_FILE)       
        self.runJob('xmipp_image_operate', params)
        
    def assignPdbStep(self):
        
        shifts= self.inputVolume.get().getOrigin(force=True).getShifts()
        originX=shifts[0]/self.inputVolume.get().getSamplingRate()
        originY=shifts[1]/self.inputVolume.get().getSamplingRate()
        originZ=shifts[2]/self.inputVolume.get().getSamplingRate()   
        print('origen = %f %f %f' %(originX, originY, originZ))    
        params = ' --pdb %s ' % self.inputPDB.get()  
        params += ' --vol %s ' % self._getFileName(RESTA_FILE) 
        params += ' --mask %s ' % self.mask_xmipp         
        params += ' -o %s ' % self._getFileName(PDB_VALUE_FILE)    
        params += ' --sampling %f' % self.inputVolume.get().getSamplingRate()
        params += ' --origin %f %f %f' %(originX, originY, originZ)
        params += ' --radius 1' 
        self.runJob('xmipp_pdb_from_volume', params)       


    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'resolution_Volume'):
            messages.append(
                'Information about the method/article in ' + VALIDATE_METHOD_URL)
        return messages

    def _validate(self):
        errors = []
        if self.inputVolume.hasValue():
            #FIX CSVList:
            # volume.hasHalfMaps() does not work unless you call
            # getHalfMaps() or something else that triggers the CSVList.get()
            # that, populates the objValue. Just a print vol.getHalfMaps() will
            # change the behaviour of hasValue()
            # To review when migrating to Scipion3
            if not self.inputVolume.get().getHalfMaps():
                errors.append("Input Volume needs to have half maps. "
                "If you have imported the volume, be sure to import the half maps.")

        try:
            import bsoft
        except ImportError as e:
            errors.append("This protocol requires bsoft plugin to run.")
            
        return errors    

    def _citations(self):
        return ['Ramirez-Aportela 2020']

