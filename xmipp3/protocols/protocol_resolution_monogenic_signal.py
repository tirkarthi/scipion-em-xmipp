# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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
from pyworkflow.protocol.params import (PointerParam, BooleanParam, FloatParam,
                                        LEVEL_ADVANCED)
from pwem.protocols import ProtAnalysis3D
from pyworkflow.object import Float
from pwem.emlib.image import ImageHandler
from pyworkflow.utils import getExt
from pwem.objects import Volume
import pwem.emlib.metadata as md
from pwem.convert import Ccp4Header


MONORES_METHOD_URL = 'http://github.com/I2PC/scipion/wiki/XmippProtMonoRes'
OUTPUT_RESOLUTION_FILE = 'monoresResolutionMap.mrc'
OUTPUT_RESOLUTION_FILE_CHIMERA = 'monoresResolutionChimera.mrc'
OUTPUT_MASK_FILE = 'refinedMask.mrc'
FN_MEAN_VOL = 'meanvol'
METADATA_MASK_FILE = 'metadataresolutions'
FN_METADATA_HISTOGRAM = 'hist.xmd'
BINARY_MASK = 'binarymask'
FN_GAUSSIAN_MAP = 'gaussianfilter'


class XmippProtMonoRes(ProtAnalysis3D):
    """    
    Given a map the protocol assigns local resolutions to each voxel of the map.
    """
    _label = 'local MonoRes'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.min_res_init = Float() 
        self.max_res_init = Float()
       
    
    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('halfVolumes', BooleanParam, default=False,
                      label="Would you like to use half volumes?",
                      help='The noise estimation for determining the '
                      'local resolution is performed via half volumes.')

        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input Volume", important=True,
                      condition = 'not halfVolumes',
                      help='Select a volume for determining its '
                      'local resolution.')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Volume Half 1", important=True,
                      condition = 'halfVolumes', 
                      help='Select a volume for determining its '
                      'local resolution.')

        form.addParam('inputVolume2', PointerParam, pointerClass='Volume',
                      label="Volume Half 2", important=True,
                      condition='halfVolumes',
                      help='Select a second volume for determining a '
                      'local resolution.')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask', 
                      condition='(halfVolumes) or (not halfVolumes)',
                      allowsNull=True,
                      label="Binary Mask", 
                      help='The mask determines which points are specimen'
                      ' and which are not')

        form.addParam('maskExcl', PointerParam, pointerClass='VolumeMask',
		      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="Exclude Area", 
                      help='The mask determines the area of the protein to be'
		      'excluded in the estimation of the local resolution')

        group = form.addGroup('Extra parameters')
        line = group.addLine('Resolution Range (Å)',
                            help="If the user knows the range of resolutions or"
                                " only a range of frequencies needs to be analysed." 
                                "If Low is empty MonoRes will try to estimate the range. "
                                "it should be better if a range is provided")
        
        group.addParam('significance', FloatParam, default=0.95, 
                       expertLevel=LEVEL_ADVANCED,
                      label="Significance",
                      help='Relution is computed using hipothesis tests, '
                      'this value determines the significance of that test')
        
        group.addParam('maskthreshold', FloatParam, default=0.5, 
                       expertLevel=LEVEL_ADVANCED,
                      label="Mask threshold",
                      help='If the provided mask is not binary. Then, MonoRes'
                      'will try to binarize it. Bmask values below the threshold'
                      'will be change to 0 and above the thresthol will be 1')
        
        form.addParam('noiseonlyinhalves', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label="Use noise inside protein?", 
                      condition = 'halfVolumes',
                      help='(Yes recommended) the noise distribution will be estimated'
		      ' in the protein region (inside the mask) by means of the '
		      'difference of both half maps.')

        form.addParam('gaussianNoise', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Consider noise gaussian?", 
                      help='It assumens the noise in the map as gaussian. '
			'Note that this assumption might not be true, despite ,it is '
			' in general.')
        
        line.addParam('minRes', FloatParam, default=0, label='High')
        line.addParam('maxRes', FloatParam, allowsNull=True, label='Low')
        line.addParam('stepSize', FloatParam, allowsNull=True,
                      expertLevel=LEVEL_ADVANCED, label='Step')
        
        form.addParallelSection(threads = 4, mpi = 0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called """
        myDict = {
                 FN_MEAN_VOL: self._getExtraPath('mean_volume.mrc'),
                 METADATA_MASK_FILE: self._getExtraPath('mask_data.xmd'),
                 BINARY_MASK: self._getExtraPath('binarized_mask.mrc'),
                 FN_GAUSSIAN_MAP: self._getExtraPath('gaussianfilted.mrc')
                 }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
            # Convert input into xmipp Metadata format
        self._createFilenameTemplates() 
        self._insertFunctionStep('convertInputStep', )
        self._insertFunctionStep('resolutionMonogenicSignalStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep("createHistogram")

    def convertInputStep(self):
        """ Read the input volume.
        """
        
        self.micsFn = self._getPath()

        if self.halfVolumes:
            self.vol1Fn = self.inputVolume.get().getFileName()
            self.vol2Fn = self.inputVolume2.get().getFileName()
            extVol1 = getExt(self.vol1Fn)
            extVol2 = getExt(self.vol2Fn)
            if (extVol1 == '.mrc') or (extVol1 == '.map'):
                self.vol1Fn = self.vol1Fn + ':mrc'
            if (extVol2 == '.mrc') or (extVol2 == '.map'):
                self.vol2Fn = self.vol2Fn + ':mrc'
                
            if not self.Mask.hasValue():
                self.ifNomask(self.vol1Fn)
            else:
                self.maskFn = self.Mask.get().getFileName()

            self.inputVolumes.set(None)
        else:
            self.vol0Fn = self.inputVolumes.get().getFileName()
            extVol0 = getExt(self.vol0Fn)
            if (extVol0 == '.mrc') or (extVol0 == '.map'):
                self.vol0Fn = self.vol0Fn + ':mrc'
                
            if not self.Mask.hasValue():
                self.ifNomask(self.vol0Fn)
            else:
                self.maskFn = self.Mask.get().getFileName()
            
            self.inputVolume.set(None)
            self.inputVolume2.set(None)


        extMask = getExt(self.maskFn)
        if (extMask == '.mrc') or (extMask == '.map'):
            self.maskFn = self.maskFn + ':mrc'
            
        if self.Mask.hasValue():
            params = ' -i %s' % self.maskFn
            params += ' -o %s' % self._getFileName(BINARY_MASK)
            params += ' --select below %f' % self.maskthreshold.get()
            params += ' --substitute binarize'
             
            self.runJob('xmipp_transform_threshold', params)

    def ifNomask(self, fnVol):
        if self.halfVolumes:
            xdim, _ydim, _zdim = self.inputVolume.get().getDim()
            params = ' -i %s' % fnVol
        else:
            xdim, _ydim, _zdim = self.inputVolumes.get().getDim()
            params = ' -i %s' % fnVol
        params += ' -o %s' % self._getFileName(FN_GAUSSIAN_MAP)
        setsize = 0.02*xdim
        params += ' --fourier real_gaussian %f' % (setsize)
     
        self.runJob('xmipp_transform_filter', params)
        img = ImageHandler().read(self._getFileName(FN_GAUSSIAN_MAP))
        imgData = img.getData()
        max_val = np.amax(imgData)*0.05
         
        params = ' -i %s' % self._getFileName(FN_GAUSSIAN_MAP)
        params += ' --select below %f' % max_val
        params += ' --substitute binarize'
        params += ' -o %s' % self._getFileName(BINARY_MASK)
     
        self.runJob('xmipp_transform_threshold', params)
        
        self.maskFn = self._getFileName(BINARY_MASK)


    def resolutionMonogenicSignalStep(self):
        if self.stepSize.hasValue():
            freq_step = self.stepSize.get()
        else:
            freq_step = 0.25

        if self.halfVolumes:
            params = ' --vol %s' % self.vol1Fn
            params += ' --vol2 %s' % self.vol2Fn
            params += ' --sampling_rate %f' % self.inputVolume.get().getSamplingRate()
            if (self.noiseonlyinhalves.get() is True):
                params += ' --noiseonlyinhalves'
        else:
            params = ' --vol %s' % self.vol0Fn
            params += ' --sampling_rate %f' % self.inputVolumes.get().getSamplingRate()
        
        params += ' --mask %s' % self._getFileName(BINARY_MASK)
        
        if self.maskExcl.hasValue():
            params += ' --maskExcl %s' % self.maskExcl.get().getFileName()

        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --step %f' % freq_step
        params += ' -o %s' % self._getExtraPath()
        if self.gaussianNoise.get() is True:
            params += ' --gaussian'
        params += ' --significance %f' % self.significance.get()
        params += ' --threads %i' % self.numberOfThreads.get()  

        self.runJob('xmipp_resolution_monogenic_signal', params)



    def createHistogram(self):

        M = float(self.max_res_init)
        m = float(self.min_res_init)
        range_res = round((M - m)*4.0)

        params = ' -i %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --mask binary_file %s' % self._getExtraPath(OUTPUT_MASK_FILE)
        params += ' --steps %f' % (range_res)
        params += ' --range %f %f' % (self.min_res_init, self.max_res_init)
        params += ' -o %s' % self._getExtraPath(FN_METADATA_HISTOGRAM)

        self.runJob('xmipp_image_histogram', params)
        
        
    def readMetaDataOutput(self):
        mData = md.MetaData(self._getFileName(METADATA_MASK_FILE))
        NvoxelsOriginalMask = float(mData.getValue(md.MDL_COUNT, mData.firstObject()))
        NvoxelsOutputMask = float(mData.getValue(md.MDL_COUNT2, mData.firstObject()))
        nvox = int(round(
                ((NvoxelsOriginalMask-NvoxelsOutputMask)/NvoxelsOriginalMask)*100))
        return nvox

    def getMinMax(self, imageFile):
        img = ImageHandler().read(imageFile)
        imgData = img.getData()
        min_res = round(np.amin(imgData) * 100) / 100
        max_res = round(np.amax(imgData) * 100) / 100
        return min_res, max_res

    def createOutputStep(self):
        volume=Volume()
        volume.setFileName(self._getExtraPath(OUTPUT_RESOLUTION_FILE))
        if (self.halfVolumes):
            volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
            volume.setOrigin(self.inputVolume.get().getOrigin(True))
            self._defineOutputs(resolution_Volume=volume)
            self._defineSourceRelation(self.inputVolume, volume)
            inputVolumeFileName = self.inputVolume.get().getFileName()
        else:
            volume.setSamplingRate(self.inputVolumes.get().getSamplingRate())
            volume.setOrigin(self.inputVolumes.get().getOrigin(True))
            self._defineOutputs(resolution_Volume=volume)
            self._defineSourceRelation(self.inputVolumes, volume)
            inputVolumeFileName = self.inputVolumes.get().getFileName()

        #Setting the min max for the summary
        imageFile = self._getExtraPath(OUTPUT_RESOLUTION_FILE_CHIMERA)
        min_, max_ = self.getMinMax(imageFile)
        self.min_res_init.set(round(min_*100)/100)
        self.max_res_init.set(round(max_*100)/100)
        self._store(self.min_res_init)
        self._store(self.max_res_init)
            

        # fill ccp4 header so we can transfer the origin to the
        # viewer (otherwise since imageFile is not a Scipion object
        # no sampling/origin information can be transfered

        #Ccp4Header(imageFile).copyCCP4Header(
        #    inputVolumeFileName, volume.getShiftsFromOrigin(),
        #    volume.getSamplingRate(), originField=Ccp4Header.START)

        # also update the output volume header. This is not needed
        # since sampling and origin is in the database but
        # it may be usefull if other programs -outside scipion-
        # require  these data.

        #Ccp4Header(volume.getFileName()).copyCCP4Header(
        #     inputVolumeFileName, volume.getShiftsFromOrigin(),
        #     volume.getSamplingRate(), originField=Ccp4Header.START)



    # --------------------------- INFO functions ------------------------------

    def _methods(self):
        messages = []
        if hasattr(self, 'resolution_Volume'):
            messages.append(
                'Information about the method/article in ' + MONORES_METHOD_URL)
        return messages
    
    def _summary(self):
        summary = []
        summary.append("Highest resolution %.2f Å,   "
                       "Lowest resolution %.2f Å. \n" % (self.min_res_init,
                                                         self.max_res_init))
        return summary

    def _citations(self):
        return ['Vilas2018']

