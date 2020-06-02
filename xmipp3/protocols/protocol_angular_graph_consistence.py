# **************************************************************************
# *
# * Authors:         Jeison Mendez (jmendez@utp.edu.co) (2020)
# *
# * Universidad Nacional Autónoma de México, UNAM
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

from os.path import join, isfile
from shutil import copyfile
import os

from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        STEPS_PARALLEL,
                                        StringParam, BooleanParam, IntParam,
                                        LEVEL_ADVANCED, USE_GPU, GPU_LIST)

from pyworkflow.utils.path import moveFile, makePath, cleanPattern
from pyworkflow.gui.plotter import Plotter

from pwem.objects import Volume
from pwem import emlib
import pwem.emlib.metadata as md
from pwem.protocols import ProtAnalysis3D
from pwem.constants import ALIGN_PROJ

from xmipp3.convert import writeSetOfParticles, writeSetOfVolumes, \
    getImageLocation, setXmippAttributes, createItemMatrix



class XmippProtAngularGraphConsistence(ProtAnalysis3D):
    """
    Performs soft alignment validation of a set of particles previously aligned
    confronting them using Graph filtered correlations representation. This
    protocol produces an histogram with two groups of particles.
    """
    _label = 'angular graph consistence'

    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input volume",
                      help='Select the input volume(s).')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignment',
                      label="Input particles", important=True,
                      help='Select the input projection images.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group",
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp')
        form.addParam('angularSampling', FloatParam, default=3,
                      expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling (degrees)",
                      help='Angular distance (in degrees) between neighboring projection points ')
        form.addParam('doNotUseWeights', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Do not use the weights",
                      help='Do not use the weights in the clustering calculation')
        form.addParam('minTilt', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Minimum allowed tilt angle",
                      help='Tilts below this value will not be considered for the alignment')
        form.addParam('maxTilt', FloatParam, default=180,
                      expertLevel=LEVEL_ADVANCED,
                      label="Maximum allowed tilt angle without mirror check",
                      help='Tilts above this value will not be considered for the alignment without mirror check')

        form.addSection(label='Preprocess')
        form.addParam('doWiener', BooleanParam, default='True',
                      label="CTF correction",
                      help='Perform CTF correction by Wiener filtering.')
        form.addParam('isIsotropic', BooleanParam, default='True',
                      label="Isotropic Correction", condition='doWiener',
                      help='If true, Consider that there is not astigmatism and then it is performed an isotropic correction.')
        form.addParam('padding_factor', IntParam, default=2,
                      expertLevel=LEVEL_ADVANCED,
                      label="Padding factor", condition='doWiener',
                      help='Padding factor for Wiener correction ')
        form.addParam('wiener_constant', FloatParam, default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      label="Wiener constant", condition='doWiener',
                      help=' Wiener-filter constant (if < 0: use FREALIGN default)')
        form.addParam('correctEnvelope', BooleanParam, default='False',
                      expertLevel=LEVEL_ADVANCED,
                      label="Correct for CTF envelope", condition='doWiener',
                      help=' Only in cases where the envelope is well estimated correct for it')
        form.addParam('targetResolution', FloatParam, default=8,
                      label='Target resolution (A)',
                      help='Low pass filter the particles to this resolution. This usually helps a lot obtaining good alignment. You should have a good'
                           ' reason to modify this value outside the range  [8-10] A')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        convertId = self._insertFunctionStep('convertInputStep',
                                             self.inputParticles.get().getObjId())
        sym = self.symmetryGroup.get() 
        for i, vol in enumerate(self._iterInputVols()):
            volName = getImageLocation(vol)
            volDir = self._getVolDir(i + 1)
            pmStepId = self._insertFunctionStep('projectionLibraryStep',
                                                volName, volDir,
                                                self.angularSampling.get(),
                                                prerequisites=[convertId])
            sigStedId = self._insertFunctionStep('assignmentStep', 
                                                volName, volDir, sym,
                                                prerequisites=[pmStepId])
        self._insertFunctionStep('createOutputStep', prerequisites=[sigStedId])

    def assignmentStep(self, volName, volDir, sym):
        fnGallery = (volDir + '/gallery.doc')
        angStep = self.angularSampling.get();
        volName = self._getExtraPath("volume.vol")    
            
        newTs, newXdim = self._getModifiedSizeAndSampling() 
                
        args = ' -i %s' % self._getExtraPath('scaled_particles.xmd')
        args += ' -o %s' % self._getExtraPath('anglesOutput.xmd')
        args += ' -ref %s' % fnGallery
        args += ' -sampling %f' % newTs
        args += ' -odir %s' % volDir
        args += ' -angleStep %f' % angStep
        args += ' --Nsimultaneous %f' % self.numberOfMpi.get() * self.numberOfThreads.get()
        args += ' --sym %s' % sym
        args += ' --useForValidation'
        args += ' --refVol %s' % volName         
        self.runJob('xmipp_angular_assignment_mag', args, numberOfMpi=self.numberOfMpi.get() * self.numberOfThreads.get())

    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file.
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """

        writeSetOfParticles(self.inputParticles.get(),
                            self._getPath('input_particles.xmd'))

        if self.doWiener.get():
            params = '  -i %s' % self._getPath('input_particles.xmd')
            params += '  -o %s' % self._getExtraPath(
                'corrected_ctf_particles.stk')
            params += '  --save_metadata_stack %s' % self._getExtraPath(
                'corrected_ctf_particles.xmd')
            params += '  --pad %s' % self.padding_factor.get()
            params += '  --wc %s' % self.wiener_constant.get()
            params += '  --sampling_rate %s' % self.inputParticles.get().getSamplingRate()

            if self.inputParticles.get().isPhaseFlipped():
                params += '  --phase_flipped '

            if self.correctEnvelope:
                params += '  --correct_envelope '

            nproc = self.numberOfMpi.get()
            nT = self.numberOfThreads.get()

            self.runJob('xmipp_ctf_correct_wiener2d',
                        params)

        newTs, newXdim = self._getModifiedSizeAndSampling()

        if self.doWiener.get():
            params =  '  -i %s' % self._getExtraPath('corrected_ctf_particles.xmd')
        else :
            params =  '  -i %s' % self._getPath('input_particles.xmd')
            
        params +=  '  -o %s' % self._getExtraPath('scaled_particles.stk')
        params +=  '  --save_metadata_stack %s' % self._getExtraPath('scaled_particles.xmd')
        params +=  '  --fourier %d' % newXdim
        
        self.runJob('xmipp_image_resize',params)
        
        from pwem.emlib.image import ImageHandler
        img = ImageHandler()
        img.convert(self.inputVolumes.get(), self._getExtraPath("volume.vol"))
        Xdim = self.inputVolumes.get().getDim()[0]
        if Xdim != newXdim:
            self.runJob("xmipp_image_resize", "-i %s --dim %d" % \
                        (self._getExtraPath("volume.vol"),
                         newXdim), numberOfMpi=1)


    def _getModifiedSizeAndSampling(self):
        Xdim = self.inputParticles.get().getDimensions()[0]
        Ts = self.inputParticles.get().getSamplingRate()
        newTs = self.targetResolution.get() * 0.4
        newTs = max(Ts, newTs)
        newXdim = Xdim * Ts / newTs
        return newTs, newXdim

    def projectionLibraryStep(self, volName, volDir, angularSampling):
        # Generate projections from this reconstruction
        nproc = self.numberOfMpi.get()
        nT = self.numberOfThreads.get()
        volName = self._getExtraPath("volume.vol")
        makePath(volDir)
        fnGallery = (volDir + '/gallery.stk')
        params = '-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %f --experimental_images %s --max_tilt_angle %f --min_tilt_angle %f' \
                 % (
                 volName, fnGallery, angularSampling, self.symmetryGroup.get(),
                 -1, self._getExtraPath('scaled_particles.xmd'),
                 self.maxTilt.get(), self.minTilt.get())

        self.runJob("xmipp_angular_project_library", params, numberOfMpi=nproc,
                    numberOfThreads=nT)


    def createOutputStep(self):
        fnAngles = self._getExtraPath('anglesOutput.xmd')
        imgSet = self.inputParticles.get()
        imgSetOut = self._createSetOfParticles()
        imgSetOut.copyInfo(imgSet)
        imgSetOut.setAlignmentProj()
        imgSetOut.copyItems(imgSet,
                            updateItemCallback=self._updateParticle,
                            itemDataIterator=md.iterRows(fnAngles,
                                                         sortByLabel=md.MDL_ITEM_ID) )
        self._defineOutputs(outputParticles=imgSetOut)
        self._defineSourceRelation(self.inputParticles, imgSetOut)
        self.createPlot2D(fnAngles)
        cleanPattern(self._getExtraPath("scaled_particles.*"))
        cleanPattern(self._getExtraPath("corrected_ctf_particles.*"))
        cleanPattern(self._getExtraPath("volume.vol"))

    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_ANGLE_ROT, md.MDL_ANGLE_TILT,
                           md.MDL_ANGLE_PSI, md.MDL_SHIFT_X, md.MDL_SHIFT_Y,
                           md.MDL_MAXCC, md.MDL_WEIGHT,
                           md.MDL_ANGULAR_GRAPHCONSISTENCE)
        createItemMatrix(item, row, align=ALIGN_PROJ)
        
    def createPlot2D(self, fnAngles):
        mdParticles = emlib.MetaData(fnAngles)
        
        ccList = mdParticles.getColumnValues(emlib.MDL_MAXCC)
        graphConsistList = mdParticles.getColumnValues(emlib.MDL_ANGULAR_GRAPHCONSISTENCE)
        
        # threshold
        th_gsp = self.otsu(graphConsistList)
        th_cc = self.otsu(ccList)
        
        cont_gsp = 0        
        cont_cc = 0
        total = 0
        for objId in mdParticles:
            total += 1
            if ( mdParticles.getValue(emlib.MDL_ANGULAR_GRAPHCONSISTENCE, objId) > th_gsp ):
                cont_gsp += 1
            if ( mdParticles.getValue(emlib.MDL_MAXCC, objId) > th_cc ):
                cont_cc += 1
        # percentage of images in reliable assignment zone
        p_gsp = cont_gsp/total * 100
        p_cc = cont_cc/total * 100
        
        parameterFile = self._getExtraPath('parameter.txt')
        fh = open(parameterFile, "w")
        fh.write("%.2f" % p_gsp)
        fh.close()
                
        figurePath = self._getExtraPath('graphConsistenceHistogram.png')
        figureSize = (8, 6)
        plotter = Plotter(*figureSize)
        figure = plotter.getFigure()
        ax = figure.add_subplot(111)
        ax.set_title('Histogram - Soft alignment validation')
        ax.set_xlabel('Modified cross-correlation based on GSP')
        ax.set_ylabel('Num. of images')
        lb = '%.2f'%p_gsp
        lb += r'$\%$ images are in the reliable assignment zone'
        # histogram
        ax.hist(graphConsistList, bins=50, label=lb)
        ax.legend()
        ax.autoscale_view(True, True, True)
        plotter.savefig(figurePath)
        plotter.show()
        
        # histogram of angular_assignment_mag
        figurePath2 = self._getExtraPath('ccMaxHistogram.png')
        figureSize2 = (8, 6)
        plotter2 = Plotter(*figureSize2)
        figure2 = plotter2.getFigure()
        ax2 = figure2.add_subplot(111)
        ax2.set_title('Histogram - Soft alignment validation')
        ax2.set_xlabel('Modified cross-correlation based on GSP')
        ax2.set_ylabel('Num. of images')
        lb = '%.2f'%p_cc
        lb += r'$\%$ images are in the reliable assignment zone'
        # histogram
        ax2.hist(ccList, bins=50, label=lb)
        ax2.legend()
        ax2.autoscale_view(True, True, True)
        plotter2.savefig(figurePath2)
        
        return plotter
     
    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolumes.get() and not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')
        return validateMsgs

    def _summary(self):
        summary = [
            "Input particles:  %s" % self.inputParticles.get().getNameId()]
        summary.append("-----------------")
        
        if (not hasattr(self, 'outputParticles')):
            summary.append("Output not ready yet.")
        else:
            parameterFile = self._getExtraPath('parameter.txt')
            fh = open(parameterFile, "r")
            val = fh.readline()
            fh.close()
            text = 'After validation, a %s' % val
            text += '%'
            text += ' of images are within the reliable assignment zone'
            summary.append(text)
            figPath = self._getExtraPath('graphConsistenceHistogram.png')
            summary.append('histogram plot in %s' % figPath)
        return summary

    def _methods(self):
        messages = []
        if (hasattr(self, 'outputParticles')):
            messages.append('correlation values of previous alignment process'
                            ' are modified according to the spherical distance from'
                            ' the assigned direction to a soft and high-valued correlation zone'
                            ' of neighboring projection directions')
        return messages

    def _citations(self):
        return ['xxxxx']

    # --------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getTmpPath('vol%03d' % volIndex)

    def _iterInputVols(self):
        """ In this function we will encapsulate the logic
        to iterate through the input volumes.
        This give the flexibility of having Volumes, SetOfVolumes or
        a combination of them as input and the protocol code
        remain the same.
        """
        inputVols = self.inputVolumes.get()

        if isinstance(inputVols, Volume):
            yield inputVols
        else:
            for vol in inputVols:
                yield vol
            
    def otsu(self, ccList):
        # method from golden highres
        import numpy as np
        
        cc_number = len(ccList)
        mean_weigth = 1.0 / cc_number
        his, bins = np.histogram(ccList, int((max(ccList)-min(ccList))/0.01))
        final_thresh = -1
        final_value = -1
        cc_arr = np.arange(min(ccList),max(ccList),0.01)
        if len(cc_arr)==(len(his)+1):
            cc_arr=cc_arr[:-1]
        for t, j in enumerate(bins[1:-1]):
            idx = t+1
            pcb = np.sum(his[:idx])
            pcf = np.sum(his[idx:])
            Wb = pcb * mean_weigth
            Wf = pcf * mean_weigth

            mub = np.sum(cc_arr[:idx] * his[:idx]) / float(pcb)
            muf = np.sum(cc_arr[idx:] * his[idx:]) / float(pcf)
            value = Wb * Wf * (mub - muf) ** 2

            if value > final_value:
                final_thresh = bins[idx]
                final_value = value

        return final_thresh
    