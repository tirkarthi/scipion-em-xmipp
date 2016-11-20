# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

#from constants import *
from convert import writeSetOfVolumes, readSetOfVolumes, getImageLocation, rowFromMd, rowToClass
from pyworkflow.em import ProtClassify3D
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from xmipp3 import getEnviron, XmippMdRow
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
import xmipp


MISSING_WEDGE_X = 0
MISSING_WEDGE_Y = 1
MISSING_PYRAMID = 2
MISSING_CONE = 3


class XmippProtMLTomo(ProtClassify3D):
    """ Align and classify 3D images with missing data regions in Fourier space,
        e.g. subtomograms or RCT reconstructions, by a 3D multi-reference refinement
        based on a maximum-likelihood (ML) target function.

        See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Ml_tomo_v3 for further documentation
    """
    _label = 'mltomo'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='General params')
        form.addParam('inputVols', params.PointerParam, pointerClass="SetOfVolumes",
                      label='Set of volumes',
                      help="Set of input volumes")
        form.addParam('generateRefs', params.BooleanParam, default=True,
                      label='Automatically generate references',
                      help="If set to true, 3D classes will be generated automatically. "
                           "Otherwise you can provide initial reference volumes youself.")
        form.addParam('numberOfReferences', params.IntParam,
                      label='Number of references', default=3, condition="generateRefs",
                      help="Number of references to generate automatically")
        form.addParam('inputRefVols', params.PointerParam, pointerClass="SetOfVolumes",
                      condition="not generateRefs",
                      label='Set of input reference volumes',
                      help="Provide a set of initial reference volumes")
        form.addParam('numberOfIterations', params.IntParam,
                      label='Number of iterations', default=25,
                      help="Maximum number of iterations to perform")
        form.addParam('symmetry', params.StringParam, default='c1', label='Symmetry group',
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry for a description "
                           "of the symmetry groups format. If no symmetry is present, give c1")
        form.addParam('missingDataType', params.EnumParam,
                      choices=['wedge_x', 'wedge_y', 'pyramid', 'cone'], default=0,
                      display=params.EnumParam.DISPLAY_COMBO,
                      label='Missing data regions',
                      help="Provide missing data region type:\n\n"
                           "a) wedge_x for a missing wedge where the tilt axis is along X\n"
                           "b) wedge_y for a missing wedge where the tilt axis is along Y\n"
                           "c) pyramid for a missing pyramid where the tilt axes are along X and Y\n"
                           "d) cone for a missing cone (pointing along Z)")
        form.addParam('missingAng', params.StringParam, default='-60 60',
                      label='Angles of missing data',
                      help='Provide angles for missing data area in the following format:\n\n'
                           'for wedge_x or wedge_y: -60 60\n'
                           'for pyramid: -60 60 -60 60 (for x and y, respectively)\n'
                           'for cone: 45')
        form.addParam('maxCC', params.BooleanParam, default=False,
                      label='Use CC instead of ML',
                      help='Use constrained cross-correlation and weighted averaging instead of ML')

        form.addSection(label='Sampling')
        form.addParam('angSampling', params.FloatParam, default=10.0,
                      label='Angular sampling (deg)',
                      help="Angular sampling rate (in degrees)")
        form.addParam('angSearch', params.FloatParam, default=-1.0,
                      label='Angular search range (deg)',
                      help="Angular search range around orientations of input particles "
                           "(by default, exhaustive searches are performed)")
        form.addParam('globalPsi', params.BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Exhaustive psi search',
                      help="Exhaustive psi searches (only for c1 symmetry)")
        form.addParam('limitTrans', params.FloatParam, default=-1.0, expertLevel=LEVEL_ADVANCED,
                      label='Max shift (px)',
                      help="Maximum allowed shifts (negative value means no restriction)")

        line = form.addLine('Tilt angle limits (deg)', expertLevel=LEVEL_ADVANCED,
                            help='Limits for tilt angle search (in degrees)')
        line.addParam('tiltMin', params.FloatParam, default=0.0, label='Min')
        line.addParam('tiltMax', params.FloatParam, default=180.0, label='Max')
        form.addParam('psiSampling', params.FloatParam, default=-1.0, expertLevel=LEVEL_ADVANCED,
                      label='Psi angle sampling (deg)',
                      help="Angular sampling rate for the in-plane rotations (in degrees)")

        form.addSection(label='Restrictions')
        form.addParam('dim', params.IntParam, default=-1,
                      label='Downscale input to (px)',
                      help="Use downscaled (in fourier space) images of this size (in pixels)")
        form.addParam('maxRes', params.FloatParam, default=0.5,
                      label='Maximum resolution (px^-1)',
                      help="Maximum resolution (in pixel^-1) to use")
        form.addParam('doPerturb', params.BooleanParam, default=False,
                      label='Perturb',
                      help="Apply random perturbations to angular sampling in each iteration")
        form.addParam('dontRotate', params.BooleanParam, default=False,
                      label='Do not rotate',
                      help="Keep orientations fixed, only translate and classify")
        form.addParam('dontAlign', params.BooleanParam, default=False,
                      label='Do not align',
                      help="Keep angles and shifts fixed (otherwise start from random)")
        form.addParam('onlyAvg', params.BooleanParam, default=False,
                      label='Only average',
                      help="Keep orientations and classes, only output weighted averages")
        form.addParam('dontImpute', params.BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Do not impute',
                      help='Use weighted averaging, rather than imputation')
        form.addParam('noImpThresh', params.FloatParam, default=1.0, expertLevel=LEVEL_ADVANCED,
                      condition='dontImpute',
                      label='Threshold for averaging',
                      help='Threshold to avoid division by zero for weighted averaging')
        form.addParam('fixSigmaNoise', params.BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Fix sigma noise',
                      help='Do not re-estimate the standard deviation in the pixel noise')
        form.addParam('fixSigmaOffset', params.BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Fix sigma offsets',
                      help='Do not re-estimate the standard deviation in the origin offsets')
        form.addParam('fixFrac', params.BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Fix model fractions',
                      help='Do not re-estimate the model fractions')

        form.addSection(label='Advanced')
        line = form.addLine('Regularization parameters',
                            help='Regularization parameters (in N/K^2)')

        line.addParam('regIni', params.FloatParam, default=0.0,
                  label='Initial')
        line.addParam('regFinal', params.FloatParam, default=0.0,
                  label='Final')

        form.addParam('numberOfImpIterations', params.IntParam, expertLevel=LEVEL_ADVANCED,
                      label='Iterations in inner imputation loop', default=1,
                      help="Number of iterations for inner imputation loop")
        form.addParam('iterStart', params.IntParam, expertLevel=LEVEL_ADVANCED,
                      label='Initial iteration', default=1,
                      help="Number of initial iteration")
        form.addParam('eps', params.FloatParam, default=5e-5, expertLevel=LEVEL_ADVANCED,
                      label='Stopping criterium',
                      help="Stopping criterium")
        form.addParam('stdNoise', params.FloatParam, default=1.0, expertLevel=LEVEL_ADVANCED,
                      label='Expected noise std',
                      help="Expected standard deviation for pixel noise")
        form.addParam('stdOrig', params.FloatParam, default=3.0, expertLevel=LEVEL_ADVANCED,
                      label='Expected origin offset std (px)',
                      help="Expected standard deviation for origin offset (in pixels)")
        #FIXME next param should be metadata file
        #form.addParam('fracFile', params.PointerParam, pointerClass='FracFile', expertLevel=LEVEL_ADVANCED,
        #              label='Model fractions file',
        #              help='Metadata with expected model fractions (default: even distribution)')
        form.addParam('maskFile', params.PointerParam, pointerClass='Mask', expertLevel=LEVEL_ADVANCED,
                      label='Mask', condition='dontAlign',
                      help='Mask particles; only valid in combination with --dont_align')

        form.addParallelSection(threads=1, mpi=3)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputs')
        self._insertFunctionStep('runMLTomo')
        self._insertFunctionStep('createOutput')
    
    #--------------------------- STEPS functions --------------------------------------------

    def convertInputs(self):
        inputVols = self.inputVols.get()
        # writeSetOfVolumes(self.inputVols.get(), fnVols)
        #TODO: convert input vols, refvols to xmipp format

        self.createMetadata(inputVols)

    def createMetadata(self, inputVols):
        fnVols = self._getExtraPath('input_volumes.xmd')
        missType = ['wedge_x', 'wedge_y', 'pyramid', 'cone']
        missNum = self.missingDataType.get()
        missAng = self.missingAng.get()


        md = xmipp.MetaData()
        for vol in inputVols:
            imgId = vol.getObjId()
            objId = md.addObject()
            row = XmippMdRow()
            row.setValue(xmipp.MDL_ITEM_ID, long(imgId))
            row.setValue(xmipp.MDL_IMAGE, getImageLocation(vol))
            row.setValue(xmipp.MDL_ENABLED, 1)
            row.setValue(xmipp.MDL_MISSINGREGION_NR, missNum + 1)
            row.setValue(xmipp.MDL_ANGLE_ROT, float(0))
            row.setValue(xmipp.MDL_ANGLE_TILT, float(0))
            row.setValue(xmipp.MDL_ANGLE_PSI, float(0))
            row.setValue(xmipp.MDL_SHIFT_X, float(0))
            row.setValue(xmipp.MDL_SHIFT_Y, float(0))
            row.setValue(xmipp.MDL_SHIFT_Z, float(0))
            row.setValue(xmipp.MDL_REF, 1)
            row.setValue(xmipp.MDL_LL, float(1))
            row.writeToMd(md, objId)

        md.write(fnVols, xmipp.MD_APPEND)

        missDataFn = self._getExtraPath('wedges.xmd')
        md = xmipp.MetaData()

        for i in range(1, inputVols.getSize() + 1):
            objId = md.addObject()
            row = XmippMdRow()
            row.setValue(xmipp.MDL_MISSINGREGION_NR, missNum + 1)
            row.setValue(xmipp.MDL_MISSINGREGION_TYPE, missType[missNum])
            row.setValue(xmipp.MDL_MISSINGREGION_THX0, float(0))
            row.setValue(xmipp.MDL_MISSINGREGION_THXF, float(0))
            row.setValue(xmipp.MDL_MISSINGREGION_THY0, float(45))
            row.setValue(xmipp.MDL_MISSINGREGION_THYF, float(0))
            row.writeToMd(md, objId)

        md.write(missDataFn, xmipp.MD_APPEND)

    def runMLTomo(self):
        fnVols = self._getExtraPath('input_volumes.xmd')
        missDataFn = self._getExtraPath('wedges.xmd')
        outDir = self._getExtraPath("results")
        pwutils.makePath(outDir)

        params = ' -i %s' % fnVols
        params += ' --oroot %s' % (outDir + '/mltomo')
        params += ' --iter %d' % self.numberOfIterations.get()
        params += ' --sym %s' % self.symmetry.get()
        params += ' --missing %s' % missDataFn
        params += ' --maxres %0.2f' % self.maxRes.get()
        params += ' --dim %d' % self.dim.get()
        params += ' --ang %0.1f' % self.angSampling.get()
        params += ' --ang_search %0.1f' % self.angSearch.get()
        params += ' --limit_trans %0.1f' % self.limitTrans.get()
        params += ' --tilt0 %0.1f --tiltF %0.1f' % (self.tiltMin.get(), self.tiltMax.get())
        params += ' --psi_sampling %0.1f' % self.psiSampling.get()
        params += ' --reg0 %0.1f --regF %0.1f' % (self.regIni.get(), self.regFinal.get())
        params += ' --impute_iter %d' % self.numberOfImpIterations.get()
        params += ' --istart %d' % self.iterStart.get()
        params += ' --eps %0.2f' % self.eps.get()
        params += ' --pixel_size %0.2f' % self.inputVols.get().getSamplingRate()
        params += ' --noise %0.1f --offset %0.1f' % (self.stdNoise.get(), self.stdOrig.get())
        params += ' --thr %d' % self.numberOfThreads.get()
        if self.generateRefs:
            params += ' --nref %d' % self.numberOfReferences.get()
        else:
            params += ' --ref %s' % self.inputRefVols.get().getFileName()  # FIXME get filename
        if self.doPerturb:
            params += ' --perturb'
        if self.dontRotate:
            params += ' --dont_rotate'
        if self.dontAlign:
            params += ' --dont_align'
        if self.onlyAvg:
            params += ' --only_average'
        if self.dontImpute:
            params += ' --dont_impute'
            params += ' --noimp_threshold %0.1f' % self.noImpThresh.get()
        if self.fixSigmaNoise:
            params += ' --fix_sigma_noise'
        if self.fixSigmaOffset:
            params += ' --fix_sigma_offset'
        if self.fixFrac:
            params += ' --fix_fractions'
        if self.globalPsi:
            params += ' --dont_limit_psirange'
        if self.maxCC:
            params += ' --maxCC'
        if self.dontAlign and self.maskFile:
            params += ' --mask %s' % self.maskFile.get().getFileName() #FIXME get filename

        self.runJob('xmipp_ml_tomo', '%s' % params,
                    env=self.getMLTomoEnviron(), numberOfMpi=self.numberOfMpi.get(),
                    numberOfThreads=self.numberOfThreads.get())
    
    def createOutput(self):
        # output files:
        #   mltomo_ref.xmd contains weight and signal change for output 3D classes
        #   mltomo_refXXXXXX.vol output 3D classes
        #   mltomo_img.xmd contains metadata for input vols
        #   mltomo.fsc
        #   mltomo_refXXXXXX_img.xmd contains metadata for each 3D class

        outputGlobalMdFn = self._getExtraPath('results/mltomo_img.xmd')
        setOfClasses = self._createSetOfClassesVol()
        setOfClasses.setImages(self.inputVols.get())
        self._readSetOfClassesVol(self.inputVols.get(), readSetOfVolumes, setOfClasses, outputGlobalMdFn)
        self._defineOutputs(outputClasses=setOfClasses)
        self._defineSourceRelation(self.inputVols, self.outputClasses)

    
    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []
        return messages

    def _validate(self):
        errors=[]
        return errors

    def _citations(self):
        return ['Scheres2009c']

    #--------------------------- UTILS functions --------------------------------------------
    def getMLTomoEnviron(self):
        env = getEnviron()
        return env

    def _readSetOfClassesVol(self, classBaseSet, readSetOfVolumes, classesSet, outputGlobalMdFn):
        classesMd = xmipp.MetaData(outputGlobalMdFn)

        for objId in classesMd:
            classItem = classesSet.ITEM_TYPE()
            classRow = rowFromMd(classesMd, objId)
            # class id should be set in rowToClass function using MDL_REF
            classItem = rowToClass(classRow, classItem)
            classBaseSet.copyInfo(classesSet.getImages())
            classesSet.append(classItem)
            ref = classItem.getObjId()
            readSetOfVolumes('noname@%s' % outputGlobalMdFn, classItem)

            # Update with new properties of classItem such as _size
            classesSet.update(classItem)
