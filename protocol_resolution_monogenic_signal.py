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



from pyworkflow.protocol.params import (PointerParam, StringParam, BooleanParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol.protocol_3d import ProtRefine3D
from convert import readSetOfVolumes
from pyworkflow.utils import getExt, removeExt
from os.path import basename, abspath

OUTPUT_RESOLUTION_FILE = 'mgresolution.vol'


class XmippProtMonoRes(ProtRefine3D):
    """    
    Given a map the protocol assigns local resolutions to each pixel of the map.
    """
    _label = 'monogenic resolution'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('halfVolums', BooleanParam, default=False,
                      label="Would you like to use half volums?",
                      help='The noise estimation for determining the local resolution '
                           'is performed via half volumes.')

        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input Volume",
                      help='Select a volume for determining its local resolucion.')

        form.addParam('inputVolume2', PointerParam, pointerClass='Volume',
                      label="Second Half Volume", condition='halfVolums',
                      help='Select a second volume for determining a local resolucion.')

        form.addParam('provideMaskInHalves', BooleanParam, default=False,
                      condition='halfVolums',
                      label="Use mask with halves volumes?",
                      help='Sometimes the volume is in an sphere, then this option ought to be selected')

        form.addParam('Mask', PointerParam, pointerClass='VolumeMask', allowsNull=True,
                      condition='(provideMaskInHalves and halfVolums) or (not halfVolums)',
                      label="Binary Mask",
                      help='The mask determines the those points where the macromolecle is')

        form.addParam('symmetry', StringParam, default='c1',
                      label="Symmetry",
                      help='Symmetry group. By default = c1')

        line = form.addLine('Resolution Range (A)',
                            help="If the user knows the range of resolutions or only a"
                                 " range of frequency needs to be analysed", expertLevel=LEVEL_ADVANCED)

        line.addParam('minRes', FloatParam, default=1, label='Min')
        line.addParam('maxRes', FloatParam, default=100, label='Max')

        form.addParam('significance', FloatParam, label="Significance", default=0.95, expertLevel=LEVEL_ADVANCED,
                      help='The resolution is computed performing hypothesis tests. This parameter determines'
                           ' the significance for that test.')
        form.addParam('exact', BooleanParam, default=True, expertLevel=LEVEL_ADVANCED,
                      label="Use gausian resolution",
                      help='The noise estimation can be performed exact (slow) or approximated (fast)'
                           'ussually there has not difference between them')
        form.addParam('filterInput', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label="Filter input volume with local resolution?",
                      help='The input map is locally filtered at the local resolution map.')
        form.addParam('trimming', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label="Remove atypical resolution values",
                      help='In some situations bad voxels appear. This option allow to remove those voxels')

        form.addParam('kValue', FloatParam, label="Trimming Value", condition='trimming',
                      default=5,
                      help='This value performs post-processing, smoothing the output resolutions.'
                           'The resolutions in this percentile, will be changed by the mean value')

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self.micsFn = self._getPath()
        self.vol1Fn = self.inputVolume.get().getFileName()
        self.maskFn = self.Mask.get().getFileName()

        if self.halfVolums.get() is True:
            self.vol2Fn = self.inputVolume2.get().getFileName()

            # Convert input into xmipp Metadata format
        convertId = self._insertFunctionStep('convertInputStep', )

        MS = self._insertFunctionStep('resolutionMonogenicSignalStep', prerequisites=[convertId])

        self._insertFunctionStep('createOutputStep', prerequisites=[MS])

        self._insertFunctionStep("createChimeraScript")

        self._insertFunctionStep("createHistrogram")

    def convertInputStep(self):
        """ Read the input volume.
        """
        extVol1 = getExt(self.vol1Fn)
        if (extVol1 == '.mrc') or (extVol1 == '.map'):
            self.vol1Fn = self.vol1Fn + ':mrc'

        if self.halfVolums.get() is True:
            extVol2 = getExt(self.vol2Fn)
            if (extVol2 == '.mrc') or (extVol2 == '.map'):
                self.vol2Fn = self.vol2Fn + ':mrc'

        extMask = getExt(self.maskFn)
        if (extMask == '.mrc') or (extMask == '.map'):
            self.maskFn = self.maskFn + ':mrc'

    def resolutionMonogenicSignalStep(self):

        if self.halfVolums.get() is False:
            params = ' --vol %s' % self.vol1Fn
            params += ' --mask %s' % self.maskFn
        else:
            params = ' --vol %s' % self.vol1Fn
            params += ' --vol2 %s' % self.vol2Fn
            if self.provideMaskInHalves.get() is True:
                params += ' --mask %s' % self.maskFn
            else:
                params += ' --mask %s' % ''

        params += ' -o %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --sampling_rate %f' % self.inputVolume.get().getSamplingRate()
        params += ' --number_frequencies %f' % 50
        params += ' --minRes %f' % self.minRes.get()
        params += ' --maxRes %f' % self.maxRes.get()
        params += ' --chimera_volume %s' % self._getExtraPath('MG_Chimera_resolution.vol')
        params += ' --linear '
        params += ' --sym %s' % self.symmetry.get()
        params += ' --significance %s' % self.significance.get()
        if self.exact.get() is False:
            params += ' --exact'
        if self.filterInput.get():
            params += ' --filtered_volume %s' % self._getExtraPath('filteredMap.vol')
        else:
            params += ' --filtered_volume %s' % ''

        if self.trimming.get() is True:
            params += ' --trimmed %f' % self.kValue.get()
        else:
            params += ' --trimmed %f' % 0

        self.runJob('xmipp_resolution_monogenic_signal', params)

    def createChimeraScript(self):
        fnRoot = "extra/"
        scriptFile = self._getPath('Chimera_resolution.cmd')

        fnbase = removeExt(self.inputVolume.get().getFileName())
        ext = getExt(self.inputVolume.get().getFileName())
        fninput = abspath(fnbase + ext[0:4])
        fhCmd = open(scriptFile, 'w')
        fhCmd.write("open %s\n" % fninput)
        fhCmd.write("open %s\n" % (fnRoot + "MG_Chimera_resolution.vol"))
        smprt = self.inputVolume.get().getSamplingRate()
        fhCmd.write("volume #1 voxelSize %s\n" % (str(smprt)))
        fhCmd.write("vol #1 hide\n")
        fhCmd.write("scolor #0 volume #1 cmap rainbow reverseColors True\n")
        fhCmd.close()

    def createHistrogram(self):

        params = ' -i %s' % self._getExtraPath(OUTPUT_RESOLUTION_FILE)
        params += ' --mask binary_file %s' % self.maskFn
        params += ' --steps %f' % 30
        params += ' --range %f %f' % (self.minRes.get(), self.maxRes.get())
        params += ' -o %s' % self._getExtraPath('hist.xmd')

        self.runJob('xmipp_image_histogram', params)

    def createOutputStep(self):
        volume_path = self._getExtraPath(OUTPUT_RESOLUTION_FILE)

        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate(self.inputVolume.get().getSamplingRate())

        readSetOfVolumes(volume_path, volumesSet)

        self._defineOutputs(outputVolume=volumesSet)
        self._defineSourceRelation(self.inputVolume, volumesSet)

        if self.filterInput.get():
            volume_filtered_path = self._getExtraPath('filteredMap.vol')
            volumesSet2 = self._createSetOfVolumes()
            volumesSet2.setSamplingRate(self.inputVolume.get().getSamplingRate())
            readSetOfVolumes(volume_filtered_path, volumesSet2)
            self._defineOutputs(outputVolume=volumesSet2)
            self._defineSourceRelation(self.inputVolume, volumesSet2)

    # --------------------------- INFO functions --------------------------------------------
    def _validate(self):

        validateMsgs = []
        if not self.inputVolume.get().hasValue():
            validateMsgs.append('Please provide input volume.')
        return validateMsgs

    def _methods(self):
        messages = []
        if hasattr(self, 'outputVolume'):
            messages.append(
                'The local resolution is performed [Publication: Not yet]')
        return messages

    def _citations(self):
        return ['Not yet']
