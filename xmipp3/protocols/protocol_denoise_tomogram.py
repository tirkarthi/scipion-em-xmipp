# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (daniel.delhoyo.gomez@alumnos.upm.es)
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

import pyworkflow.em.metadata as md
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.object import String
from pyworkflow.protocol.params import IntParam, EnumParam, LEVEL_ADVANCED, FloatParam, BooleanParam, PointerParam


from pyworkflow.utils.properties import Message

from xmipp3.convert import writeSetOfParticles, writeSetOfClasses2D, xmippToLocation


class XmippProtDenoiseTomogram(EMProtocol):
    """ Remove particles noise by filtering them.
    This filtering process is based on a projection over a basis created
    from some averages (extracted from classes). This filtering is not
    intended for processing particles. The huge filtering they will be
    passed through is known to remove part of the signal with the noise.
    However this is a good method for clearly see which particle are we
    going to process before it's done.
    """
    _label = 'denoise tomogram'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputTomogram', PointerParam, pointerClass='Tomogram',
                      label=Message.LABEL_INPUT_VOLS,
                      help='Select one tomogram')
        form.addParam('method', EnumParam,
                      choices=['Anistropic Non-linear Diffusion','BFlow', 'Edge Enhancing Diffusion'], default=0,
                      label='Denoising method',
                      help='Denoising method to use.')

        form.addParam('SigmaGaussian', FloatParam, default=0.5,
                      label='Sigma Gaussian Filter',
                      help='Sigma for initial gaussian filtering.')
        #TODO: nIter just different for the default, better way to write it? wizard
        form.addParam('nIter', IntParam, default=40,
                      label='Number of Iterations',
                      help='Number of Iterations of denoising.')
        #
        form.addParam('TimeStep', FloatParam, default=0.1,
                      label='Time Step',
                      help='Time Step for Iterations (max 0.15)')

        #form.addSection(label='Method parameters')
        #form.addParam('Lambda', IntParam, default=-1, condition='method==2' or 'method==2',
         #             label='Lambda EED',# expertLevel=LEVEL_ADVANCED,
          #            help='Apply a density scaling.')


    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Insert every step of the protocol"""

        # Convert input images if necessary
        self._insertFunctionStep('denoiseTomogram')

        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions --------------------------------------------
    def denoiseTomogram(self):
        # We start preparing writing those elements we're using as input to keep them untouched
        args=''
        outfile='outputtomogram.mrc'
        call='tomoeed {} inputtomogram.mrc {}'.format(args,outfile)

        self.runJob(call)

        self.outputMd = String(outfile)



        imagesMd = self._getPath('images.xmd')
        writeSetOfParticles(self.inputParticles.get(), imagesMd)
        classesMd = self._getPath('classes.xmd')
        writeSetOfClasses2D(self.inputClasses.get(), classesMd)

        fnRoot = self._getExtraPath('pca')
        fnRootDenoised = self._getExtraPath('imagesDenoised')

        args = '-i Particles@%s --oroot %s --eigenvectors %d --maxImages %d' % (
        imagesMd, fnRoot, self.maxPCABases.get(), self.maxClasses.get())
        self.runJob("xmipp_image_rotational_pca", args)

        N = min(self.maxPCABases.get(), self.PCABases2Project.get())
        args = '-i %s -o %s.stk --save_metadata_stack %s.xmd --basis %s.stk %d' \
               % (imagesMd, fnRootDenoised, fnRootDenoised, fnRoot, N)

        self.runJob("xmipp_transform_filter", args)

        self.outputMd = String('%s.stk' % fnRootDenoised)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()

        partSet.copyInfo(imgSet)
        partSet.copyItems(imgSet,
                          updateItemCallback=self._updateLocation,
                          itemDataIterator=md.iterRows(self.outputMd.get(), sortByLabel=md.MDL_ITEM_ID))

        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            summary.append('PCA basis created by using %d classes' % len(self.inputClasses.get()))
            summary.append('Max. number of classes defined for PCA basis creation: %d' % self.maxClasses.get())
            summary.append('Max. number of PCA bases defined for PCA basis creation: %d' % self.maxPCABases.get())
            summary.append('PCA basis on which to project for denoising: %d' % self.PCABases2Project.get())
        return summary

    def _validate(self):
        pass

    def _methods(self):
        methods = []
        if not hasattr(self, 'outputParticles'):
            methods.append("Output particles not ready yet.")
        else:
            methods.append('An input dataset of %d particles was filtered creating a PCA basis (%d components) with '
                           'xmipp_image_rotational_pca and projecting the dataset into that base with xmipp_transform_filter.' \
                           % (len(self.inputParticles.get()), len(self.inputClasses.get())))
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _updateLocation(self, item, row):
        index, filename = xmippToLocation(row.getValue(md.MDL_IMAGE))
        item.setLocation(index, filename)

