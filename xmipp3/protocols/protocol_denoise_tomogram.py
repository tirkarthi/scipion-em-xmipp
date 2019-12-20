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

from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import IntParam, EnumParam, LEVEL_ADVANCED, FloatParam, BooleanParam, PointerParam

from pyworkflow.utils.properties import Message
from tomo.objects import Tomogram, SetOfTomograms

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

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.bins_dir='/home/daniel/tomobins/'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        # First we customize the inputParticles param to fit our needs in this protocol
        form.addSection(label='Input')
        form.addParam('inputSetTomograms', PointerParam, pointerClass='SetOfTomograms',
                      label='Set Of Tomograms',
                      help='Select one set of tomograms')
        form.addParam('method', EnumParam,
                      choices=['Anistropic Non-linear Diffusion','BFlow', 'Edge Enhancing Diffusion'], default=0,
                      label='Denoising method',
                      help='Denoising method to use.')
        form.addSection(label='Parameters')
        form.addParam('SigmaGaussian', FloatParam, default=0.5,
                      label='Sigma Gaussian Filter',
                      help='Sigma for initial gaussian filtering.')
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
        inputTomos = self.inputSetTomograms.get()
        self.outputFiles = []
        pre = []
        iter = 1
        for tomo in inputTomos.iterItems():
            stepId= self._insertFunctionStep('denoiseTomogramStep', tomo.getFileName(), iter)
            pre.append(stepId)
            iter += 1

        self._insertFunctionStep('createOutputStep', prerequisites=pre)

    # --------------------------- STEPS functions --------------------------------------------
    def denoiseTomogramStep(self, inp_tomo_path, iter):
        # We start preparing writing those elements we're using as input to keep them untouched
        if self.method.get() == 0:
            print('Denoising by Anistropic Non-linear Diffusion')
            out_tomo_path = self.call_AND(inp_tomo_path, iter)

        elif  self.method.get() == 1:
            print('Denoising by BFlow')
            # call BFlow
            out_tomo_path = self.call_BFlow(inp_tomo_path)

        elif self.method.get() == 2:
            print('Denoising by Edge Enhancing Diffusion')
            out_tomo_path = self.call_EED(inp_tomo_path)

        self.outputFiles.append(out_tomo_path)


    def createOutputStep(self):
        inputTomos = self.inputSetTomograms.get()
        outputTomos = self._createSetOfTomograms()
        outputTomos.copyInfo(inputTomos)
        for i, inp_tomo in enumerate(inputTomos):
            tomo_path = self.outputFiles[i]
            tomo = Tomogram()
            tomo.setLocation(tomo_path)
            outputTomos.append(tomo)

        self._store()
        self._defineOutputs(outputTomograms=outputTomos)
        self._defineSourceRelation(self.inputSetTomograms, outputTomos)


    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        pass

    def _methods(self):
        pass

    # --------------------------- UTILS functions --------------------------------------------
    def create_AND_file(self, input_tomo_path, iter):
        '''Create the csh file that stores the parameters'''
        input_csh = self.bins_dir + 'tomoand.csh'
        output_csh = self._getExtraPath('tomoand_{}.csh'.format(iter))
        output_tomo_path = self._getExtraPath('denoisedEED_'+input_tomo_path.split('/')[-1])
        and_bin = self.bins_dir+'tomoand'

        call = 'sed'
        args = "'{}' {} > {}".format('s/Input_Tomogram.mrc/{}/g'.format(input_tomo_path.replace('/','\/')),
                                     input_csh, output_csh)
        self.runJob(call, args)
        args = "-i '{}' {}".format('s/Output_Tomogram.mrc/{}/g'.format(output_tomo_path.replace('/','\/')),
                                   output_csh)
        self.runJob(call, args)
        args = "-i '{}' {}".format('s/.\/tomoand/{}/g'.format(and_bin.replace('/','\/')), output_csh)
        self.runJob(call, args)
        args=''
        call='chmod 777 {}'.format(output_csh)
        self.runJob(call, args)
        return output_csh, output_tomo_path

    def call_AND(self, inp_tomo_path, iter):
        '''Denoises de tomogram using the AND method'''
        and_csh, out_tomo_path = self.create_AND_file(inp_tomo_path, iter)
        args=''
        self.runJob(and_csh,args)
        return out_tomo_path

    def call_BFlow(self, inp_tomo_path):
        '''Denoises de tomogram using the AND method'''
        params = '-g {} -i {} -s {}'.format(self.SigmaGaussian.get(), self.nIter.get(), self.TimeStep.get())
        out_tomo_path = self._getExtraPath('denoisedBflow_'+inp_tomo_path.split('/')[-1])
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob(self.bins_dir+'tomobflow', args)
        return out_tomo_path

    def call_EED(self, inp_tomo_path):
        '''Denoises de tomogram using the AND method'''
        params = '-g {} -i {} -s {}'.format(self.SigmaGaussian.get(), self.nIter.get(), self.TimeStep.get())
        out_tomo_path = self._getExtraPath('denoisedEED_'+inp_tomo_path.split('/')[-1])
        args = '{} {} {}'.format(params, inp_tomo_path, out_tomo_path)
        self.runJob(self.bins_dir+'tomoeed', args)
        return out_tomo_path


