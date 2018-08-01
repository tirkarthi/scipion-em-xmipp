# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *              David Maluenda Niubo (dmaluenda@cnb.csic.es) [2]
# *
# * [1] SciLifeLab, Stockholm University
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import pyworkflow.em
import pyworkflow.utils as pwutils

from .base import *
from .constants import XMIPP_HOME


_logo = "xmipp_logo.png"
_references = ['delaRosaTrevin2013', 'Sorzano2013']


class Plugin(pyworkflow.em.Plugin):
    _homeVar = XMIPP_HOME
    _pathVars = [XMIPP_HOME]
    _supportedVersions = []

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(XMIPP_HOME, 'xmipp-18.05')

    @classmethod
    def getEnviron(cls, xmippFirst=True):
        """ Create the needed environment for Xmipp programs. """
        environ = pwutils.Environ(os.environ)
        pos = pwutils.Environ.BEGIN if xmippFirst else pwutils.Environ.END
        environ.update({
            'PATH': getXmippPath('bin'),
            'LD_LIBRARY_PATH': getXmippPath('lib'),
        }, position=pos)

        if os.environ['CUDA'] != 'False':  # environ variables are strings not booleans
            environ.update({
                'LD_LIBRARY_PATH': os.environ['NVCC_LIBDIR']
            }, position=pos)

        return environ

    #TODO: Standardize to just: runProgram
    @classmethod
    def runXmippProgram(cls, program, args=""):
        """ Internal shortcut function to launch a Xmipp program. """
        pwutils.runJob(None, program, args, env=cls.getEnviron())

    @classmethod
    def getMatlabEnviron(cls, *toolPaths):
        """ Return an Environment prepared for launching Matlab
        scripts using the Xmipp binding.
        """
        env = pwutils.getEnviron()
        env.set('PATH', os.environ['MATLAB_BINDIR'], pwutils.Environ.BEGIN)
        env.set('LD_LIBRARY_PATH', os.environ['MATLAB_LIBDIR'], pwutils.Environ.BEGIN)
        for toolpath in toolPaths:
            env.set('MATLABPATH', toolpath, pwutils.Environ.BEGIN)
        env.set('MATLABPATH', os.path.join(os.environ[XMIPP_HOME], 'libraries', 'bindings', 'matlab'),
                pwutils.Environ.BEGIN)

        return env



    @classmethod
    def registerPluginBinaries(cls, env):
        scons = env.addPipModule('scons', '2.3.6', target='scons-2.3.6',
                                 default=True, ignoreDefaultDeps=True)

        installCmd = "src/xmipp/xmipp config ; src/xmipp/xmipp check_config;" \
                     " src/xmipp/xmipp compile %d ; src/xmipp/xmipp install"\
                     % env.getProcessors()

        print ("cls.getHome():" + cls.getHome())

        target = "%s/build/bin/xmipp_reconstruct_significant" % cls.getHome()

        # sconsArgs = " ".join([a for a in sys.argv[2:] if not a in env.getTargetNames()])
        xmipp = env.addPackage('xmipp', version='18.05',
                               tar='xmipp-18.05.tgz',
                               commands=[(installCmd, target)],
                               default=True,
                               deps=[scons])  # Old dependencies. Now we check it inside xmipp script:
                                              #   fftw3, scikit, nma, tiff, sqlite, opencv, sh_alignment, hdf5

        xmippBin = env.addPackage('xmippBin', version='18.05',
                                  tar='xmippBinaries-18.05.tgz')

pyworkflow.em.Domain.registerPlugin(__name__)