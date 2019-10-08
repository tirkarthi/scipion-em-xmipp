# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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
from pyworkflow.protocol.params import PointerParam, IntParam, Positive
from pyworkflow.em.protocol import ProtProcessMovies
from pyworkflow.em.data import Image, Movie, SetOfMovies, SetOfImages
from pyworkflow.protocol.constants import STEPS_SERIAL
import pyworkflow.protocol.constants as cons

import xmippLib
import numpy as np
from scipy.stats import chisquare, poisson


class XmippProtMoviePoisson(ProtProcessMovies):
    """Particle polishing."""
    _label = 'movie poisson'
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtProcessMovies.__init__(self, **args)
        self.stepsExecutionMode = STEPS_SERIAL

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.movieFolder = self._getTmpPath('movie_%(movieId)06d/')
        self.frameRoot = self.movieFolder + 'frame_%(frame)02d'

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMovies', PointerParam, label="Input movies",
                      important=True,
                      pointerClass='SetOfMovies', help='Set of aligned movies.')


    # --------------------------- STEPS functions ---------------------------------------------------
    def _insertNewMoviesSteps(self, insertedDict, inputMovies):
        """ Insert steps to process new movies (from streaming)
        Params:
            insertedDict: contains already processed movies
            inputMovies: input movies set to be check
        """
        self.numBins = 20
        deps = []
        if isinstance(self.inputMovies.get(), Movie):
            movie = self.inputMovies.get()
            if movie.getObjId() not in insertedDict:
                stepId = self._insertMovieStep(movie)
                deps.append(stepId)
                insertedDict[movie.getObjId()] = stepId
        else:
            # For each movie insert the step to process it
            for idx, movie in enumerate(self.inputMovies.get()):
                if movie.getObjId() not in insertedDict:
                    stepId = self._insertMovieStep(movie)
                    deps.append(stepId)
                    insertedDict[movie.getObjId()] = stepId
        return deps

    def _convertInputStep(self):
        pass

    def _processMovie(self, movie):

        fnMovie = movie.getFileName()
        movieId = movie.getObjId()
        mov = xmippLib.Image()
        mov.read(fnMovie)
        movnp = mov.getData()
        frames, _, x, y = movnp.shape
        from os.path import exists
        if not exists(self._getExtraPath("frames.txt")):
            f = open(self._getExtraPath("frames.txt"), "w")
            f.write(str(frames))
            f.close()
        histTot = np.zeros((self.numBins, frames),dtype=int)
        Inp=np.zeros((x,y),dtype=int)
        lambdaEst = movie.getAcquisition().getDosePerFrame() * ((movie.getSamplingRate()) ** 2)
        # pix = np.zeros((x*y),dtype=float)
        for f in range(frames):
            Inp[:,:] = movnp[f,:,:,:]
            #if f==0 and int(movieId)==1: #REMOVEEEEE
            #    np.savetxt(self._getExtraPath('frame_%d_movie_%d.csv' % (f,movieId)), Inp, delimiter=' ')
            hist, bins = np.histogram(Inp, bins=range(0, self.numBins))
            histTot[0:len(hist),f] = hist
            lambdaExp = float(sum(hist*bins[0:-1]))/float(sum(hist))
            h, p = chisquare(hist / float(sum(hist)), f_exp=poisson.pmf(range(self.numBins - 1), lambdaExp))
            if lambdaExp<lambdaEst-(lambdaEst*0.25) or lambdaExp>lambdaEst+(lambdaEst*0.25):
                print("Anormal dose: check frame %i in movie %i " %(f,movieId))
                print("Estimated lambda %f, experimental lambda %f" % (lambdaEst, lambdaExp))
            if p < 0.05:
                print( "The experimental data does not follow a Poisson distribution. Frame %i in movie %i " % (f, movieId))
                print("h %f, p %f" % (h, p))
        np.savetxt(self._getExtraPath('hist_%d.csv'%(int(movieId))), histTot, delimiter=' ')


    def _checkNewInput(self):
        if isinstance(self.inputMovies.get(), SetOfMovies):
            ProtProcessMovies._checkNewInput(self)

    def createOutputStep(self):
        from os import remove
        import glob
        listHist = glob.glob(self._getExtraPath('hist*csv'))
        f = open(self._getExtraPath("frames.txt"), "r")
        frames = f.readline()
        f.close()
        frames=int(frames)
        histTot = np.zeros((self.numBins, frames * len(listHist)),dtype=int)
        for i,fnHist in enumerate(listHist):
            histAux = np.loadtxt(fnHist)
            _, frames = histAux.shape
            for f in range(frames):
                histTot[0:len(histAux), i*frames+f] = histAux[:,f]
            remove(fnHist)
        np.savetxt(self._getExtraPath('histMatrix.csv'), histTot, delimiter=' ')

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m.clone() for m in self.listOfMovies
                   if int(m.getObjId()) not in doneList and
                   self._isMovieDone(m)]

        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input movies
        # (stream closed) and the number of processed movies is
        # equal to the number of inputs
        self.finished = self.streamClosed and \
                        allDone == len(self.listOfMovies)

        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)


    # -------------------------- UTILS functions ------------------------------
    def _getFnRelated(self, keyFile, movId, frameIndex):
        return self._getFileName(keyFile, movieId=movId, frame=frameIndex)



    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods
