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
from pyworkflow.protocol.params import (PointerParam, IntParam, Positive,
                                        BooleanParam)
from pwem.protocols import ProtProcessMovies
from pwem.objects import (Movie, SetOfMovies, Set)
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
from pwem import emlib
import numpy as np
from scipy.stats import chisquare, poisson
import time
import pwem.emlib.metadata as md
from math import ceil
import os


class XmippProtMoviePoisson(ProtProcessMovies):
    """Particle polishing."""
    _label = 'movie poisson'
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtProcessMovies.__init__(self, **args)
        self.stepsExecutionMode = cons.STEPS_SERIAL

    # --------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMovies', PointerParam, label="Input movies",
                      important=True,
                      pointerClass='SetOfMovies', help='Set of aligned movies.')
        form.addParam('autoRej', BooleanParam, label='Automatic rejection',
                      default=False, help='Automatic rejection of movies with '
                                          'experimental doses far from the estimated one '
                                          'and distributions not following a Poisson. '
                                          'The rejected movies will appear with the '
                                          'disabled flag in the output set.')

    # --------------------------- STEPS functions ---------------------------------------------------
    def _insertNewMoviesSteps(self, insertedDict, inputMovies):
        """ Insert steps to process new movies (from streaming)
        Params:
            insertedDict: contains already processed movies
            inputMovies: input movies set to be check
        """
        self.numBins = 20
        self.movieMd = md.MetaData()
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

        bin = 2
        start = time.time()
        fnMovie = movie.getFileName()
        movieId = movie.getObjId()
        mov = emlib.Image()
        mov.read(fnMovie)
        movnp = mov.getData()
        frames, _, x, y = movnp.shape
        from os.path import exists
        if not exists(self._getExtraPath("frames.txt")):
            f = open(self._getExtraPath("frames.txt"), "w")
            f.write(str(frames))
            f.close()
        histTot = np.zeros((self.numBins, frames), dtype=int)
        Inp = np.zeros((int(ceil(float(x) / bin)), int(ceil(float(y) / bin))),
                       dtype=int)
        lambdaEst = movie.getAcquisition().getDosePerFrame() * (
                (movie.getSamplingRate()) ** 2)
        lambdaNp = np.zeros((frames), dtype=float)
        countNoPoiss = 0
        end = time.time()
        print('Time loading: %f' % (end - start))
        startA = time.time()
        for f in range(frames):
            Inp[:, :] = movnp[f, :, 0:x:bin, 0:y:bin]
            hist, bins = np.histogram(Inp, bins=range(0, self.numBins))
            histTot[0:len(hist), f] = hist
            lambdaExp = float(sum(hist * bins[0:-1])) / float(sum(hist))
            h, p = chisquare(hist / float(sum(hist)),
                             f_exp=poisson.pmf(range(self.numBins - 1),
                                               lambdaExp))
            lambdaNp[f] = lambdaExp
            if lambdaExp < lambdaEst - (
                    lambdaEst * 0.25) or lambdaExp > lambdaEst + (
                    lambdaEst * 0.10):
                print(
                        "Abnormal dose in frame %i in movie %i - Estimated lambda %f, experimental lambda %f" % (
                f, movieId, lambdaEst, lambdaExp))
            if p < 0.05:
                countNoPoiss += 1
                print(
                        "The experimental data does not follow a Poisson distribution in frame %i in movie %i - h %f, p %f" % (
                f, movieId, h, p))
        np.savetxt(self._getExtraPath('hist_%d.csv' % (int(movieId))), histTot,
                   delimiter=' ')
        endA = time.time()
        print('Time calculating ALL histograms: %f' % (endA - startA))
        lambdaAvg = np.mean(lambdaNp)
        lambdaStd = np.std(lambdaNp)

        mdRow = md.Row()
        mdRow.setValue(md.MDL_ITEM_ID, long(movieId))
        mdRow.setValue(md.MDL_ENABLED, 1)
        mdRow.setValue(md.MDL_IMAGE, fnMovie)
        mdRow.setValue(md.MDL_SAMPLINGRATE, movie.getSamplingRate())
        mdRow.setValue(md.MDL_POISSON_LAMBDA_MEAN, lambdaAvg)
        mdRow.setValue(md.MDL_POISSON_LAMBDA_STD, lambdaStd)
        mdRow.setValue(md.MDL_POISSON_REJECTED_COUNT, countNoPoiss)
        if self.autoRej:
            if countNoPoiss > frames / 2 or lambdaAvg < lambdaEst - (
                    lambdaEst * 0.70) or lambdaAvg > lambdaEst + (
                    lambdaEst * 0.10) or lambdaStd > lambdaAvg * 0.1:
                mdRow.setValue(md.MDL_ENABLED, 0)
                movie.setEnabled(0)
        mdRow.writeToMd(self.movieMd, self.movieMd.addObject())
        movieFn = self._getMovieFn(movie)
        with open(movieFn, 'a') as f:
            f.write('%d %f %f %d \n' % (
            movie.isEnabled(), lambdaAvg, lambdaStd, countNoPoiss))
            f.close()

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
        frames = int(frames)
        histTot = np.zeros((self.numBins, frames * len(listHist)), dtype=int)
        for i, fnHist in enumerate(listHist):
            histAux = np.loadtxt(fnHist)
            _, frames = histAux.shape
            for f in range(frames):
                histTot[0:len(histAux), i * frames + f] = histAux[:, f]
            remove(fnHist)
        np.savetxt(self._getExtraPath('histMatrix.csv'), histTot, delimiter=' ')
        self.movieMd.write(self._getExtraPath('input_movies.xmd'))

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return
        # Load previously done items (from text file)
        doneList = self._readDoneList()
        # Check for newly done items
        newDone = [m.clone() for m in self.listOfMovies
                   if int(m.getObjId()) not in doneList and
                   self._isMovieDone(m)]

        firstTime = len(doneList) == 0
        allDone = len(doneList) + len(newDone)
        # We have finished when there is not more input movies
        # (stream closed) and the number of processed movies is
        # equal to the number of inputs
        self.finished = self.streamClosed and \
                        allDone == len(self.listOfMovies)
        streamMode = Set.STREAM_CLOSED if self.finished else Set.STREAM_OPEN

        if newDone:
            self._writeDoneList(newDone)
        elif not self.finished:
            # If we are not finished and no new output have been produced
            # it does not make sense to proceed and updated the outputs
            # so we exit from the function here
            return

        movieSet = self._loadOutputSet(SetOfMovies, 'movies.sqlite')

        for movie in newDone:
            newMovie = movie.clone()
            movieFn = self._getMovieFn(movie)
            if os.path.exists(movieFn):
                with open(movieFn, 'r') as f:
                    dataMv = f.read()
                    dataMv = dataMv.split()
                    enabled = int(dataMv[0])
                    f.close()
            newMovie.setEnabled(enabled)
            pwutils.cleanPath(movieFn)
            movieSet.append(newMovie)
        self._updateOutputSet('outputMovies', movieSet, streamMode)

        if firstTime:
            self._defineSourceRelation(self.inputMovies, movieSet)

        if self.finished:  # Unlock createOutputStep if finished all jobs
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    # -------------------------- UTILS functions ------------------------------

    def _getFnRelated(self, keyFile, movId, frameIndex):
        return self._getFileName(keyFile, movieId=movId, frame=frameIndex)

    def _loadOutputSet(self, SetClass, baseName):
        """
        Load the output set if it exists or create a new one.
        fixSampling: correct the output sampling rate if binning was used,
        except for the case when the original movies are kept and shifts
        refers to that one.
        """
        setFile = self._getPath(baseName)

        if os.path.exists(setFile) and os.path.getsize(setFile) > 0:
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputMovies = self.inputMovies.get()
        outputSet.copyInfo(inputMovies)

        return outputSet

    def _getMovieFn(self, movie):
        return self._getExtraPath('DATA_movie_%06d.TXT' % movie.getObjId())

    # --------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods



