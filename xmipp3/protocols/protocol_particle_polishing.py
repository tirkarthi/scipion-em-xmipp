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

import os

from pyworkflow import VERSION_2_0
from pyworkflow.protocol.params import PointerParam, IntParam, Positive
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.protocol import ProtAnalysis3D, ProtExtractMovieParticles
from pyworkflow.em.data import Image
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md
from xmipp3.base import XmippMdRow

import xmippLib
from xmipp3.convert import setXmippAttributes, xmippToLocation, coordinateToRow, \
    readSetOfMovieParticles

        
class XmippProtParticlePolishing(ProtAnalysis3D, ProtExtractMovieParticles):
    """Particle polishing."""
    _label = 'particle polishing'
    _lastUpdateVersion = VERSION_2_0
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        self.movieFolder = self._getTmpPath('movie_%(movieId)06d/')
        self.frameRoot = self.movieFolder + 'frame_%(frame)02d'
        myDict = {
            'frameImages': self.frameRoot + '_images',
            'frameMic': self.frameRoot + '.mrc',
            'frameMdFile': self.frameRoot + '_images.xmd',
            'frameCoords': self.frameRoot + '_coordinates.xmd',
            'frameStk': self.frameRoot + '_images.stk'
        }

        self._updateFilenamesDict(myDict)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMovies', PointerParam, label="Input movies", important=True,
                      pointerClass='SetOfMovies', help='Set of aligned movies.')
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles', help='Set of 3D aligned particles.')
        form.addParam('inputVolume', PointerParam, label="Input volume", important=True,
                      pointerClass='Volume',
                      help='Volume to generate reprojections and compare with frame particles.')
        form.addParam('inputCoords', PointerParam, label="Input coordinates", important=True,
                      pointerClass='SetOfCoordinates',
                      help='Set of coordinates to read the input particles from the input movies.')
        form.addParam('inputCTF', PointerParam, pointerClass='SetOfCTF',
                      label="Input CTF", important=True)
        form.addParam('boxSize', IntParam, default=0,
                      label='Particle box size (px)', validators=[Positive],
                      help='In pixels. The box size is the size of the boxed '
                           'particles, actual particles may be smaller than '
                           'this.')
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('input_imgs.xmd')
        self.volFn = self._getTmpPath("volume.vol")

        self._insertFunctionStep("convertStep")
        self._insertFunctionStep("processMovie")
        self._insertFunctionStep("createOutputStep")

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self):
        from xmipp3.convert import writeSetOfParticles
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, self.imgsFn)

        img = ImageHandler()
        img.convert(self.inputVolume.get(), self.volFn)
        xdim=self.inputVolume.get().getDim()[0]

        imgXdim=imgSet.getDim()[0]
        if xdim!=imgXdim:
            self.runJob("xmipp_image_resize","-i %s --dim %d"%(self.volFn,imgXdim),numberOfMpi=1)
    

    def processMovie(self):

        inputMovies = self.inputMovies.get()
        fnCTF = self.inputCTF.get().getFileName()
        for movie in inputMovies:
            movId = movie.getObjId()
            x, y, n = movie.getDim()
            iniFrame, lastFrame, _ = movie.getFramesRange()
            boxSize = self.boxSize.get()

            if movie.hasAlignment():
                shiftX, shiftY = movie.getAlignment().getShifts()
            else:
                shiftX = [0] * (lastFrame - iniFrame + 1)
                shiftY = shiftX

            stkIndex = 0
            movieStk = self._getMovieName(movie, '.stk')
            movieMdFile = self._getMovieName(movie, '.xmd')
            movieMd = md.MetaData()
            frameMd = md.MetaData()
            frameMdImages = md.MetaData()
            frameRow = md.Row()

            if self._hasCoordinates(movie):
                imgh = ImageHandler()

                for frame in range(iniFrame, lastFrame + 1):
                    indx = frame - iniFrame
                    frameName = self._getFnRelated('frameMic', movId, frame)
                    frameMdFile = self._getFnRelated('frameMdFile', movId, frame)
                    coordinatesName = self._getFnRelated('frameCoords', movId,
                                                         frame)
                    frameImages = self._getFnRelated('frameImages', movId, frame)
                    frameStk = self._getFnRelated('frameStk', movId, frame)

                    self._writeXmippPosFile(movie, coordinatesName,
                                            shiftX[indx], shiftY[indx])

                    movieName = imgh.fixXmippVolumeFileName(movie)
                    imgh.convert((frame, movieName), frameName)

                    self.info("Extracting particles")
                    args = '-i %(frameName)s --pos %(coordinatesName)s ' \
                           '-o %(frameImages)s --Xdim %(boxSize)d' % locals()
                    args += " --invert"
                    args += " --fillBorders"
                    args += " --ctfparam %s"%fnCTF

                    self.runJob('xmipp_micrograph_scissor', args)
                    cleanPath(frameName)

                    self.info("Combining particles into one stack.")

                    frameMdImages.read(frameMdFile)
                    frameMd.read('particles@%s' % coordinatesName)
                    frameMd.merge(frameMdImages)

                    for objId in frameMd:
                        stkIndex += 1
                        frameRow.readFromMd(frameMd, objId)
                        location = xmippToLocation(frameRow.getValue(md.MDL_IMAGE))
                        newLocation = (stkIndex, movieStk)
                        imgh.convert(location, newLocation)

                        # Fix the name to be accessible from the Project directory
                        # so we know that the movie stack file will be moved
                        # to final particles folder
                        newImageName = '%d@%s' % newLocation
                        frameRow.setValue(md.MDL_IMAGE, newImageName)
                        frameRow.setValue(md.MDL_MICROGRAPH_ID, long(movId))
                        frameRow.setValue(md.MDL_MICROGRAPH, str(movId))
                        frameRow.setValue(md.MDL_FRAME_ID, long(frame))
                        frameRow.setValue(md.MDL_PARTICLE_ID,
                                          frameRow.getValue(md.MDL_ITEM_ID))
                        frameRow.writeToMd(movieMd, movieMd.addObject())
                    movieMd.addItemId()
                    movieMd.write(movieMdFile)
                    cleanPath(frameStk)


    def createOutputStep(self):

        #For the extracted particles
        inputMovies = self.inputMovies.get()
        particleSet = self._createSetOfMovieParticles()
        particleSet.copyInfo(inputMovies)

        mData = md.MetaData()
        mdAll = md.MetaData()

        self._micNameDict = {}

        for movie in inputMovies:
            self._micNameDict[movie.getObjId()] = movie.getMicName()
            movieName = self._getMovieName(movie)
            movieStk = movieName.replace('.mrc', '.stk')
            movieMdFile = movieName.replace('.mrc', '.xmd')

            # Store particle stack and metadata files in final particles folder
            if os.path.exists(movieStk):
                mData.read(movieMdFile)
                mdAll.unionAll(mData)

        particleMd = self._getPath('movie_particles.xmd')
        mdAll.addItemId()
        mdAll.write(particleMd)
        readSetOfMovieParticles(particleMd, particleSet, removeDisabled=False,
                                postprocessImageRow=self._postprocessImageRow)

        self._defineOutputs(outputParticles=particleSet)
        self._defineSourceRelation(self.inputMovies, particleSet)


    # -------------------------- UTILS functions ------------------------------

    def _writeXmippPosFile(self, movie, coordinatesName, shiftX, shiftY):
        """ Create Xmipp coordinate files to be extracted
        from the frames of the movie.
        """
        coordSet = self.getCoords()

        mData = md.MetaData()
        coordRow = XmippMdRow()

        for coord in coordSet.iterCoordinates(movie.getObjId()):
            coord.shiftX(int(round(float(shiftX))))
            coord.shiftY(int(round(float(shiftY))))
            coordinateToRow(coord, coordRow)
            coordRow.writeToMd(mData, mData.addObject())

        self.info("Writing coordinates metadata: %s, with shifts: %s %s"
                  % (coordinatesName, shiftX, shiftY))
        mData.write('particles@' + coordinatesName)

    def _hasCoordinates(self, movie):
        coordSet = self.getCoords()

        len = 0
        for coord in coordSet.iterCoordinates(movie.getObjId()):
            len  += 1
            break
        if len > 0:
            return True
        else:
            return False

    def _getFnRelated(self, keyFile, movId, frameIndex):
        return self._getFileName(keyFile, movieId=movId, frame=frameIndex)

    def _postprocessImageRow(self, img, imgRow):
        img.setFrameId(imgRow.getValue(md.MDL_FRAME_ID))
        img.setParticleId(imgRow.getValue(md.MDL_PARTICLE_ID))
        micName = self._micNameDict[imgRow.getValue(md.MDL_MICROGRAPH_ID)]
        img.getCoordinate().setMicName(micName)


    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Images evaluated: %i" % self.inputSet.get().getSize())
        summary.append("Volume: %s" % self.inputVolume.getNameId())
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputParticles'):
            methods.append(
                "We projected the volume %s following the directions in %i input images %s." \
                % (
                self.getObjectTag('inputVolume'), self.inputSet.get().getSize(),
                self.getObjectTag('inputSet')))
        return methods

