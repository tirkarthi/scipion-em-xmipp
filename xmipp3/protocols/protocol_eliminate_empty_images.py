# *****************************************************************************
# *
# * Authors:     Tomas Majtner         tmajtner@cnb.csic.es (2017)
# *              David Maluenda        dmaluenda@cnb.csic.es
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
# *****************************************************************************

import os
from datetime import datetime

import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow.protocol.constants as cons
import pyworkflow.protocol.params as param
from pyworkflow import VERSION_2_0
from pyworkflow.em.protocol import ProtClassify2D
from pyworkflow.object import Set
from pyworkflow.utils import cleanPath
from pyworkflow.utils.properties import Message

from xmipp3.convert import writeSetOfParticles, \
    readSetOfParticles, setXmippAttributes


class XmippProtEliminateEmptyBase(ProtClassify2D):
    """ Base to eliminate images using statistical methods
    (variance of variances of sub-parts of input image) eliminates those samples,
    where there is no object/particle (only noise is presented there).
    Threshold parameter can be used for fine-tuning the algorithm for type of data.
    """
    _lastUpdateVersion = VERSION_2_0

    def __init__(self, **args):
        ProtClassify2D.__init__(self, **args)
        self.stepsExecutionMode = em.STEPS_PARALLEL

    def addAdvancedParams(self, form):
        form.addParam('addFeatures', param.BooleanParam, default=False,
                      label='Add features', expertLevel=param.LEVEL_ADVANCED,
                      help='Add features used for the ranking to each '
                           'one of the input particles.')
        form.addParam('useDenoising', param.BooleanParam, default=True,
                      label='Turning on denoising',
                      expertLevel=param.LEVEL_ADVANCED,
                      help='Option for turning on denoising method '
                           'while computing emptiness feature.')
        form.addParam('denoising', param.IntParam, default=5,
                      expertLevel=param.LEVEL_ADVANCED,
                      condition='useDenoising',
                      label='Denoising factor:',
                      help='Factor to be used during Gaussian blurring. '
                           'Higher value applies stronger denoising, '
                           'could be more precise but also slower.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.outputSize = 0
        self.check = None
        self.stepCount = 0
        self.fnInputMd = self._getExtraPath("input%d.xmd")
        self.fnOutputMd = self._getExtraPath("output.xmd")
        self.fnElimMd = self._getExtraPath("eliminated.xmd")
        self.fnOutMdTmp = self._getExtraPath("outTemp.xmd")
        self.fnElimMdTmp = self._getExtraPath("elimTemp.xmd")
        checkStep = self._insertNewPartsSteps()
        self._insertFunctionStep('createOutputStep',
                                 prerequisites=checkStep, wait=True)

    def _insertNewPartsSteps(self):
        deps = []
        self.stepCount += 1
        stepId = self._insertFunctionStep('eliminationStep',
                                          self.stepCount,
                                          prerequisites=[])
        deps.append(stepId)
        return deps

    def eliminationStep(self, stepId):
        fnInputMd = self.fnInputMd % stepId
        partsSet = self.prepareImages()

        if self.check == None:
            writeSetOfParticles(partsSet, fnInputMd,
                                alignType=em.ALIGN_NONE, orderBy='creation')
        else:
            writeSetOfParticles(partsSet, fnInputMd,
                                alignType=em.ALIGN_NONE, orderBy='creation',
                                where='creation>"' + str(self.check) + '"')

        self.specialBehavoir(partsSet)

        self.lenPartsSet = len(partsSet)

        args = "-i %s -o %s -e %s -t %f" % (fnInputMd, self.fnOutputMd,
                                            self.fnElimMd, self.threshold.get())
        if self.addFeatures:
            args += " --addFeatures"
        if self.useDenoising:
            args += " --useDenoising -d %f" % self.denoising.get()
        self.runJob("xmipp_image_eliminate_empty_particles", args)

    def specialBehavoir(self, inSet):
        pass

    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None

    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all particles
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'

    def createOutputStep(self):
        pass

    def _stepsCheck(self):
        # Input particles set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _checkNewInput(self):
        # Check if there are new particles to process from the input set
        partsFile = self.getInput().getFileName()
        self.lastCheck = getattr(self, 'lastCheck', datetime.now())
        mTime = datetime.fromtimestamp(os.path.getmtime(partsFile))
        # If the input movies.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime:
            return None
        self.lastCheck = datetime.now()

        outputStep = self._getFirstJoinStep()

        self.prepareImages()

        fDeps = self._insertNewPartsSteps()
        if outputStep is not None:
            outputStep.addPrerequisites(*fDeps)
        self.updateSteps()

    def prepareImages(self):
        pass

    def _checkNewOutput(self):
        if getattr(self, 'finished', False):
            return

        self.finished = self.streamClosed and self.outputSize == self.lenPartsSet

        self.createOutputs()

        if self.finished:  # Unlock createOutputStep if finished all jobs
            cleanPath(self._getPath('particlesAUX.sqlite'))
            cleanPath(self._getPath('averagesAUX.sqlite'))
            outputStep = self._getFirstJoinStep()
            if outputStep and outputStep.isWaiting():
                outputStep.setStatus(cons.STATUS_NEW)

    def createOutputs(self):
        """ To be implemented by childs. (create, fill and close the outputSet)
        """
        pass

    def _loadOutputSet(self, SetClass, baseName):
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = SetClass(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = SetClass(filename=setFile)
            outputSet.setStreamState(outputSet.STREAM_OPEN)

        inputs = self.inputImages
        outputSet.copyInfo(inputs)

        return outputSet

    def _updateOutputSet(self, outputName, outputSet, state=Set.STREAM_OPEN):
        outputSet.setStreamState(state)
        if self.hasAttribute(outputName):
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protocol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            self._defineOutputs(**{outputName: outputSet})
            self._defineTransformRelation(self.getInput(), outputSet)
            self._store(outputSet)

        # Close set databaset to avoid locking it
        outputSet.close()

    # --------------------------- UTILS functions -----------------------------
    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_SCORE_BY_EMPTINESS)
        if row.getValue(md.MDL_ENABLED) <= 0:
            item._appendItem = False
        else:
            item._appendItem = True

    def getInput(self):
        pass


class XmippProtEliminateEmptyParticles(XmippProtEliminateEmptyBase):
    """ Takes a set of particles and using statistical methods
    (variance of variances of sub-parts of input image) eliminates those samples,
    where there is no object/particle (only noise is presented there).
    Threshold parameter can be used for fine-tuning the algorithm for type of data.
    """

    _label = 'eliminate empty particles'

    def __init__(self, **args):
        XmippProtEliminateEmptyBase.__init__(self, **args)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        # - - - F O R   P A R T I C L E S - - -
        form.addParam('inputParticles', param.PointerParam, important=True,
                      label="Input particles", pointerClass='SetOfParticles',
                      help='Select the input particles to be classified.')
        form.addParam('threshold', param.FloatParam, default=1.1,
                      label='Threshold used in elimination:',
                      help='Higher threshold => more particles will be '
                           'eliminated. Set to -1 for no elimination, even so '
                           'the "xmipp_scoreEmptiness" value will be attached to '
                           'every paricle for a posterior inspection.')

        self.addAdvancedParams(form)

    # --------------------------- INSERT steps functions ----------------------
    def specialBehavoir(self, partsSet):
        for p in partsSet.iterItems(orderBy='creation', direction='DESC'):
            self.check = p.getObjCreation()
            break
        partsSet.close()

    def createOutputs(self):
        streamMode = Set.STREAM_CLOSED if getattr(self, 'finished', False) \
            else Set.STREAM_OPEN

        def updateOutputs(mdFn, suffix):
            newData = os.path.exists(mdFn)
            lastToClose = getattr(self, 'finished', False) and \
                          hasattr(self, '%sParticles' % suffix)
            if newData or lastToClose:
                outSet = self._loadOutputSet(em.SetOfParticles,
                                             '%sParticles.sqlite' % suffix)
                if newData:
                    partsSet = self._createSetOfParticles("AUX")
                    readSetOfParticles(mdFn, partsSet)
                    outSet.copyItems(partsSet,
                                     updateItemCallback=self._updateParticle,
                                     itemDataIterator=md.iterRows(mdFn,
                                                    sortByLabel=md.MDL_ITEM_ID))
                    self.outputSize = self.outputSize + len(partsSet)
                self._updateOutputSet('%sParticles'%suffix, outSet, streamMode)
                cleanPath(mdFn)

        updateOutputs(self.fnOutMdTmp, 'output')
        updateOutputs(self.fnElimMdTmp, 'eliminated')

    def getInput(self):
        return self.inputParticles.get()

    def prepareImages(self):
        self.inputImages = self.getInput()
        partsSet = self.inputImages
        # partsFile = self.inputImages.getFileName()
        # partsSet = SetOfParticles(filename=partsFile)
        partsSet.loadAllProperties()
        self.streamClosed = partsSet.isStreamClosed()

        return partsSet


DISCARDED = 0
ACCEPTED = 1
ATTACHED = 2

class XmippProtEliminateEmptyClasses(XmippProtEliminateEmptyBase):
    """ Takes a set of classes (or averages) and using statistical methods
    (variances of sub-parts of input image) eliminates those samples,
    where there is no object/particle (only noise is presented there).
    Threshold parameter can be used for fine-tuning the algorithm for
    type of data. Also discards classes with less population than a given
    persentage.
    """

    _label = 'eliminate empty classes'

    def __init__(self, **args):
        XmippProtEliminateEmptyBase.__init__(self, **args)
        self.enableCls = {}

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        # - - - F O R   C L A S S E S - - -
        form.addParam('inputClasses', param.PointerParam, important=True,
                      pointerClass='SetOfClasses, SetOfAverages',
                      label="Input classes",
                      help='Select the input averages to be classified.')
        form.addParam('threshold', param.FloatParam, default=8.0,
                      label='Threshold used in elimination',
                      help='Higher threshold => more particles will be '
                           'eliminated. Set to -1 for no elimination, even so '
                           'the "xmipp_scoreEmptiness" value will be attached to '
                           'every paricle for a posterior inspection.')
        form.addParam('usePopulation', param.BooleanParam, default=True,
                      label='Use class population',
                      help="Use class population to reject a class.")
        form.addParam('minPopulation', param.FloatParam, default=20,
                      label='Min. population (%)',
                      condition="usePopulation",
                      expertLevel=param.LEVEL_ADVANCED,
                      help="Minimum population to accept a class.\n"
                           "Classes with less population than the mean population "
                           "times this value will be rejected.")

        self.addAdvancedParams(form)

    # --------------------------- INSERT steps functions ----------------------
    def _validate(self):
        errors=[]
        if isinstance(self.getInput(), em.SetOfImages) and self.usePopulation.get:
            errors.append("Using population to reject classes is not possible "
                          "with Averages as input.\nPlease, introduce a "
                          "setOfClasses to analyse or disable the _use class "
                          "population_ option.")

    def specialBehavoir(self, partSet):
        idsToCheck = []
        for p in partSet.iterItems(orderBy='creation', direction='ASC'):
            self.check = p.getObjCreation()
            idsToCheck.append(p.getObjId())
        partSet.close()

        self.rejectByPopulation(idsToCheck)

    def createOutputs(self):
        streamMode = Set.STREAM_CLOSED if getattr(self, 'finished', False) \
            else (Set.STREAM_CLOSED if self.streamClosed else Set.STREAM_OPEN)

        def updateOutputs(mdFn, suffix, newData):
            lastToClose = getattr(self, 'finished', False) and \
                          hasattr(self, '%sClasses' % suffix)
            enableOut = {}
            if newData or lastToClose:
                outSet = self._loadOutputSet(em.SetOfAverages,
                                             '%sAverages.sqlite' % suffix)
                if newData:
                    # if new data, we read it
                    partsSet = self._createSetOfParticles("AUX")
                    readSetOfParticles(mdFn, partsSet)
                    # updating the enableCls dictionary
                    print(" - %s Averages:" % ("ACCEPTED" if suffix == 'output' else "DISCARTDED"))
                    for part in partsSet:
                        partId = part.getObjId()
                        if partId not in self.enableCls:
                            # this happends when a classifier give an empty class
                            continue
                        # - accept if we are in accepted and the current is accepted
                        # - discard if we are in the discarted scope and any
                        currentStatus = self.enableCls[partId]
                        decision = suffix == 'output' and currentStatus == ACCEPTED
                        enableOut[partId] = ACCEPTED if decision else DISCARDED
                        print("%d: %s -> %s" % (partId, currentStatus, decision))
                    # updating the Averages set
                    outSet.copyItems(partsSet,
                                     updateItemCallback=self._updateParticle,
                                     itemDataIterator=md.iterRows(mdFn,
                                                    sortByLabel=md.MDL_ITEM_ID))
                    self.outputSize = self.outputSize + len(partsSet)

                self._updateOutputSet('%sAverages' % suffix, outSet, streamMode)
                cleanPath(mdFn)

            return enableOut

        accOut = updateOutputs(self.fnOutMdTmp, 'output',
                               os.path.exists(self.fnOutMdTmp))

        discOut = updateOutputs(self.fnElimMdTmp, 'eliminated',
                                os.path.exists(self.fnElimMdTmp))
        print(type(accOut))
        self.createOutputClasses('output', streamMode, accOut)
        self.createOutputClasses('eliminated', streamMode, discOut)

    # ------------- UTILS Fuctions ------------------------------------
    def prepareImages(self):
        inSet = self.getInput()

        if isinstance(inSet, em.SetOfImages):
            firstRep = inSet.getFirstItem()
            getImage = lambda item: item.clone()
            self.classesDict = None
        else:  # SetOfClasses
            firstRep = inSet.getFirstItem().getFirstItem()
            getImage = lambda item: item.getRepresentative().clone()
            self.classesDict = {cls.getObjId(): cls.getSize() for cls in inSet}

        self.inputImages = self._createSetOfAverages("AUX")
        self.inputImages.enableAppend()
        self.inputImages.copyAttributes(firstRep, '_samplingRate')
        self.inputImages.copyAttributes(inSet, '_streamState')
        self.streamClosed = self.inputImages.isStreamClosed()

        for item in inSet:
            self.inputImages.append(getImage(item))

        self.inputImages.write()
        self._store(self.inputImages)

        return self.inputImages

    def rejectByPopulation(self, ids):
        if self.usePopulation.get() and self.classesDict is not None:
            sizeDict = {clsId: size for clsId, size
                        in self.classesDict.iteritems()
                        if clsId in ids}

            meanPop = sum(sizeDict.values())/len(sizeDict)

            for clsId, size in sizeDict.iteritems():
                decision = int(size*100 > meanPop * self.minPopulation.get())
                self.enableCls[clsId] = ACCEPTED if decision else DISCARDED
        else:
            self.enableCls = {clsId: ACCEPTED for clsId in ids}

    def createOutputClasses(self, suffix, streamingState, enableDict):
        if not self.classesDict:
            # If there are no classes, nothing to do
            return

        baseName = '%sClasses.sqlite' % suffix
        setFile = self._getPath(baseName)
        if os.path.exists(setFile):
            outputSet = em.SetOfClasses2D(filename=setFile)
            outputSet.loadAllProperties()
            outputSet.enableAppend()
        else:
            outputSet = em.SetOfClasses2D(filename=setFile)
            outputSet.setStreamState(streamingState)

        outputSet.copyInfo(self.getInput())  # if fails, delete

        if enableDict:
            # FIXME: Review this !!!
            decision = ACCEPTED if suffix == 'output' else DISCARDED
            print("in createOutputClasses... %s" % ('ACCEPTED' if suffix == 'output' else 'DISCARDED'))
            desiredIds = [ids for ids, enable in enableDict.iteritems()
                          if enable == decision]
            print("self.enableCls: %s" % enableDict)
            print("desiredIds: %s" % desiredIds)

            for cls in self.getInput():
                repId = cls.getObjId()
                if repId in desiredIds:
                    newClass = cls.clone()

                    outputSet.append(newClass)

            for cls in outputSet:
                repId = cls.getObjId()
                cls.copyItems(self.getInput()[repId])

                outputSet.update(cls)

            self.enableCls.update({idItem: ATTACHED for idItem in desiredIds})

        outputSet.setStreamState(streamingState)
        outputName = '%sClasses' % suffix
        if self.hasAttribute(outputName):
            outputSet.write()  # Write to commit changes
            outputAttr = getattr(self, outputName)
            # Copy the properties to the object contained in the protocol
            outputAttr.copy(outputSet, copyId=False)
            # Persist changes
            self._store(outputAttr)
        else:
            self._defineOutputs(**{outputName: outputSet})
            self._defineSourceRelation(self.inputClasses, outputSet)
            self._store(outputSet)

        # Close set databaset to avoid locking it
        outputSet.close()

    def getInput(self):
        return self.inputClasses.get()
