# **************************************************************************
# *
# * Authors:     David Maluenda (dmaluenda@cnb.csic.es)
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

import pickle

from pyworkflow import VERSION_2_0
from pyworkflow.em import Class3D, String
from pyworkflow.em.protocol.protocol import EMProtocol
from pyworkflow.protocol.params import MultiPointerParam
import pyworkflow.utils as pwutils


class XmippProtConsensusClasses3D(EMProtocol):
    """ Compare several SetOfClasses3D.
        Return the intersection of the input classes.
    """
    _label = 'consensus classes 3D'
    _lastUpdateVersion = VERSION_2_0

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMultiClasses', MultiPointerParam, important=True,
                      label="Input Classes", pointerClass='SetOfClasses3D',
                      help='Select several sets of classes where '
                           'to evaluate the intersections.')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        """ Inserting one step for each intersections analisis
        """
        self.intersectsList = []

        self._insertFunctionStep('compareFirstStep', 
                                 self.inputMultiClasses[0].get().getObjId(),
                                 self.inputMultiClasses[1].get().getObjId())

        if len(self.inputMultiClasses) > 2:
            for i in range(2, len(self.inputMultiClasses)):
                self._insertFunctionStep('compareOthersStep', i,
                                     self.inputMultiClasses[i].get().getObjId())

        self._insertFunctionStep('computingDistances')

        self._insertFunctionStep('createOutputStep')

    def compareFirstStep(self, objId1, objId2):
        """ We find the intersections for the two firsts sets of classes
        """
        set1Id = 0
        set2Id = 1
        set1 = self.inputMultiClasses[set1Id].get()
        set2 = self.inputMultiClasses[set2Id].get()

        print('Computing intersections between classes form set %s and set %s:'
               % (set1.getNameId(), set2.getNameId()))

        newList = []
        for cls1 in set1:
            cls1Id = cls1.getObjId()
            ids1 = cls1.getIdSet()

            for cls2 in set2:
                cls2Id = cls2.getObjId()
                ids2 = cls2.getIdSet()

                interTuple = intersectClasses(set1Id, cls1Id, ids1,
                                              set2Id, cls2Id, ids2)

                newList.append(interTuple)

        self.intersectsList = newList

        # saving data just in case of a failure in the next step
        with open(self._getTmpPath('intersection1.pkl'), 'wb') as fileW1:
            pickle.dump(newList, fileW1, pickle.HIGHEST_PROTOCOL)

    def compareOthersStep(self, set1Id, objId):
        """ We find the intersections for the rest of sets of classes
        """
        set1 = self.inputMultiClasses[set1Id].get()

        print('Computing intersections between classes form set %s and '
              'the previous ones:' % (set1.getNameId()))

        currDB = self.intersectsList
        if not currDB:  # rescuing data after a failure/continue
            prevInter = 'intersection%d.pkl' % (set1Id-1)
            with open(self._getTmpPath(prevInter), 'rb') as fileR:
                currDB = pickle.load(fileR)

        newList = []
        for cls1 in set1:
            cls1Id = cls1.getObjId()
            ids1 = cls1.getIdSet()

            for currTuple in currDB:
                ids2 = currTuple[1]
                set2Id = currTuple[2]
                cls2Id = currTuple[3]
                clSize = currTuple[4]
                label  = currTuple[5]

                interTuple = intersectClasses(set1Id, cls1Id, ids1,
                                              set2Id, cls2Id, ids2, clSize, label)
                newList.append(interTuple)
                
        self.intersectsList = newList
        intersectFn = 'intersection%d.pkl' % set1Id
        with open(self._getTmpPath(intersectFn), 'wb') as fileW2:
            pickle.dump(newList, fileW2, pickle.HIGHEST_PROTOCOL)

    def computingDistances(self):

        inputClasses = []
        inClassesLabels = []
        for clsSetId, clsSet in enumerate(self.inputMultiClasses):
            for cls in clsSet.get():
                clsId = cls.getObjId()
                inClassesLabels.append("Set%d-Cls%d" % (clsSetId, clsId))
                idsInClass = cls.getIdSet()
                inputClasses.append(idsInClass)

        currDB = self.intersectsList
        if not currDB:  # rescuing data after a failure/continue
            prevInter = 'intersection%d.pkl' % (len(self.inputMultiClasses)-1)
            with open(self._getTmpPath(prevInter), 'rb') as fileR:
                currDB = pickle.load(fileR)

        nodes = [tup[1] for tup in currDB if tup[0] > 0]
        nodesLabels = [tup[5] for tup in currDB if tup[0] > 0]

        classDistances = getNodeDistances(inputClasses, nodes, range(1, len(nodes)+1))


        # Creating the dendogram. FIXME: Take this and put it in the viewer!!

        from scipy.cluster import hierarchy
        import matplotlib.pyplot as plt
        import numpy as np

        Y = classDistances.values()
        Z = hierarchy.linkage(np.asarray(Y), 'single')
        plt.figure()
        nplabels = np.asarray([x.replace(':', '\n') for x in nodesLabels])
        dn = hierarchy.dendrogram(Z, labels=nplabels)

        plt.style.use("seaborn-whitegrid")
        # plt.title("Dendogram to find clusters")
        plt.ylabel("Distance")
        plt.title(u"Node: cls1\u2229cls2\u2229cls3...")
        plt.savefig(self._getExtraPath("dendogram.png"))


    def createOutputStep(self):

        inputParticles = self.inputMultiClasses[0].get().getImages()
        outputClasses = self._createSetOfClasses3D(inputParticles)

        currDB = self.intersectsList
        if not currDB:  # rescuing data after a failure/continue
            prevInter = 'intersection%d.pkl' % (len(self.inputMultiClasses) - 1)
            with open(self._getTmpPath(prevInter), 'rb') as fileR:
                currDB = pickle.load(fileR)

        for classItem in currDB:
            numOfPart = classItem[0]
            partIds = classItem[1]
            setRepId = classItem[2]
            clsRepId = classItem[3]
            nodeLabel = classItem[5]

            if not partIds:
                continue

            setRep = self.inputMultiClasses[setRepId].get()
            clRep = setRep[clsRepId]

            newClass = Class3D()
            # newClass.copyInfo(clRep)
            newClass.setAcquisition(clRep.getAcquisition())
            newClass.setRepresentative(clRep.getRepresentative())
            newClass.nodeLabel = String(nodeLabel)

            outputClasses.append(newClass)

            enabledClass = outputClasses[newClass.getObjId()]
            enabledClass.enableAppend()
            for itemId in partIds:
                enabledClass.append(inputParticles[itemId].clone())

            outputClasses.update(enabledClass)

        self._defineOutputs(nodeIntersectionClasses=outputClasses)
        for item in self.inputMultiClasses:
            self._defineSourceRelation(item, outputClasses)


    # --------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = [] if len(self.inputMultiClasses)>1 else \
                 ["More than one Input Classes is needed to compute the consensus."]
        return errors

# --------------------------- WORKERS functions ------------------------------
def intersectClasses(setId1, clId1, ids1,
                     setId2, clId2, ids2, clsSize2=None, label=''):
    size1 = len(ids1)
    size2 = len(ids2) if clsSize2 is None else clsSize2

    inter = ids1.intersection(ids2)

    if size1 < size2:
        setId = setId1
        clsId = clId1
        clsSize = size1
    else:
        setId = setId2
        clsId = clId2
        clsSize = size2

    if pwutils.envVarOn('SCIPION_DEBUG'):
        print(" ")
        print(" - Intersection of cl%d of set%d (%d part.) and "
                                 "cl%d of set%d (%d part.):"
               % (clId1, setId1, len(ids1), clId2, setId2, len(ids2)))
        print("    Size1=%d < Size2=%d = %s"
               % (size1, size2, size1<size2))
        print("      -> from set %d calss %d, with %d part. in the intersection."
               % (setId, clsId, len(inter)))
        print(" -  -  -  -  -  -  -  -  -  -")

    if label == '':  # u'\u2229'
        label = '%d:%d' % (clId1, clId2)
    else:
        label += ':%d' % (clId1)

    return len(inter), inter, setId, clsId, clsSize, label


def distanceClassToNode(clsIds, nodeIds, i, k, doPrint=True):

    value = 1 - len(clsIds.intersection(nodeIds)) / float(len(clsIds))
    if doPrint and pwutils.envVarOn('SCIPION_DEBUG'):
        if False:  # Verbose debug
            print("C_%s = %s" % (i, clsIds))
            print("N_%s = %s" % (k, nodeIds))
        print("d(C_%s, N_%s) = 1 - %s / %s = %f" % (i, k,
                len(clsIds.intersection(nodeIds)), float(len(clsIds)), value))

    return value


def getClassDistances(classes, nodes, classLabels=None):
    classDistances = {}

    for i, inClass_i in enumerate(classes):
        for j, inClass_j in enumerate(classes):
            if j <= i:
                continue
            sum_k = 0
            for k, node in enumerate(nodes):
                d_ik = distanceClassToNode(inClass_i, node, i, k, False)
                d_jk = distanceClassToNode(inClass_j, node, j, k, False)
                sum_k += abs(d_ik - d_jk) * abs(d_ik - d_jk) / float(len(node))
            classDistances['%d,%d' % (i, j)] = float(sum_k)  # passing value

    if classLabels:
        print('\nDistances:')
        print(' ' * 23 + '      '.join(classLabels[:0:-1]))
        print('  ' * 9 + '-' * 19 * len(classes))
        for i in range(0, len(classes) - 1):
            line = ["%s   | " % classLabels[i]]
            for j in range(len(classes) - 1, 0, -1):
                if j <= i:
                    continue
                line.append("   %.6f   " % classDistances['%d,%d' % (i, j)])
            print('\t'.join(line))
        print(' ' * len(classLabels[-1] + '\t'.join(classLabels[:])))

    return classDistances


def getNodeDistances(classes, nodes, nodeLabels=None):

    nodeDistances = {}
    for i, node_i in enumerate(nodes):
        for j, node_j in enumerate(nodes):
            if j <= i:
                continue
            sum_k = 0
            for k, cls in enumerate(classes):
                d_ik = distanceClassToNode(cls, node_i, i, k)
                d_jk = distanceClassToNode(cls, node_j, j, k, False)
                sum_k += abs(d_ik - d_jk) * abs(d_ik - d_jk) / float(len(cls))
            nodeDistances['%d,%d' % (i, j)] = float(sum_k)  # passing value

    factor = 1/(float(max(nodeDistances.values())))
    nodeDistances = {k: v*factor for k, v in nodeDistances.iteritems()}

    if nodeLabels and len(nodeLabels) < 30:
        print('\nDistance between nodes:')
        print(' ' * 12 + '      '.join([str(x) for x in nodeLabels[:0:-1]]))
        print(' ' * 6 + '-' * 8 * len(nodes))
        for i in range(0, len(nodes) - 1):
            if i < 9:  # nodeLabel[i]=i+1
                line = ["  %d  | " % nodeLabels[i]]
            else:
                line = ["%d  | " % nodeLabels[i]]
            for j in range(len(nodes) - 1, 0, -1):
                if j <= i:
                    continue
                line.append(" %.2f " % nodeDistances['%d,%d' % (i, j)])
            print(' '.join(line))
        print("")

    return nodeDistances
