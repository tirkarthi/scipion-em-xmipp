# **************************************************************************
# *
# * Authors:    Jose Luis Vilas Prieto (jlvilas@cnb.csic.es)
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

import unittest, sys

from pyworkflow.em import exists
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.em.protocol import ProtImportMicrographs, ProtImportMovies

from xmipp3.protocols import XmippProtSplitOddEven


class TestSplitOddEvenBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='splitoddeven'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micsFn = cls.dataset.getFile('allMics')

    @classmethod
    def runImportMicrographs(cls, pattern, samplingRate):
        """ Run an Import volumes protocol. """
        cls.protImport = cls.newProtocol(ProtImportMicrographs,
                                      filesPath=cls.micsFn,
                                      samplingRate=1.0, voltage=300)
        cls.launchProtocol(cls.protImport)

        return cls.protImport

class TestOddEven(TestSplitOddEvenBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestSplitOddEvenBase.setData()
        cls.protImportMics = cls.runImportMicrographs(cls.micsFn, 1.0)
 
    def testOddEven1(self):
        OddEven = self.newProtocol(XmippProtSplitOddEven,
                                   objLabel = 'split chess',
                                   inputMicrographs = self.protImportMics.outputMicrographs,
                                   typeSplit = 0)
           
        self.launchProtocol(OddEven)
        self.assertTrue(exists(OddEven._getExtraPath('mgresolution.vol')),
                        "Fail in the splitting step")
