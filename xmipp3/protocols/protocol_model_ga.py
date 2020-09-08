# *****************************************************************************
# *
# * Authors:     David Herreros Calero         dherreros@cnb.csic.es
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

import numpy as np
import sys
from joblib import delayed, Parallel
from numba import njit

from pwem.objects import Volume, Float
from pwem.protocols import ProtAnalysis3D
from pwem.emlib.image import ImageHandler

import pyworkflow.protocol.params as params

# bulk = {
#     "A": 11.500,
#     "R": 14.280,
#     "N": 12.820,
#     "D": 11.680,
#     "C": 13.460,
#     "Q": 14.450,
#     "E": 13.570,
#     "G": 3.400,
#     "H": 13.690,
#     "I": 21.400,
#     "L": 21.400,
#     "K": 15.710,
#     "M": 16.250,
#     "F": 19.800,
#     "P": 17.430,
#     "S": 9.470,
#     "T": 15.770,
#     "W": 21.670,
#     "Y": 18.030,
#     "V": 21.570,
# }

bulk = {
    "A": 87.8,
    "R": 192.9,
    "N": 124.7,
    "D": 125.5,
    "C": 105.4,
    "Q": 147.3,
    "E": 13.570,
    "G": 148.0,
    "H": 156.3,
    "I": 166.1,
    "L": 168,
    "K": 184.5,
    "M": 165.2,
    "F": 189.7,
    "P": 123.3,
    "S": 91.7,
    "T": 118.3,
    "W": 227.9,
    "Y": 191.2,
    "V": 138.8,
}

hidroScale = {
    "A": 5.300,
    "R": 4.180,
    "N": 3.710,
    "D": 3.590,
    "C": 7.930,
    "Q": 3.870,
    "E": 3.650,
    "G": 4.480,
    "H": 5.100,
    "I": 8.830,
    "L": 4.470,
    "K": 2.950,
    "M": 8.950,
    "F": 9.030,
    "P": 3.870,
    "S": 4.090,
    "T": 4.490,
    "W": 7.660,
    "Y": 5.890,
    "V": 7.630,
}

@njit
def connectivityIndividual(seqs, individual, num_regions, cMat, idi):
    # score = 0
    score = np.zeros(len(seqs))
    for idx in range(len(seqs)):
        chain_regions = np.where(individual == (idx + 1))
        dMat = dijkstraMatrix(chain_regions[0], cMat)
        score[idx] = connectivityMap(chain_regions[0], dMat)
    return np.sum(score) / (num_regions ** 2), idi

@njit
def dijkstraMatrix(chain_regions, cMat):
    num_regions = len(chain_regions)
    dMat = np.zeros((num_regions, num_regions))
    for idr in range(num_regions):
        row = dijkstra(idr, chain_regions, cMat)
        dMat[idr,:] = row
    return dMat

@njit
def dijkstra(src, chain_regions, cMat):
    num_regions = len(chain_regions)
    dist = [sys.maxsize] * num_regions
    dist[src] = 0
    sptSet = [False] * num_regions

    for cout in range(num_regions):
        u = minDistance(dist, sptSet, num_regions)
        sptSet[u] = True
        for v in range(num_regions):
            if cMat[chain_regions[u], chain_regions[v]] > 0 and \
                    sptSet[v] == False and \
                    dist[v] > dist[u] + cMat[chain_regions[u], chain_regions[v]]:
                dist[v] = dist[u] + cMat[chain_regions[u], chain_regions[v]]
    return np.asarray(dist)

@njit
def minDistance(dist, sptSet, num_regions):
    min = sys.maxsize
    for v in range(num_regions):
        if dist[v] <= min and sptSet[v] == False:
            min = dist[v]
            min_index = v
    return min_index

@njit
def connectivityMap(chain_regions, dMat):
    score = 0
    for idm in range(len(chain_regions)):
        min_dist = sys.maxsize
        max_dist = 0
        for idn in range(len(chain_regions)):
            aux = dMat[idm,idn]
            if idm != idn and aux < min_dist:
                min_dist = aux
            elif idm != idn and aux > max_dist:
                max_dist = aux
        score += min_dist + max_dist
    return score

@njit
def massIndividual(individual, seqs, submass, map_region_mass, idi):
    score = np.zeros(len(seqs))
    for idx in range((len(seqs))):
        map_regions = np.where(individual == (idx + 1))
        score[idx] = np.abs(submass[idx] - np.sum(map_region_mass[map_regions]))
    # return np.sum(score)/len(score), idi
    return np.amax(score), idi

@njit
def connectivityChains(seqs, individual, num_regions, cMat, idi):
    score = np.ones(len(seqs) - 1)
    for idx in range((len(seqs)) - 1):
        firstChain = np.where(individual == (idx + 1))
        secondChain = np.where(individual == (idx + 2))
        if len(firstChain[0]) > 0 and len(secondChain[0]) > 0:
            aux = 0
            count = 0
            for idm in firstChain[0]:
                for idn in secondChain[0]:
                    # if cMat[idm, idn] < interRegionScore:
                    if cMat[idm, idn] != num_regions:
                        aux += cMat[idm, idn]
                        count += 1
            # # Comentar el if para buscar las zonas sin ninguna conexion
            # if count == 0:
            #     aux = sys.maxsize
            # else:
            #     aux /= count * num_regions
        else:
            aux = sys.maxsize
        score[idx] = aux
    return np.sum(score) / (len(score)), idi

@njit
def connectivityChains2(seqs, individual, num_regions, cMat, dMat, idi):
    cc_mat = np.asarray([[-1, 1, 1, 0], [1, -1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1]])
    score = np.ones(len(seqs) - 1)
    for idr in range(len(cc_mat)):
        scores_row = []
        row = cc_mat[idr]
        for idc in range(len(row)):
            # if row[idc] == 1:
                firstChain = np.where(individual == (idr + 1))
                secondChain = np.where(individual == (idc + 1))
                if len(firstChain[0]) > 0 and len(secondChain[0]) > 0:
                    aux = 0
                    count = 0
                    for idm in firstChain[0]:
                        for idn in secondChain[0]:
                            if row[idc] == 1:
                                #if cMat[idm, idn] != num_regions:
                                    aux += cMat[idm, idn]
                                    count += 1
                            else:
                                #if dMat[idm, idn] != num_regions:
                                    aux += dMat[idm, idn]
                                    count += 1
                    if count == 0:
                        aux = 0
                    else:
                        aux /= 10*count * num_regions
                else:
                    aux = sys.maxsize
                scores_row.append(aux)
        score[idr] = max(scores_row)
    return np.sum(score) / (len(score)), idi

@njit
def hidrophobicityIndividual(chainHidro, individual, regionsContactBg, idi):
    score = np.ones(len(chainHidro))
    for idx in range(len(chainHidro)):
        chainRegions = np.where(individual == (idx+1))
        score[idx] = np.abs(chainHidro[idx] - np.sum(regionsContactBg[chainRegions]))
    return np.amax(score), idi


class XmippProtModelGA(ProtAnalysis3D):
    """Modeling implemented through genetic algorithm"""
    _label = 'genetic algortihm modeling'

    # --------------------------- DEFINE param functions ---------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMask', params.PointerParam, pointerClass='Volume', label='Input Masks',
                      important=True, help='Input mask with identifiers (not binary).')
        form.addParam('inputSeqs', params.MultiPointerParam, pointerClass='Sequence', label='Input Sequence', important=True,
                      help='Secuences to be modeled in the EM map (it can be a protein a chain...)')
        form.addParam('population', params.IntParam, label='Population Size', default=100, expertlevel=params.LEVEL_ADVANCED,
                      help='Number of indivuals in the population for the genetic algorithm')
        form.addParam('generations', params.IntParam, label='Number of Generations', default=5, expertlevel=params.LEVEL_ADVANCED,
                      help='Number of generations to be computed by the genetic algorithm')
        form.addParam('parents', params.IntParam, label='Number of Parents', default=5, expertlevel=params.LEVEL_ADVANCED,
                      help='Number of parents to be mated')
        form.addParam('p_mutation', params.FloatParam, label='Mutation Probability', default=0.3, expertlevel=params.LEVEL_ADVANCED,
                      help='Probability of introducing a mutation in an individual from the offspring')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- Steps functions ---------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('geneticAlgorithm')
        self._insertFunctionStep('createOutputStep')

    def geneticAlgorithm(self):
        ih = ImageHandler()
        self.seqs = [seq.get().getSequence() for seq in self.inputSeqs]
        self.seqs = np.asarray(self.seqs)
        self.idMask = ih.read(self.inputMask.get().getFileName()).getData()
        self.idMask = np.squeeze(self.idMask)

        sol_per_population = self.population.get()
        self.regions_id = np.unique(np.reshape(self.idMask, -1))
        self.regions_id = np.delete(self.regions_id, 0)
        self.num_regions = len(self.regions_id)
        pop_size = (sol_per_population, self.num_regions)
        new_population = np.random.random_integers(low=1, high=len(self.seqs), size=pop_size)
        num_generations = self.generations.get()
        num_parents = self.parents.get()
        self.cMat, self.dMat, self.regionsContactBg = self.connectivityMatrix()
        # print(self.regionsContactBg)

        mean_density_prot = 8.1325e-04  # KDa / (A^3)
        mean_mass_aa = 0.110  # KDa
        sampling_rate = self.inputMask.get().getSamplingRate() ** 3  # A^3 / voxel
        mean_density_prot *= sampling_rate

        self.submass = [mean_mass_aa * len(subseq) for subseq in self.seqs]
        print(sum([mean_mass_aa * len(subseq) for subseq in self.seqs]))
        # self.submass = []
        # for seq in self.seqs:
        #     vol = 0
        #     for aa in seq:
        #         vol += bulk[aa]
        #     self.submass.append(vol)
        # print(self.submass)
        self.submass = np.asarray(self.submass)
        # self.submass /= (mean_density_prot / sampling_rate)

        self.map_region_mass = [mean_density_prot * np.sum(self.idMask == idr) for idr in self.regions_id]
        self.map_region_mass = np.asarray(self.map_region_mass)
        print(mean_density_prot * np.sum(self.idMask != 0))
        # factor = np.sum(self.submass) / np.sum(self.map_region_mass)
        # self.map_region_mass *= factor

        self.map_region_mass /= np.sum(self.map_region_mass)
        self.submass /= np.sum(self.submass)

        self.chainHidro = []
        for seq in self.seqs:
            hidro = 0
            for aa in seq:
                hidro += hidroScale[aa]
            self.chainHidro.append(hidro)
        self.chainHidro = np.asarray(self.chainHidro)
        self.chainHidro /= np.sum(self.chainHidro)
        # print(self.chainHidro)

        for generation in range(num_generations):
            print('Generation: ', (generation+1))
            score_population = self.massScore(new_population)
            score_population += self.connectivityScore(new_population)
            score_population += self.connectivityScoreChain(new_population)
            score_population += self.hidrophobicityScore(new_population)
            parents = self.matingPool(new_population, score_population, num_parents)
            offspring_size = (pop_size[0] - parents.shape[0], self.num_regions)
            offspring_crossover = self.crossover(new_population, offspring_size)
            offspring_mutation = self.mutation(offspring_crossover)
            new_population[0:parents.shape[0], :] = parents
            new_population[parents.shape[0]:, :] = offspring_mutation

            # FIXME: Probably this can be removed
            score_population = self.massScore(new_population)
            score_population += self.connectivityScore(new_population)
            score_population += self.connectivityScoreChain(new_population)
            score_population += self.hidrophobicityScore(new_population)
            print('Best result after generation %d: %f' % ((generation+1), np.amin(score_population)))
            sys.stdout.flush()

        best_individuals = np.argsort(score_population)
        print(new_population[best_individuals[0]])
        idx = np.round(np.linspace(0, num_parents, 20)).astype(int)
        self.bestIndividuals = new_population[best_individuals[idx]]
        self.bestScores = score_population[best_individuals[idx]]

    def createOutputStep(self):
        ih = ImageHandler()
        outputMasks = self._createSetOfVolumes()
        outputMasks.setSamplingRate(self.inputMask.get().getSamplingRate())
        score_id = 1
        for idi, individual in enumerate(self.bestIndividuals):
            outMask = ih.createImage()
            outData = np.zeros(self.idMask.shape, float)
            for pos, idm in enumerate(self.regions_id):
                logic_mask = self.idMask == idm
                outData += individual[pos] * (self.idMask * logic_mask / idm)
            outMask.setData(outData)
            ih.write(outMask, self._getExtraPath('outMask_%d.mrc' % score_id))
            volume = Volume()
            volume.setSamplingRate(self.inputMask.get().getSamplingRate())
            volume.setLocation(self._getExtraPath('outMask_%d.mrc' % score_id))
            volume.score = Float(self.bestScores[idi])
            outputMasks.append(volume)
            score_id += 1
        self._defineOutputs(outputMasks=outputMasks)
        self._defineSourceRelation(self.inputMask, outputMasks)


    # --------------------------- Utils functions ----------------------
    def massScore(self, population):
        score_population = np.zeros(len(population))
        out = Parallel(n_jobs=self.numberOfThreads.get()) \
            (delayed(massIndividual)(individual, self.seqs, self.submass, self.map_region_mass, idi)
             for idi, individual in enumerate(population))

        for score, pos in out:
            score_population[pos] = score

        return score_population

    def connectivityScore(self, population):
        score_population = np.zeros(len(population))
        # for idi, individual in enumerate(population):
        out = Parallel(n_jobs=self.numberOfThreads.get())\
            (delayed(connectivityIndividual)(self.seqs, individual, self.num_regions, self.cMat, idi)
             for idi, individual in enumerate(population))

        for score, pos in out:
            score_population[pos] = score

        return score_population

    def connectivityScoreChain(self, population):
        score_population = np.zeros(len(population))
        # for idi, individual in enumerate(population):
        out = Parallel(n_jobs=self.numberOfThreads.get())\
            (delayed(connectivityChains2)(self.seqs, individual, self.num_regions, self.cMat, self.dMat, idi)
             for idi, individual in enumerate(population))

        for score, pos in out:
            score_population[pos] = score

        return score_population

    def hidrophobicityScore(self, population):
        score_population = np.zeros(len(population))
        # for idi, individual in enumerate(population):
        out = Parallel(n_jobs=self.numberOfThreads.get())\
            (delayed(hidrophobicityIndividual)(self.chainHidro, individual, self.regionsContactBg, idi)
             for idi, individual in enumerate(population))

        for score, pos in out:
            score_population[pos] = score

        return score_population

    # def connectivityMap(self, chain_regions, dMat):
    #     score = 0
    #     for idm in range(len(chain_regions)):
    #         min_dist = sys.maxsize
    #         max_dist = 0
    #         for idn in range(len(chain_regions)):
    #             aux = dMat[idm,idn]
    #             if idm != idn and aux < min_dist:
    #                 min_dist = aux
    #             elif idm != idn and aux > max_dist:
    #                 max_dist = aux
    #         score += min_dist + max_dist
    #
    #     # score = 0
    #     # size = len(chain_regions)
    #     # connected = np.zeros(size)
    #     # for idx in range(size - 1):
    #     #     min_dist = np.inf
    #     #     for idy in range(size):
    #     #         aux = dMat[idx,idy]
    #     #         if idx != idy and aux < min_dist and connected[idy] == 0:
    #     #             min_dist = aux
    #     #     score += min_dist
    #     #     connected[idx] = 1
    #
    #     return score

    def matingPool(self, population, score, num_parents):
        parents = np.empty((num_parents, population.shape[1]))

        for idp in range(num_parents):
            parents[idp, :] = population[np.argmin(score), :]
            score[np.argmin(score)] = np.inf

        return parents

    def crossover(self, population, offspring_size):
        offspring = np.empty(offspring_size)

        for k in range(offspring_size[0]):
            crossover_point = np.random.random_integers(0, offspring_size[1] - 1)
            parent1_idx = np.random.randint(0, population.shape[0])
            parent2_idx = np.random.randint(0, population.shape[0])
            offspring[k, 0:crossover_point] = population[parent1_idx, 0:crossover_point]
            offspring[k, crossover_point:] = population[parent2_idx, crossover_point:]

        return offspring

    def mutation(self, offspring):
        p_mutation = self.p_mutation.get()
        for idx in range(offspring.shape[0]):
            num_rnd = np.random.uniform()
            idg_rnd = np.random.random_integers(0, offspring.shape[1] - 1)
            chain_rnd = np.random.random_integers(1, len(self.seqs))
            if num_rnd <= p_mutation:
                offspring[idx, idg_rnd] = float(chain_rnd)

        return offspring

    def connectivityMatrix(self):
        regionsContactBg = np.zeros(self.num_regions)
        voxelsRegion = np.zeros(self.num_regions)
        cMat = np.zeros((self.num_regions, self.num_regions))
        for idr in range(self.num_regions):
            row, boundary_voxels, contactBgVoxels = self.neighbours(self.regions_id[idr])
            voxelsRegion[idr] = boundary_voxels
            regionsContactBg[idr] = contactBgVoxels
            cMat[idr] = row

        # Connectivity matrix normalization (Dice coefficient)
        for idm in range(self.num_regions):
            for idn in range(self.num_regions):
                sumVoxels = voxelsRegion[idm] + voxelsRegion[idn]
                cMat[idm,idn] = 1 - (2 * cMat[idm,idn] / sumVoxels)
        dMat = 1 - cMat
        dMat[dMat == 1] = self.num_regions
        cMat[cMat == 1] = self.num_regions
        return cMat, dMat, regionsContactBg / np.sum(regionsContactBg)

    def neighbours(self, region_id):
        boundary_voxels = 0
        contactBgVoxels = 0
        voxels = np.asarray(np.where(self.idMask == region_id))
        row = np.zeros(self.num_regions)
        for idv in range(voxels.shape[1]):
            coords = voxels[:,idv]
            submat = self.idMask[coords[0]-1:coords[0]+2, coords[1]-1:coords[1]+2, coords[2]-1:coords[2]+2]
            submat = submat.reshape(-1)
            touchesBg = True if np.sum(submat == 0) > 0 else False
            submat = np.unique(submat)
            submat = submat[(submat != 0) * (submat != region_id)]
            if len(submat) > 0:
                row[submat.astype('int') - 1] += 1
                boundary_voxels += 1
            if touchesBg:
                contactBgVoxels += 1
        return row, boundary_voxels, contactBgVoxels

    # def dijkstraMatrix(self, chain_regions):
    #     num_regions = len(chain_regions)
    #     dMat = np.zeros((num_regions, num_regions))
    #     for idr in range(num_regions):
    #         row = self.dijkstra(idr, chain_regions)
    #         dMat[idr,:] = row
    #     return dMat
    #
    # def dijkstra(self, src, chain_regions):
    #     num_regions = len(chain_regions)
    #     dist = [sys.maxsize] * num_regions
    #     dist[src] = 0
    #     sptSet = [False] * num_regions
    #
    #     for cout in range(num_regions):
    #         u = self.minDistance(dist, sptSet, num_regions)
    #         sptSet[u] = True
    #         for v in range(num_regions):
    #             if self.cMat[chain_regions[u], chain_regions[v]] > 0 and \
    #                     sptSet[v] == False and \
    #                     dist[v] > dist[u] + self.cMat[chain_regions[u], chain_regions[v]]:
    #                 dist[v] = dist[u] + self.cMat[chain_regions[u], chain_regions[v]]
    #     return np.asarray(dist)
    #
    # def minDistance(self, dist, sptSet, num_regions):
    #     min = sys.maxsize
    #     for v in range(num_regions):
    #         if dist[v] <= min and sptSet[v] == False:
    #             min = dist[v]
    #             min_index = v
    #     return min_index

    # --------------------------- Info functions ----------------------
    def _methods(self):
        pass

    def _summary(self):
        summary = []
        if self.getOutputsSize() >= 1:
            summary.append("Number of individuals saved: %d" % self.outputMasks.getSize())
        else:
            summary.append("Output masks not ready yet.")
        return summary
