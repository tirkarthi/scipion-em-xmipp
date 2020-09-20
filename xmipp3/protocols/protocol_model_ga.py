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


def contactVoxels(idMask):
    num_regions = int(np.amax(idMask))
    cc_mat = np.zeros((num_regions, num_regions))
    bgVoxels = np.zeros(num_regions)
    for rid in range(1, num_regions + 1):
        row, boundary_voxels = neighbours(rid, idMask, num_regions)
        bgVoxels[rid-1] = boundary_voxels
        cc_mat[rid-1] = row
    return cc_mat, bgVoxels

def diceBoundaryMatrix(cc_mat, bgVoxels):
    num_regions = len(bgVoxels)
    dc_mat = np.zeros((num_regions, num_regions))
    for idm in range(num_regions):
        for idn in range(num_regions):
            # sumVoxels = bgVoxels[idm] + bgVoxels[idn] + np.sum(cc_mat[idm]) + np.sum(cc_mat[idn])
            sumVoxels = np.sum(cc_mat[idm]) + np.sum(cc_mat[idn])
            dc_mat[idm,idn] = 2 * cc_mat[idm,idn] / sumVoxels
    return dc_mat

@njit
def neighbours(region_id, idMask, num_regions):
    bgVoxels = 0
    voxels = np.where(idMask == region_id)
    row = np.zeros(num_regions)
    for idv in range(len(voxels[0])):
        coords = [voxels[0][idv], voxels[1][idv], voxels[2][idv]]
        submat = idMask[coords[0]-1:coords[0]+2, coords[1]-1:coords[1]+2, coords[2]-1:coords[2]+2]
        submat = submat.copy().reshape(-1)
        # submat = np.unique(submat)
        submat = submat[(submat != region_id)].astype(np.int32)
        if len(submat) > 0:
            for rid in submat:
                if rid == 0:
                    bgVoxels += 1
                else:
                    row[rid - 1] += 1
    return row, bgVoxels

def costCCChains(cc_chain, max_cc, cc_mat, bgVoxels, numChains, individual, idi):
    # limit = 0.5 * max_cc
    limit = 0.1
    if len(np.unique(individual)) == numChains:
        chain_cc_mat, chainBgVoxels = computeChains(cc_mat, bgVoxels, numChains, individual)
        dbc_chains = diceBoundaryMatrix(chain_cc_mat, chainBgVoxels)
        dbc_chains[dbc_chains >= limit] = 1
        dbc_chains[dbc_chains < limit] = 0
        score = np.sum(np.abs(cc_chain - dbc_chains)) / (2 * chain_cc_mat.size)
    else:
        score = sys.maxsize
    return score, idi

def computeChains(cc_mat, bgVoxels, numChains, individual):
    chain_cc_mat = np.zeros((numChains, numChains))
    chainBgVoxels = np.zeros(numChains)
    for idm in range(numChains):
        cid_1 = idm + 1
        cc_idm = np.where(individual == cid_1)
        for idn in range(idm+1, numChains):
            cid_2 = idn + 1
            cc_idn = np.where(individual == cid_2)
            chain_cc_mat[idm, idn] += np.sum(cc_mat[cc_idm][:, cc_idn]) / 2
        chainBgVoxels[idm] += np.sum(bgVoxels[cc_idm])
    return chain_cc_mat + chain_cc_mat.T, chainBgVoxels

@njit
def connectivityIndividual(seqs, individual, num_regions, cMat, idi):
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

# @njit
def massIndividual(individual, submass, map_region_mass, idi):
    score = np.zeros(len(submass))
    for idx in range((len(submass))):
        map_regions = np.where(individual == (idx + 1))
        score[idx] = np.abs((submass[idx] - np.sum(map_region_mass[map_regions])) / submass[idx])
    return np.sum(score) / len(score), idi

# @njit
# def connectivityChains(seqs, individual, num_regions, cMat, idi):
#     score = np.ones(len(seqs) - 1)
#     for idx in range((len(seqs)) - 1):
#         firstChain = np.where(individual == (idx + 1))
#         secondChain = np.where(individual == (idx + 2))
#         if len(firstChain[0]) > 0 and len(secondChain[0]) > 0:
#             aux = 0
#             count = 0
#             for idm in firstChain[0]:
#                 for idn in secondChain[0]:
#                     # if cMat[idm, idn] < interRegionScore:
#                     if cMat[idm, idn] != num_regions:
#                         aux += cMat[idm, idn]
#                         count += 1
#             # # Comentar el if para buscar las zonas sin ninguna conexion
#             # if count == 0:
#             #     aux = sys.maxsize
#             # else:
#             #     aux /= count * num_regions
#         else:
#             aux = sys.maxsize
#         score[idx] = aux
#     return np.sum(score) / (len(score)), idi

# @njit
# def connectivityChains2(seqs, individual, cMat, dMat, idi):
#     cc_mat = np.asarray([[-1, 1, 1, 0], [1, -1, 0, 1], [1, 0, -1, 0], [0, 1, 0, -1]])
#     score = np.ones(len(seqs) - 1)
#     for idr in range(len(cc_mat)):
#         scores_row = []
#         row = cc_mat[idr]
#         for idc in range(len(row)):
#             firstChain = np.where(individual == (idr + 1))
#             secondChain = np.where(individual == (idc + 1))
#             if len(firstChain[0]) > 0 and len(secondChain[0]) > 0:
#                 score_total = 0
#                 for idm in firstChain[0]:
#                     aux = sys.maxsize
#                     count = 0
#                     for idn in secondChain[0]:
#                         count += 1
#                         if row[idc] == 1:
#                             if cMat[idm, idn] <= 0.5:
#                                 aux = 0
#                             elif aux > cMat[idm, idn]:
#                                 aux = cMat[idm, idn]
#                         elif row[idc] == 0:
#                             if aux > dMat[idm, idn]:
#                                 aux = dMat[idm, idn]
#                         else:
#                             aux = 0
#                     score_total += aux / count
#                 scores_row.append(score_total)
#             else:
#                 aux = sys.maxsize
#                 scores_row.append(aux)
#         scores_row = np.asarray(scores_row, dtype=np.float64)
#         score[idr] = np.sum(scores_row) / len(scores_row)
#     return np.sum(score) / (len(score) * 40), idi


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
        self.cc_mat, self.bgVoxels = contactVoxels(self.idMask)
        self.cMat = diceBoundaryMatrix(self.cc_mat, self.bgVoxels)
        self.cMat = 1 - self.cMat
        self.cMat[self.cMat == 1] = self.num_regions
        # self.cMat, self.dMat, self.regionsContactBg = self.connectivityMatrix()
        # print(self.regionsContactBg)

        mean_density_prot = 8.1325e-04  # KDa / (A^3)
        mean_mass_aa = 0.110  # KDa
        sampling_rate = self.inputMask.get().getSamplingRate() ** 3  # A^3 / voxel
        mean_density_prot *= sampling_rate

        self.submass = [mean_mass_aa * len(subseq) for subseq in self.seqs]
        print(sum(self.submass))
        self.submass = np.asarray(self.submass)

        self.map_region_mass = [mean_density_prot * np.sum(self.idMask == rid) for rid in self.regions_id]
        self.map_region_mass = np.asarray(self.map_region_mass)
        print(np.sum(self.map_region_mass))
        # factor = np.sum(self.submass) / np.sum(self.map_region_mass)
        # self.map_region_mass *= factor

        for generation in range(num_generations):
            print('Generation: ', (generation+1))
            score_population = self.massScore(new_population)
            score_population += self.connectivityScore(new_population)
            score_population += self.connectivityScoreChain(new_population)
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
            print('Best result after generation %d: %f' % ((generation+1), np.amin(score_population)))
            sys.stdout.flush()

        best_individuals = np.argsort(score_population)
        print(new_population[best_individuals[0]])
        idx = np.round(np.linspace(0, num_parents, 20)).astype(int)
        # idx = np.round(np.linspace(0, sol_per_population-1, 20)).astype(int)
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
            (delayed(massIndividual)(individual, self.submass, self.map_region_mass, idi)
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
        max_cc = np.amax(self.cc_mat)
        cc_chain = np.asarray([[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
        # for idi, individual in enumerate(population):
        #     out = costCCChains(individual, cc_chain, self.idMask, idi)
        out = Parallel(n_jobs=self.numberOfThreads.get())\
            (delayed(costCCChains)(cc_chain, max_cc, self.cc_mat, self.bgVoxels, len(self.seqs), individual, idi)
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

    # def connectivityMatrix(self):
    #     regionsContactBg = np.zeros(self.num_regions)
    #     voxelsRegion = np.zeros(self.num_regions)
    #     cMat = np.zeros((self.num_regions, self.num_regions))
    #     for idr in range(self.num_regions):
    #         row, boundary_voxels, contactBgVoxels = self.neighbours(self.regions_id[idr])
    #         voxelsRegion[idr] = boundary_voxels
    #         regionsContactBg[idr] = contactBgVoxels
    #         cMat[idr] = row
    #
    #     # Connectivity matrix normalization (Dice coefficient)
    #     for idm in range(self.num_regions):
    #         for idn in range(self.num_regions):
    #             sumVoxels = voxelsRegion[idm] + voxelsRegion[idn]
    #             cMat[idm,idn] = 1 - (2 * cMat[idm,idn] / sumVoxels)
    #     dMat = 1 - cMat
    #     dMat[dMat == 1] = self.num_regions
    #     cMat[cMat == 1] = self.num_regions
    #     return cMat, dMat, regionsContactBg / np.sum(regionsContactBg)
    #
    # def neighbours(self, region_id):
    #     boundary_voxels = 0
    #     contactBgVoxels = 0
    #     voxels = np.asarray(np.where(self.idMask == region_id))
    #     row = np.zeros(self.num_regions)
    #     for idv in range(voxels.shape[1]):
    #         coords = voxels[:,idv]
    #         submat = self.idMask[coords[0]-1:coords[0]+2, coords[1]-1:coords[1]+2, coords[2]-1:coords[2]+2]
    #         submat = submat.reshape(-1)
    #         touchesBg = True if np.sum(submat == 0) > 0 else False
    #         submat = np.unique(submat)
    #         submat = submat[(submat != 0) * (submat != region_id)]
    #         if len(submat) > 0:
    #             row[submat.astype('int') - 1] += 1
    #             boundary_voxels += 1
    #         if touchesBg:
    #             contactBgVoxels += 1
    #     return row, boundary_voxels, contactBgVoxels

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
