from enum import Enum
import InputData as Data
import numpy as np
import scipy.stats as stat
import math as math
import InputData as Data
import scr.MarkovClasses as MarkovCls
import scr.RandomVariantGenerators as Random
import scr.FittingProbDist_MM as Est

class HealthStats(Enum):
    """ health states of patients with HIV """
    WELL = 0
    STROKE = 1
    POST_STROKE = 2
    DEATH = 3
    BACKGROUND_DEATH = 4

class Therapies(Enum):
    """ mono vs. combination therapy """
    NONE = 0
    ANTICOAG = 1


class ParametersFixed():
    def __init__(self, therapy):

        # selected therapy
        self._therapy = therapy

        # simulation time step
        self._delta_t = Data.DELTA_T

        # initial health state
        self._initialHealthState = HealthStats.WELL

        # transition probability matrix of the selected therapy
        self._prob_matrix = []
        # treatment relative risk
        self._treatmentRR = 0

        # calculate transition probabilities depending of which therapy options is in use
        if therapy == Therapies.NONE:
            self._prob_matrix = Data.TRANS_MATRIX
        else:
            self._prob_matrix = calculate_prob_matrix_anticoag()

    def get_initial_health_state(self):
        return self._initialHealthState

    def get_delta_t(self):
        return self._delta_t

    def get_transition_prob(self, state):
        return self._prob_matrix[state.value]

    def get_annual_state_cost(self, state):
        if state == HealthStats.STROKE_DEATH or state == HealthStats.BACKGROUND_DEATH:
            return 0
        else:
            return self._annualStateCosts[state.value]

    def get_annual_state_utility(self, state):
        if state == HealthStats.STROKE_DEATH or state == HealthStats.BACKGROUND_DEATH:
            return 0
        else:
            return self._annualStateUtilities[state.value]

    def get_annual_treatment_cost(self):
        return self._annualTreatmentCost

class ParametersFixed(_Parameters):
    def __init__(self, therapy):

        # initialize the base class
        _Parameters.__init__(self, therapy)

        # calculate transition probabilities between stroke states
        self._prob_matrix = calculate_prob_matrix_NONE()
        # add background mortality if needed
        if Data.ADD_BACKGROUND_MORT:
            add_background_mortality(self._prob_matrix)

        # update the transition probability matrix if combination therapy is being used
        if self._therapy == Therapies.ANTICOAG:
            # treatment relative risk
            self._treatmentRR = Data.TREATMENT_RR
            # calculate transition probability matrix for the combination therapy
            self._prob_matrix = calculate_prob_matrix_ANTICOAG(
                matrix_none=self._prob_matrix, anticoag_rr=Data.TREATMENT_RR)

        # annual state costs and utilities
        self._annualStateCosts = Data.ANNUAL_STATE_COST
        self._annualStateUtilities = Data.ANNUAL_STATE_UTILITY

class ParametersProbabilistic(_Parameters):
    def __init__(self, seed, therapy):

        # initializing the base class
        _Parameters.__init__(self, therapy)

        self._rng = Random.RNG(seed)    # random number generator to sample from parameter distributions
        self._strokeProbMatrixRVG = []  # list of dirichlet distributions for transition probabilities
        self._lnRelativeRiskRVG = None  # random variate generator for the natural log of the treatment relative risk
        self._annualStateCostRVG = []       # list of random variate generators for the annual cost of states
        self._annualStateUtilityRVG = []    # list of random variate generators for the annual utility of states

        # HIV transition probabilities
        j = 0
        for prob in Data.TRANS_MATRIX:
            self._strokeProbMatrixRVG.append(Random.Dirichlet(prob[j:]))
            j += 1

        # treatment relative risk
        # find the mean and st_dev of the normal distribution assumed for ln(RR)
        sample_mean_lnRR = math.log(Data.TREATMENT_RR)
        sample_std_lnRR = \
            (math.log(Data.TREATMENT_RR_CI[1])-math.log(Data.TREATMENT_RR_CI[0]))/(2*stat.norm.ppf(1-0.05/2))
        self._lnRelativeRiskRVG = Random.Normal(loc=sample_mean_lnRR, scale=sample_std_lnRR)

        # annual state cost
        for cost in Data.ANNUAL_STATE_COST:
            # find shape and scale of the assumed gamma distribution
            estDic = Est.get_gamma_params(mean=cost, st_dev=cost/4)
            # append the distribution
            self._annualStateCostRVG.append(
                Random.Gamma(a=estDic["a"], loc=0, scale=estDic["scale"]))

        # annual state utility
        for utility in Data.ANNUAL_STATE_UTILITY:
            # find alpha and beta of the assumed beta distribution
            estDic = Est.get_beta_params(mean=utility, st_dev=utility/4)
            # append the distribution
            self._annualStateUtilityRVG.append(
                Random.Beta(a=estDic["a"], b=estDic["b"]))

        # resample parameters
        self.__resample()

    def __resample(self):

        # calculate transition probabilities
        # create an empty matrix populated with zeroes
        self._prob_matrix = []
        for s in HealthStats:
            self._prob_matrix.append([0] * len(HealthStats))

        # for all health states
        for s in HealthStats:
            # if the current state is death
            if s in [HealthStats.STROKE_DEATH, HealthStats.BACKGROUND_DEATH]:
                # the probability of staying in this state is 1
                self._prob_matrix[s.value][s.value] = 1
            else:
                # sample from the dirichlet distribution to find the transition probabilities between hiv states
                sample = self._strokeProbMatrixRVG[s.value].sample(self._rng)
                for j in range(len(sample)):
                    self._prob_matrix[s.value][s.value+j] = sample[j]

        # add background mortality if needed
        if Data.ADD_BACKGROUND_MORT:
            add_background_mortality(self._prob_matrix)

        # update the transition probability matrix if combination therapy is being used
        if self._therapy == Therapies.ANTICOAG:
            # treatment relative risk
            self._treatmentRR = math.exp(self._lnRelativeRiskRVG.sample(self._rng))
            # calculate transition probability matrix for the combination therapy
            self._prob_matrix = calculate_prob_matrix_anticoag(
                matrix_none=self._prob_matrix, anticoag_rr=self._treatmentRR)

        # sample from gamma distributions that are assumed for annual state costs
        self._annualStateCosts = []
        for dist in self._annualStateCostRVG:
            self._annualStateCosts.append(dist.sample(self._rng))

        # sample from beta distributions that are assumed for annual state utilities
        self._annualStateUtilities = []
        for dist in self._annualStateUtilityRVG:
            self._annualStateUtilities.append(dist.sample(self._rng))

def calculate_prob_matrix_anticoag():
    """ :returns transition probability matrix under anticoagulation use"""

    # create an empty matrix populated with zeroes
    prob_matrix = []
    for s in HealthStats:
        prob_matrix.append([0] * len(HealthStats))

    # for all health states
    for s in HealthStats:
        # if the current state is post-stroke
        if s == HealthStats.POST_STROKE:
            # post-stoke to stroke
            prob_matrix[s.value][HealthStats.STROKE.value]\
                = Data.RR_STROKE*Data.TRANS_MATRIX[s.value][HealthStats.STROKE.value]
            # post-stroke to death
            prob_matrix[s.value][HealthStats.DEATH.value] \
                = Data.RR_STROKE * Data .RR_BLEEDING * Data.TRANS_MATRIX[s.value][HealthStats.DEATH.value]
            # staying in post-stroke
            prob_matrix[s.value][s.value]\
                = 1 -prob_matrix[s.value][HealthStats.STROKE.value] -prob_matrix[s.value][HealthStats.DEATH.value]
        else:
            prob_matrix[s.value] = Data.TRANS_MATRIX[s.value]

    return prob_matrix
