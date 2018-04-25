POP_SIZE = 2000     # cohort population size
SIM_LENGTH = 15   # length of simulation (years)
ALPHA = 0.05        # significance level for calculating confidence intervals
DELTA_T = 1         # years (length of time step, how frequently you look at the patient)
DISCOUNT = 0.03     # annual discount rate

ADD_BACKGROUND_MORT = True  # if background mortality should be added
DELTA_T = 1/4       # years
PSA_ON = False      # if probabilistic sensitivity analysis is on


# transition matrix
TRANS_MATRIX = [
    [0.75,  0.15,   0.0,    0.1],   # Well
    [0,     0.0,    0.90,   0.10],   # Stroke
    [0,     0.17,   0.80,   0.03],   # Post-Stroke
    [0.0,   0.0,    0.0,    1.0],   # Dead
    ]

# annual cost of each health state
ANNUAL_STATE_COST = [
    0000.0,   # Well
    5000.0,   # Stroke
    200.0,    # Post-Stroke
    0000.0    # Dead
    ]

# annual health utility of each health state
ANNUAL_STATE_UTILITY = [
    1.00,   # Well
    0.20, # Stroke
    0.90,   # Post-Stroke
    0.00    # Dead
    ]

# annual drug costs
ANTICOAG_COST = 2000.0

# anticoagulation relative risk in reducing stroke incidence and stroke death while in “Post-Stroke”
RR_STROKE = 0.65
# anticoagulation relative risk in increasing mortality due to bleeding is 1.05.
RR_BLEEDING = 1.05

# annual probability of background mortality (number per year per 1,000 population)
ANNUAL_PROB_BACKGROUND_MORT = 17.638 / 1000
