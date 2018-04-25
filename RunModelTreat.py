import ParameterClasses as P
import MarkovModel as MarkovCls
import SupportMarkovModel as SupportMarkov
import scr.SamplePathClasses as PathCls
import scr.FigureSupport as Figs


# simulating none therapy
# create a cohort
cohort_none = MarkovCls.Cohort(
    id=0,
    therapy=P.Therapies.NONE)
# simulate the cohort
simOutputs_none = cohort_none.simulate()

# simulating anticoag therapy
# create a cohort
cohort_anticoag = MarkovCls.Cohort(
    id=0,
    therapy=P.Therapies.ANTICOAG)
# simulate the cohort
simOutputs_anticoag = cohort_anticoag.simulate()


# graph survival curve
PathCls.graph_sample_path(
    sample_path=simOutputs_anticoag.get_survival_curve(),
    title='Survival curve',
    x_label='Simulation time step',
    y_label='Number of alive patients'
    )

# graph histogram of survival times
Figs.graph_histogram(
    data=simOutputs_anticoag.get_survival_times(),
    title='Survival times of patients with Stroke',
    x_label='Survival time (years)',
    y_label='Counts',
    bin_width=1
)

# graph histogram of number of strokes
Figs.graph_histogram(
    data=simOutputs_anticoag.get_if_developed_stroke(),
    title='Number of Strokes per Patient',
    x_label='Strokes',
    y_label='Counts',
    bin_width=1
)

# print outcomes (means and CIs)
SupportMarkov.print_outcomes(simOutputs_anticoag, 'Treatment:')

# draw survival curves and histograms
SupportMarkov.draw_survival_curves_and_histograms(simOutputs_none, simOutputs_anticoag)

# print the estimates for the mean survival time and mean time to Stroke
SupportMarkov.print_outcomes(simOutputs_none, "None Therapy:")
SupportMarkov.print_outcomes(simOutputs_anticoag, "Anticoag Therapy:")

# print comparative outcomes
SupportMarkov.print_comparative_outcomes(simOutputs_none, simOutputs_anticoag)

# report the CEA results
SupportMarkov.report_CEA_CBA(simOutputs_none, simOutputs_anticoag)




