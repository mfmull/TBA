# ------------------------------------------------------------------------------
# Script purpose (preamble)
# ------------------------------------------------------------------------------
# This script produces standard diagnostics for the matching step used in the
# second-stage analysis. It loads a previously saved MatchIt object and generates
# (i) a covariate balance plot (Love plot) with an absolute standardized-mean-
# difference threshold and (ii) a histogram diagnostic of the propensity score /
# distance measure from the MatchIt output.
#
# Inputs
# - ../3_secondstage/match.rds
#     A saved MatchIt object (from matchit()), containing the matching
#     specification, weights, and balance information for TB ~ covariates.
#
# Diagnostics produced
# 1) Covariate balance (cobalt::love.plot):
#    - Displays absolute standardized mean differences before and after matching.
#    - Uses a conventional balance threshold of 0.1 (thresholds = c(m = 0.1)),
#      allowing quick identification of covariates that remain imbalanced.
#    - abs = TRUE plots absolute SMDs (direction suppressed) to focus on magnitude.
# 2) Matching distribution check (plot(m.out, type = "hist")):
#    - Produces MatchIt’s built-in histogram diagnostic for the “distance” used
#      in matching (typically the propensity score or linear predictor, depending
#      on the method/settings).
#    - Useful to assess overlap/common support between TB and non-TB units and to
#      see how matching changes the distribution.
#
# Output
# - Two figures in the active graphics device:
#    (a) Love plot of covariate balance with a 0.1 threshold line
#    (b) Histogram diagnostic from the MatchIt object
# ------------------------------------------------------------------------------
library(cobalt)
m.out=readRDS('../3_secondstage/match.rds')

love.plot(m.out, abs = TRUE, thresholds = c(m = .1))

plot(m.out, type = "hist")
