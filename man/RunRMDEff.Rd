 
\name{RunRMDEff}
\alias{RunRMDEff}
\title{Toxicity and Efficacy-based Dose estimation for the next patient cohort.}
\description{
  \code{RunRMDEff} estimate the dosage for the next patient cohort using TTP measured multidimensional toxicity data collected from multiple treatemnt cycles of the existing patients.
}
\usage{
  RunRMDEff(toxicity.dat, efficacy.dat=NULL, tox.target=0.28, sdose=1:6, MaxCycle=6)
}
\arguments{
  \item{toxicity.dat}{Toxicty data collected from the multiple treatment cycles of all patient cohorts. }
  \item{efficacy.dat}{Efficacy data collected from all patient cohorts (default=NULL).}
  \item{tox.target}{Target toxicity score (default = 0.28).}
  \item{sdose}{Dose range (default = 1:6).}
  \item{MaxCycle}{Maximum treatment cycle (default = 6).}
}
\value{
  A list containing the following elements:
  \item{DOSE-RECOMMENDED}{Recommended dose for the next patient cohort. }
  \item{ToxEstimate}{Mean estimates of the toxicity score at each dose level.}
  \item{EffEstimate}{Mean estimates of the efficacy score at each dose level.}
}
\details{
The \code{RunRMDEff} function uses the patient data structure \code{patdata} to fit the linear mixed model using both toxicity and efficacy scores.

}
\author{
  Jun Yin
}
\examples{
data(eff_dat);
data(tox_dat);

# Estimate the dose-recommendation of the next patient cohort based on the toxicity data:
RunRMDEff(toxicity.dat=tox_dat, efficacy.dat=eff_dat, tox.target=0.28, sdose=1:6, MaxCycle=2)
 
}
