 
\name{RunRMD}
\alias{RunRMD}
\title{Toxicity-based Dose estimation for the next patient cohort.}
\description{
  \code{RunRMD} estimate the dosage for the next patient cohort using TTP measured multidimensional toxicity data collected from multiple treatemnt cycles of the existing patients.
}
\usage{
  RunRMD(data, control, iter=6000, burnin=1000, thin=1, chains=1, tox.target=0.28)
}
\arguments{
  \item{data}{Toxicty data collected from the multiple treatment cycles of existing patient cohorts. }
  \item{control}{Prior distributions as specified by \code{parm} function.}
  \item{iter}{Total number of MCMC simulations (default = 4000).}
  \item{burnin}{Number of burn-ins in the MCMC simulation (default = 1000).}
  \item{thin}{Thinning parameter (default = 1).}
  \item{chains}{Number of chains in the MCMC simulation (default = 1).}
  \item{tox.target}{Target toxicity score (default = 0.28).}
  
}
\value{
  A list containing the following elements:
  \item{DOSE-RECOMMENDED}{Recommended dose for the next patient cohort. }
  \item{Estimate}{Mean, SD, and median estimates of the toxicity score at each dose level.}
  \item{Quantiles}{2.5\%, 25\%, 50\%, 75\%, 97.5\% quantile estimates of the toxicity score at each dose level.}
}
\details{
The \code{RunRMD} function uses the patient data structure \code{patdata} to fit the linear mixed model as specified by \code{formula} based on the prior distributions as specified by \code{control}. The parameters of the MCMC simulation are specified by \code{iter}: the number of iterations, \code{burnin}: burn-ins, \code{thin}: thinning parameter and \code{chains}: number of chains in the MCMC simulation. The target nTTP score is defined by \code{tox.target}, default = 0.28.

}
\author{
  Jun Yin
}
\examples{
# Setup the prior distributions for the phase I dose-finding model:
control <- list(
 beta.dose = parm("normal", mean = 0, var = 1000),
 beta.other = parm("normal", mean = 0, var = 1000 ),
 gamma = parm("normal", mean = 0, var = 100 ),
 s2.gamma = parm("invgamma", shape = 0.001, scale = 0.001),
 s2.epsilon = parm("invgamma", shape = 0.001, scale = 0.001)
 )
 
# Load the example patient toxicity data:
data(patdata)

# Estimate the dose-recommendation of the next patient cohort based on the toxicity data:
RunRMD(data=patdata, control=control, tox.target=0.28, 
  iter=2000, burnin=500, thin=1, chains=1)
 
}
