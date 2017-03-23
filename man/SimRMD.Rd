 
\name{SimRMD}
\alias{SimRMD}
\title{Simulate trial characters of the phase I dose-finding design.}
\description{
  \code{SimRMD} use the design parameters to simulate the trial characters of phase I dose-finding design.
}
\usage{
SimRMD(seed=2441, strDose=1, chSize=3, trlSize=36,
   numTrials=1000, sdose=1:6, MaxCycle=6, tox.target=0.28,
   control, iter=10000, burnin=4000, thin=1, chains=1,    
   pathout='./', sc, trend=0.1)
}

\arguments{
  \item{seed}{Random seed of the simulation (default = 2014).}
  \item{strDose}{Start dose (default = 1).}
  \item{chSize}{Size of each patient cohort (default = 3).}
  \item{trlSize}{Total number of patients in the trial (default = 36).}
  \item{numTrials}{Number of simulated trials (default = 1000).}
  \item{sdose}{Dose range (default = 1:6).}
  \item{MaxCycle}{Maximum treatment cycle (default = 6).}
  \item{tox.target}{Target toxicity score (default = 0.28).}
  \item{control}{Prior distributions as specified by \code{parm} function.}
  \item{iter}{Total number of MCMC simulations (default = 4000).}
  \item{burnin}{Number of burn-ins in the MCMC simulation (default = 1000).}
  \item{thin}{Thinning parameter (default = 1).}
  \item{chains}{Number of chains in the MCMC simulation (default = 1).} 
  \item{pathout}{The output directory of simulation results (default = './').}
  \item{sc}{A list containing the Dose, Cycle toxicity score matrix of each toxicity type: 'probaT1', 'probaT2', 'probaT3'.}
  \item{trend}{Time trend of toxicity scores over multiple treatment cycles (default = 0.1).}
 
}

\value{
  A list containing the following elements:
  \item{DOSE-RECOMMENDED}{Recommended dose for the next patient cohort. }
  \item{Estimate}{Mean, SD, and median estimates of the toxicity score at each dose level.}
  \item{Quantiles}{2.5\%, 25\%, 50\%, 75\%, 97.5\% quantile estimates of the toxicity score at each dose level.}
}

\details{
The \code{SimRMD} function uses the patient data structure \code{patdata} to fit the linear mixed model as specified by \code{formula} based on the prior distributions as specified by \code{control}. The parameters of the MCMC simulation are specified by \code{iter}: the number of iterations, \code{burnin}: burn-ins, \code{thin}: thinning parameter and \code{chains}: number of chains in the MCMC simulation. The target nTTP score is defined by \code{tox.target}, default = 0.28.
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
data(probaT)

# Estimate the dose-recommendation of the next patient cohort based on the toxicity data:
SimRMD(seed=2441, strDose=1, chSize=3, trlSize=9,
  numTrials=1, sdose=1:6, MaxCycle=5, tox.target=0.28,
  control=control, iter=500, burnin=100, thin=1, chains=1,    
  pathout='./', sc=probaT, trend=0.1)
 
}
