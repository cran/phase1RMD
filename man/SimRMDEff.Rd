 
\name{SimRMDEff}
\alias{SimRMDEff}
\title{Simulate trial characters of the phase I dose-finding design using toxicity and efficacy data.}
\description{
  \code{SimRMDEff} use the design parameters to simulate the trial characters of phase I dose-finding design using both toxicity and efficacy data.
}
\usage{
SimRMDEff(seed=2441, tox_matrix, 
   strDose=1, chSize=3, trlSize=36,
   sdose=1:6, MaxCycle=6, tox.target=0.28, 
   eff.structure=c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9))
}

\arguments{
  \item{seed}{Random seed of the simulation (default = 2441).}
  \item{tox_matrix}{Four-dimension array contains Dose, Cycle, Toxicity Type, Toxicity Grade.}
  \item{strDose}{Start dose (default = 1).}
  \item{chSize}{Size of each patient cohort (default = 3).}
  \item{trlSize}{Total number of patients in the trial (default=36).}
  \item{sdose}{Dose range (default = 1:6).}
  \item{MaxCycle}{Maximum treatment cycle (default = 6).}
  \item{tox.target}{Target toxicity score (default = 0.28).}  
  \item{eff.structure}{The efficacy outcome for each dose (default = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9)).}
}
\value{
  A list containing the following elements:
  \item{sc}{Mean nTTP of each cycle, mean efficacy and mean DLT at all dose levels.}
  \item{opt.dose}{Dose selection.}
}
\details{
The \code{SimRMDEff} function simulate the operation characteristics.

}
\author{
  Jun Yin
}
\examples{
data(tox_matrix)

# Simulate the operating characteristics of clinical trials:
SimRMDEff(seed = 2441,
 tox_matrix = tox_matrix,
 strDose = 1, chSize = 3, trlSize = 12, 
 sdose = 1:6, MaxCycle = 6, tox.target = 0.28,
 eff.structure = c(0.1, 0.2, 0.3, 0.4, 0.7, 0.9))
 
 
}
