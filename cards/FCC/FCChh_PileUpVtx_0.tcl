#
# Official Delphes card prepared by FCC-hh collaboration
#
#  Main authors:  Michele Selvaggi (CERN)
#
#  Released on: Dec. 1st, 2017
#
#  Configuration: FCC-hh baseline detector
#
#######################################
# Order of execution of various modules
#######################################

set MaxEvents 10000

set ExecutionPath {

  PileUpMerger
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  EnergyLoss

  TrackMerger
  TrackSmearing

  ECal
  HCal

  TimeSmearingMIP
  TimeSmearingPhotons
  TimeSmearingNH

  TimeSmearingEcal
  TimeSmearingHcal

  VertexFinderDA4D
  PileUpSubtractor4D

  ElectronFilter
  MuonFilter
  ChargedHadronFilter

  Calorimeter
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  ElectronIsolation

  MuonIsolation

  NeutrinoFilter

  MissingET
  GenMissingET

  GenJetFinder04
  FastJetFinder04
  CaloJetFinder04
  TrackJetFinder04

  JetEnergyScale

  JetFlavorAssociation

  BTagging
  CTagging
  TauTagging

  ScalarHT

  UniqueObjectFinder

  TreeWriter
}

###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  set InputArray Delphes/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  # pre-generated minbias input file
  set PileUpFile MinBias.pileup

  # average expected pile up
  set MeanPileUp 0

  # 0-poisson, 1-uniform, 2-delta
  set PileUpDistribution 2

  # maximum spread in the beam direction in m
  set ZVertexSpread 0.25

  # maximum spread in time in s
  set TVertexSpread 800E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)

  #set VertexDistributionFormula {exp(-(t^2/(2*(0.063/2.99792458E8*exp(-(z^2/(2*(0.063)^2))))^2)))}
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}

  # taking 5.3 cm x 160 ps

  #set VertexDistributionFormula { (abs(t) <= 160e-12) * (abs(z) <= 0.053) * (1.00) +
  #                                (abs(t) >  160e-12) * (abs(z) <= 0.053) * (0.00) +
  #                               (abs(t) <= 160e-12) * (abs(z) > 0.053)  * (0.00) +
  #                           (abs(t) >  160e-12) * (abs(z) > 0.053)  * (0.00)}

}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles
  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.5
  # half-length of the magnetic field coverage, in m
  set HalfLength 5

  # magnetic field
  set Bz 4.0
}


####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

 # TBC (which eta_max ? which pT min?)

 # tracking efficiency formula for charged hadrons

  set EfficiencyFormula { (pt <= 0.5) * (0.00) + \
(abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.90) + \
(abs(eta) <= 2.5) * (pt > 1) * (0.95) + \
(abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.85) + \
(abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (0.90) + \
(abs(eta) > 4 && abs(eta) <= 6) * (pt > 0.5 && pt <= 1) * (0.80) + \
(abs(eta) > 4 && abs(eta) <= 6) * (pt > 1.0) * (0.85) + \
(abs(eta) > 6.0) * (0.00)}

}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

# TBC (which eta_max ?)
# putting same as charged hadrons for now...

  set EfficiencyFormula { (pt <= 0.5) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.90) + \
  (abs(eta) <= 2.5) * (pt > 1) * (0.95) + \
  (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.85) + \
  (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (0.90) + \
  (abs(eta) > 4 && abs(eta) <= 6) * (pt > 0.5 && pt <= 1) * (0.80) + \
  (abs(eta) > 4 && abs(eta) <= 6) * (pt > 1.0) * (0.85) + \
  (abs(eta) > 6.0) * (0.00)}

}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

# TBC (which eta_max ? why eff = 0 for 4 < eta < 6 ? for now put the same as central)
# what about high pT ?
 # tracking efficiency formula for muons
  set EfficiencyFormula { (pt <= 0.5) * (0.00) + \
  (abs(eta) <= 6.0) * (pt > 0.5 && pt <= 1) * (0.90) + \
  (abs(eta) <= 6.0) * (pt > 1) * (0.99) + \
  (abs(eta) > 6.0) * (0.00)}

}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  source momentumResolutionVsP.tcl
}


###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  source momentumResolutionVsP.tcl
}


###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # TBC for just putting tracker resolution/ need to add improvement at high pT

  source muonMomentumResolutionVsP.tcl
}

#################################
# Charged Energy Loss (dE/dx)
#################################

module EnergyLoss EnergyLoss {
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons


  # absolute resolution per measurement (normalized in MeV/cm)
  # CMS pixel detector performance is reproduceable with r = 0.4
  # dedicated dEdX detector can achieve r = 0.0 or below (i.e better than Landau)
  
  #set Resolution 0.4
  set Resolution 0.2

  # fraction of measurements to ignore when computing truncated mean
  # suggested range [0.4-0.6]
  
  set TruncatedMeanFraction 0.5

  # detector properties (active fraction = nhits*thickness/L)
  set Thickness 100E-6
  set ActiveFraction 0.0006666

  # Silicon properties, for other materials:
  # cf. http://pdg.lbl.gov/2014/AtomicNuclearProperties/properties8.dat

  set Z 14.
  set A 28.0855
  set rho 2.329

  # material polarisation correction parameters
  set a 0.1492
  set m 3.2546
  set x0 0.2015
  set x1 2.8716
  set I 173.0
  set c0 4.4355

}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}



########################################
#   Smear tracks
########################################

module TrackSmearing TrackSmearing {
  set InputArray TrackMerger/tracks
  set OutputArray tracks
  set ApplyToPileUp true

  # this is CMS resolution is d0 , dZ on top of the smearing in P applied before
  source trackResolutionFCChh.tcl
  # FIXME !!!! we need to add track resolution of FCC-hh baseline detector !!!!!
}


#############
#   ECAL
#############

# TBC : calos seems ok, check eta max value though.

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackSmearing/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true

  set EnergyMin 0.5
  set EnergySignificanceMin 2.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.012 rad towers up to eta = 2.5 (barrel)
   set PhiBins {}
  for {set i -256} {$i <= 256} {incr i} {
    add PhiBins [expr {$i * $pi/256.0}]
  }

  # 0.01 unit in eta up to eta = 2.5 (barrel)
  for {set i -249} {$i <= 250} {incr i} {
    set eta [expr {$i * 0.01}]
    add EtaPhiBins $eta $PhiBins
  }

  # 0.025 rad between 2.5 and 6.0
  set PhiBins {}
  for {set i -128} {$i <= 128} {incr i} {
    add PhiBins [expr {$i * $pi/128.0}]
  }

  # 0.025 unit in eta between 2.5 and 6.0
  for {set i 0} {$i <= 140} {incr i} {
    set eta [expr { -6.00 + $i * 0.025}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 0} {$i <= 139} {incr i} {
    set eta [expr { 2.525 + $i * 0.025}]
    add EtaPhiBins $eta $PhiBins
  }

  # default energy fractions {abs(PDG code)} {fraction of energy deposited in ECAL}

  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  # add EnergyFraction {310} {0.3}
  # add EnergyFraction {3122} {0.3}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {                     (abs(eta) <= 4.0) * sqrt(energy^2*0.007^2 + energy*0.10^2) + \
                            (abs(eta) > 4.0 && abs(eta) <= 6.0) * sqrt(energy^2*0.035^2  + energy*0.30^2)}


}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false

  set EnergyMin 1.0
  set EnergySignificanceMin 2.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.025 rad towers up to eta = 2.5 (barrel)
   set PhiBins {}
  for {set i -128} {$i <= 128} {incr i} {
    add PhiBins [expr {$i * $pi/128.0}]
  }

  # 0.025 unit in eta up to eta = 2.5 (barrel)
  for {set i -99} {$i <= 100} {incr i} {
    set eta [expr {$i * 0.025}]
    add EtaPhiBins $eta $PhiBins
  }

  # 0.05 rad between 2.5 and 6.0
  set PhiBins {}
  for {set i -64} {$i <= 64} {incr i} {
    add PhiBins [expr {$i * $pi/64.0}]
  }

  # 0.05 unit in eta between 2.5 and 6.0
  for {set i 0} {$i <= 70} {incr i} {
    set eta [expr { -6.00 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 0} {$i <= 69} {incr i} {
    set eta [expr { 2.55 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }


  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  # add EnergyFraction {310} {0.7}
  # add EnergyFraction {3122} {0.7}

   # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {                     (abs(eta) <= 1.7) * sqrt(energy^2*0.03^2 + energy*0.50^2) + \
                            (abs(eta) > 1.7 && abs(eta) <= 4.0) * sqrt(energy^2*0.03^2 + energy*0.60^2) + \
                            (abs(eta) > 4.0 && abs(eta) <= 6.0) * sqrt(energy^2*0.10^2 + energy*1.00^2)}
}


########################################
#   Time Smearing Neutral MIP
########################################

module TimeSmearing TimeSmearingMIP {
  set InputArray HCal/eflowTracks
  set OutputArray tracks

  # assume constant 30 ps resolution for now
  set TimeResolution {
                       (abs(eta) > 0.0 && abs(eta) <= 6.0)* 30E-12
                     }
}

########################################
#   Time Smearing Neutral Photons
########################################

module TimeSmearing TimeSmearingPhotons {
  set InputArray ECal/eflowPhotons
  set OutputArray photons
  set TimeResolution {
                        (abs(eta) > 0.0 && abs(eta) <= 6.0) * sqrt(20e-12^2 + 150e-12^2)/energy^2
                     }
}

########################################
#   Time Smearing Neutral NeutralHadrons
########################################
#
module TimeSmearing TimeSmearingNH {
  set InputArray HCal/eflowNeutralHadrons
  set OutputArray neutralhadrons

  # assume resolution for now
  set TimeResolution {
                        (abs(eta) > 0.0 && abs(eta) <= 6.0) * sqrt(20e-12^2 + 150e-12^2)/energy^2
                     }
}

########################################
#   Time Smearing Neutral Ecal
########################################

module TimeSmearing TimeSmearingEcal {
  set InputArray ECal/ecalTowers
  set OutputArray towers
  set TimeResolution {
                        (abs(eta) > 0.0 && abs(eta) <= 6.0) * sqrt(20e-12^2 + 150e-12^2)/energy^2
                     }
}


########################################
#   Time Smearing Neutral Hcal
########################################

module TimeSmearing TimeSmearingHcal {
  set InputArray HCal/hcalTowers
  set OutputArray towers
  set TimeResolution {
                        (abs(eta) > 0.0 && abs(eta) <= 6.0) * sqrt(20e-12^2 + 150e-12^2)/energy^2
                     }
}




##################################
# Primary vertex reconstruction
##################################

module VertexFinderDA4D VertexFinderDA4D {
  set InputArray TimeSmearingMIP/tracks

  set OutputArray tracks
  set VertexOutputArray vertices

  set Verbose 0
  set MinPT 1.0

  # in mm
  set VertexSpaceSize 0.5

  # in s
  set VertexTimeSize 10E-12

  set UseTc 1
  set BetaMax 0.1
  set BetaStop 1.0
  set CoolingFactor 0.8
  set MaxIterations 100

  # in mm
  set DzCutOff 40
  set D0CutOff 30

}

##########################
# PileUpSubtractor4D
##########################

module  PileUpSubtractor4D PileUpSubtractor4D {
# add InputArray InputArray OutputArray

  add InputArray VertexFinderDA4D/tracks tracks
  add InputArray TimeSmearingPhotons/photons photons
  add InputArray TimeSmearingNH/neutralhadrons neutralhadrons
  add InputArray TimeSmearingEcal/towers towers
  add InputArray TimeSmearingHcal/towers towers

  set VertexInputArray VertexFinderDA4D/vertices

  set fChargedMinSignificance 10000.0
  set fNeutralMinSignificance 10000.0
}


#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray PileUpSubtractor4D/tracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray PileUpSubtractor4D/tracks
  set OutputArray muons
  set Invert true
  add PdgCode {13}
  add PdgCode {-13}
}


######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  set InputArray PileUpSubtractor4D/tracks
  set OutputArray chargedHadrons

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
# add InputArray InputArray
  add InputArray PileUpSubtractor4D/towers
  add InputArray MuonFilter/muons
  set OutputArray towers
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray PileUpSubtractor4D/tracks
  add InputArray PileUpSubtractor4D/photons
  add InputArray PileUpSubtractor4D/neutralhadrons
  set OutputArray eflow
}

##############
# EFlowFilter
##############

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray eflow

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}



###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set EnergyOutputArray energy
}

#################
# Neutrino Filter
#################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}


#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {

# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}


#####################
# MC truth jet finder
#####################

# TBC: is jet radius fine?

module FastJetFinder GenJetFinder04 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}



##################
# Fast Jet finder
##################

module FastJetFinder FastJetFinder04 {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder CaloJetFinder04 {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder TrackJetFinder04 {
  set InputArray TrackMerger/tracks

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}


##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder04/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}


########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 5.0
  set PartonEtaMax 6.0

}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray PileUpSubtractor4D/photons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  set EfficiencyFormula {
  (pt <= 1.0) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 1.0 && pt < 5.0)  * (0.70) +
  (abs(eta) <= 2.5) * (pt > 5.0 && pt < 10.0) * (0.85) +
  (abs(eta) <= 2.5) * (pt > 10.0)             * (0.95) +

  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt < 5.0)  * (0.60) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 5.0 && pt < 10.0) * (0.80) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0)   * (0.90) +

  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1.0 && pt < 5.0)  * (0.50) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 5.0 && pt < 10.0) * (0.70) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 10.0)             * (0.80) + \
  (abs(eta) > 6.0) * (0.00)}

}

##################
# Photon isolation
##################

# TBC: check values for iso cuts

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}


####################
# Electron isolation
####################

# TBC: check values for iso cuts

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronFilter/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}


################
# Muon isolation
################

# TBC: check values for iso cuts

module Isolation MuonIsolation {
  set CandidateInputArray MuonMomentumSmearing/muons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray muons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.2
}



###########
# b-tagging
###########

module BTagging BTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  add EfficiencyFormula {0} {

  (pt <= 10.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 500)      * (0.01) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 15000.0) * (0.01)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5)                   * (pt > 15000.0)               * (0.00) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 500)      * (0.0075) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 15000.0) * (0.0075)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 15000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 10.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 500)      * (0.05) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 15000.0) * (0.05)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5)                   * (pt > 15000.0)               * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 500)      * (0.03) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 15000.0) * (0.03)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 15000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                    * (pt > 10.0 && pt < 500)      * (0.85) +
  (abs(eta) < 2.5)                    * (pt > 500.0 && pt < 15000.0) * (0.85)*(1.0 - pt/15000.) +
  (abs(eta) < 2.5)                    * (pt > 15000.0)               * (0.000) +
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 500)      * (0.64) +
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 15000.0) * (0.64)*(1.0 - pt/15000.) +
  (abs(eta) <= 2.5 && abs(eta) < 4.0) * (pt > 15000.0)               * (0.000) +
  (abs(eta) >= 4.0) * (0.00)}

}

###########
# c-tagging
###########

module BTagging CTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 1

  add EfficiencyFormula {0} {

  (pt <= 10.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 10.0) * (0.01) + \
  (abs(eta) > 4.0) * (pt > 10.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 10.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 10.0) * (0.25) + \
  (abs(eta) > 4.0) * (pt > 10.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 10.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 10.0) * (0.03) + \
  (abs(eta) > 4.0) * (pt > 10.0) * (0.00)}

}


#############
# tau-tagging
#############


module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  add EfficiencyFormula {0} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 5000.0)    * (0.01) + \
  (abs(eta) < 2.5)                   * (pt > 5000.0 && pt < 34000.0) * (0.01)  *(8./9. - pt/30000.) + \
  (abs(eta) < 2.5)                   * (pt > 34000.0)                * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 5000.0)    * (0.0075) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 5000.0 && pt < 34000.0) * (0.0075)*(8./9. - pt/30000.) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 34000.0)                * (0.00) + \
  (abs(eta) > 4.0)                   * (0.00)}

  add EfficiencyFormula {11} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 5000.0)    * (0.005) + \
  (abs(eta) < 2.5)                   * (pt > 5000.0 && pt < 34000.0) * (0.005)  *(8./9. - pt/30000.) + \
  (abs(eta) < 2.5)                   * (pt > 34000.0)                * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 5000.0)    * (0.00375) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 5000.0 && pt < 34000.0) * (0.00375)*(8./9. - pt/30000.) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 34000.0)                * (0.00) + \
  (abs(eta) > 4.0)                   * (0.00)}

  add EfficiencyFormula {15} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5) * (pt > 10.0 && pt < 5000.0)                      * (0.8)                      + \
  (abs(eta) < 2.5) * (pt > 5000.0 && pt < 34000.0)                   * (0.8) *(8./9. - pt/30000.) + \
  (abs(eta) < 2.5)                   * (pt > 34000.0)                * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 5000.0)    * (0.65) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 5000.0 && pt < 34000.0) * (0.65)*(8./9. - pt/30000.) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 34000.0)                * (0.00) + \
  (abs(eta) > 4.0)                                                   * (0.00)}

}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale/jets jets

}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch PileUpMerger/stableParticles Particle GenParticle
  add Branch PileUpMerger/vertices VertexGen Vertex
  add Branch VertexFinderDA4D/vertices Vertex4D Vertex
  add Branch VertexFinderDA4D/tracks Tracks4D Track
  add Branch TimeSmearingMIP/tracks  MIPTracks Track

  add Branch UniqueObjectFinder/photons Photon Photon
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/muons Muon Muon
  add Branch UniqueObjectFinder/jets Jet Jet

  add Branch GenJetFinder04/jets GenJet04 Jet
  add Branch FastJetFinder04/jets Jet04 Jet
}
