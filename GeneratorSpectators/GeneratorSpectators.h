#ifndef O2_GENERATORSPECTATORS_H
#define O2_GENERATORSPECTATORS_H
// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Class to generate spectators for ZDC simulations

#include <TF1.h>
#include <TGenerator.h>

class GeneratorSpectators : public TGenerator {

public:
  GeneratorSpectators();

  virtual ~GeneratorSpectators() {};

  virtual void Init();
  virtual void GenerateEvent();
  virtual int ImportParticles(TClonesArray *particles, Option_t *option);

  // Parameters that could be set for generation
  void SetDebug() { fDebug = kTRUE; }
  void SetImpactParameter(Float_t b = -1.) { fImpactParameter = b; }
  void SetNpartFluctuation(Bool_t flag = kTRUE) { fNpartFluctuation = flag; }
  void SetNpart(Int_t npart = -1) { fNpart = npart; }
  void SetParticle(Int_t pdgcode = 2112) { fPDGcode = pdgcode; }
  void SetMomentum(Float_t ptot = 2680.) { fPtot = ptot; }
  void SetDirection(Float_t eta = 0, Float_t cosx = 0, Float_t cosy = 0,
                    Float_t cosz = 1) {
    fPseudoRapidity = eta;
    fCosx = cosx;
    fCosy = cosy;
    fCosz = cosz;
  }
  void SetFermi(Bool_t flag = kTRUE) { fFermiflag = flag; }
  void SetDivergence(Float_t bmdiv = 0.) { fBeamDiv = bmdiv; }
  void SetCrossing(Float_t xingangle = 0., Int_t xingplane = 2) {
    fBeamCrossAngle = xingangle;
    fBeamCrossPlane = xingplane;
  }

  // Getters
  Double_t GetFermi2p(Int_t key) const { return fProbintp[key]; }
  Double_t GetFermi2n(Int_t key) const { return fProbintn[key]; }
  Float_t GetZDirection() const { return fCosz; }

protected:
  Bool_t   fDebug;              // Debugging flag
  Float_t  fImpactParameter;    // Impact parameter, if <=0 sample from realistic distribution
  Bool_t   fNpartFluctuation;   // Enable gaussian fluctuation of sampled number of particles
  Int_t    fNpart;              // Number of particles to be generated
                                // if>0 overwrite the impact-parameter information
  Int_t    fPDGcode;            // Particle to be generated - can be n (2112) or p (2212)
  Float_t  fPtot;               // Nucleon momentum
  Float_t  fPseudoRapidity;     // Pseudorapidity: =0->track director cosines, !=0->pc.eta
  Float_t  fCosx;               // Director cos of the track - x direction
  Float_t  fCosy;               // Director cos of the track - y direction
  Float_t  fCosz;               // Director cos of the track - z direction
  Bool_t   fFermiflag;          // Fermi momentum flag (true -> Fermi smearing)
  Float_t  fBeamDiv;            // Beam divergence (angle in rad)
  Float_t  fBeamCrossAngle;     // Beam crossing angle (angle in rad)
  Int_t    fBeamCrossPlane;     // Beam crossing plane
                                // (=1 -> horizontal, =2 -> vertical plane)

  Double_t fProbintp[201];      // Protons momentum distribution due to Fermi
  Double_t fProbintn[201];      // Neutrons momentum distribution due to Fermi
  Double_t fPp[201];            // Spectator momenta
  TF1      fNeutronsVsImpPar;   // Parameterization of number of neutrons vs impact parameter
  TF1      fWidthNeutVsImpPar;  // Parameterization of fluctuation of neutrons (Gaussian assumed)
  TF1      fImpParDistr;        // Parameterization of impact parameter distribution

private:
  GeneratorSpectators(const GeneratorSpectators &gen);
  GeneratorSpectators & operator=(const GeneratorSpectators &gen);

  // Sampling of impact parameter and number of particles
  Float_t SampleImpPar();
  Int_t SampleNpart(Float_t impactParameter);

  // Fermi smearing, beam divergence and crossing angle
  void FermiTwoGaussian(Float_t A);
  void ExtractFermi(Int_t id, Double_t *ddp);
  void BeamDivergence(Double_t *pLab);
  void BeamCrossing(Double_t *pLab);
  void AddAngle(Double_t theta1, Double_t phi1, Double_t theta2, Double_t phi2,
                Double_t *angle);

  // Functions for parameterization of impact parameter and number of neutrons
  void InitParameterizations();
  static Double_t ImpParFunc(Double_t *x, Double_t *par);

  ClassDef(GeneratorSpectators, 2) // Generator for spectators
};

#endif
