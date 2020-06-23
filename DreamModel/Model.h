#ifndef Model_H_
#define Model_H_
#include <vector>
#include <iostream>
#include "TMath.h"
#include "TObject.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "WaveFunction.h"

class Model {
 public:
  Model();

  virtual ~Model();

  void SetrRange(float rmin, float rmax) {
    frmin = rmin;
    frmax = rmax;
  }
  ;
  void SetDrawWF(bool draw) {
    fDrawWF = draw;
  }
  ;
  void SetGaussWF(bool wf) {
    fWFGauss = wf;
  }
  ;
  void SetTwoGaussWF(bool wf) {
    fWFTwoGauss = wf;
  }
  ;
  void SetHulthenWF(bool wf) {
    fWFHulthen = wf;
  }
  ;
  void SetChiEFTWF(bool wf) {
    fWFchiEFT = wf;
  }
  ;

  double GetWaveFunction(const unsigned &WhichMomBin);
  double GetB2_Classical(const unsigned &WhichMomBin);
  double GetB2_kfir(const unsigned &WhichMomBin);
  double GetB2_HulthenKfir(const unsigned &WhichMomBin);
  double GetB2Hulthen(const unsigned &WhichMomBin);
  double GetB2TwoGaussian(const unsigned &WhichMomBin);
  double GetB2ChiEFT(const unsigned &WhichMomBin);
  double GetMomentum(const unsigned &WhichMomBin) const;
  double GetMomBinLowEdge(const unsigned &WhichMomBin) const;
  double GetMomBinUpEdge(const unsigned &WhichMomBin) const;
  void DelAllMom();
  void DelMom();
  void SetMomBins(const unsigned &nummombins, const double *mombins,
                  const double *bincenter = NULL);
  void SetMomBins(const unsigned &nummombins, const double &MinMom,
                  const double &MaxMom);
  static double dB2Hulthen(double *x, double *par);
  static double dB2TwoGaussian(double *x, double *par);
  double GetEvaluate_d_u(double r);
  double GetEvaluate_d_w(double r);
  void LoadWF();
  static double dB2ChiEFT(double *x, double *par);


  enum NotificationOptions {
    nSilent,
    nError,
    nWarning,
    nAll
  };

 private:
  WaveFunction Wavefun;
  double fWaveF;
  double fr;
  float frmin;
  float frmax;
  unsigned NumMomBins;
  double *MomBin;
  double *MomBinCenter;
  double *rWF;
  double *kWFErr;
  short Notifications;
  bool fDrawWF;
  bool fWFGauss;
  bool fWFTwoGauss;
  bool fWFHulthen;
  bool fWFchiEFT;

};
#endif

