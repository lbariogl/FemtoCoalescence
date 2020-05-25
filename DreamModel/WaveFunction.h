#ifndef WaveFunction_H_
#define WaveFunction_H_
#include <vector>
#include "DLM_Integration.h"
#include "TMath.h"

class Model;
class WaveFunction {
  friend class Model;
 public:
  WaveFunction();
  virtual ~WaveFunction();
  void LoadTheWF();
  double WFGauss(double r);
  double WFTwoGauss(double r);
  double WFHulthen(double r);
  double WaveChiEFT(double r);
  static double uWaveChiEFT(double r);
  static double wWaveChiEFT(double r);
  static double DensityChiEFT(double r);
  double WF(double r);
 private:
};
#endif /* WaveFunction_H_ */
