#include "WaveFunction.h"
#include "ReadWaveFunction.h"
WaveFunction::WaveFunction() {

}
WaveFunction::~WaveFunction() {

}
double WaveFunction::WFGauss(double r) {
  float d = 3.2;
  return pow(TMath::Pi() * d * d, -0.75) * (TMath::Exp(-r * r / 2. / d / d));
}
double WaveFunction::WFTwoGauss(double r) {
  float del = 0.581;
  float d1 = 3.979;
  float d2 = 0.890;
  return pow(TMath::Pi(), -0.75)
      * (pow(del, 0.5) / (pow(d1, 1.5)) * TMath::Exp(-r * r / 2. / d1 / d1 / 2)
          + pow(1 - del, 0.5) / (pow(d1, 1.5))
              * TMath::Exp(-r * r / 2. / d1 / d1 / 2));
}
double WaveFunction::WFHulthen(double r) {
  float a = 0.23;
  float b = 1.61;
  return pow(a * b * (a + b) / (2 * TMath::Pi() * pow((a - b), 2.)), 0.5)
      * (TMath::Exp(-a * r) - TMath::Exp(-b * r)) / r;

}

//it important to read the wave funtion to use in computation
double WaveFunction::WaveChiEFT(double r) {
  double u = Evaluate_d_u(r);
  double w = Evaluate_d_w(r);
  return 1 / (pow(4 * TMath::Pi() * r * r, 0.5)) * u;

}
double WaveFunction::uWaveChiEFT(double r) {
  double u = Evaluate_d_u(r);
  return u;
}
double WaveFunction::wWaveChiEFT(double r) {
  double w = Evaluate_d_w(r);
  return w;

}
void WaveFunction::LoadTheWF() {
  TestDLMhisto();
}
/*double WaveFunction::DensityChiEFT(double r) {
  ReadDeuteronWF("/home/sbhawani/Downloads/Coalescence_Task/Script/TwoGaussian/WFS/deuwaves_n4lo500_rspace.d", huWF, hwWF);
  double u = Evaluate_d_u(r);
  double w = Evaluate_d_w(r);
  double density =  u * u+w * w;
  return  density;

}*/
