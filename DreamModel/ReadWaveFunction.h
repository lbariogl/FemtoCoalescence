#ifndef ReadWaveFunction_H_
#define ReadWaveFunction_H_
#include <vector>
#include "DLM_Histo.h"
#include "DLM_CppTools.h"
#include "TMath.h"
#endif
void ReadDeuteronWF(const char *InputFileName,
                                  DLM_Histo<float> &OutputU,
                                  DLM_Histo<float> &OutputW) {
  FILE *InFile;
  InFile = fopen(InputFileName, "r");
  if (!InFile) {
    printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n",
           InputFileName);
    return;
  }

  char *cdummy = new char[1024];
  unsigned CurrentLine = 0;
  float fRadius, fuWF, fwWF;
  double *Radius = new double[1024];
  double *RadiusBins = new double[1025];

  unsigned NumBins = 0;
  RadiusBins[0] = 0;
  //read line-by-line until the end of the file
  while (!feof(InFile)) {
    CurrentLine++;
    //read the next line
    if (!fgets(cdummy, 1023, InFile))
      continue;
    if (CurrentLine <= 3)
      continue;
    //printf("Line#%u: %s\n",CurrentLine,cdummy);
    sscanf(cdummy, "%f %f %f", &fRadius, &fuWF, &fwWF);
    //printf("Line#%u: %.2e; %.2e; %.2e;\n",CurrentLine,fRadius,fuWF,fwWF);

    //if((CurrentLine)%int(NumRadBins)==0){
    //   sscanf(cdummy, "%f",&fMomentum);
    Radius[NumBins] = fRadius;
    if (NumBins) {
      //set the bin range in between the last two bin centers
      RadiusBins[NumBins] = 0.5 * (Radius[NumBins] + Radius[NumBins - 1]);
    }
    //}

    NumBins++;
  }
  fclose(InFile);

  //set the upper edge of the last bin, where we just add the bin width of the last bin
  //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
  RadiusBins[NumBins] = 2. * Radius[NumBins - 1] - RadiusBins[NumBins - 1];

  //for(unsigned uBin=0; uBin<NumBins; uBin++){
  //    printf("#%u: %f %f%f\n",uBin,RadiusBins[uBin],Radius[uBin],RadiusBins[uBin+1]);
  //}

  //dimension (1D)
  OutputU.SetUp(1);
  OutputU.SetUp(0, NumBins, RadiusBins, Radius);
  OutputU.Initialize();

  OutputW.SetUp(1);
  OutputW.SetUp(0, NumBins, RadiusBins, Radius);
  OutputW.Initialize();

  CurrentLine = 0;
  NumBins = 0;
  InFile = fopen(InputFileName, "r");
  while (!feof(InFile)) {
    CurrentLine++;
    //read the next line
    if (!fgets(cdummy, 1023, InFile))
      continue;
    if (CurrentLine <= 3)
      continue;
    sscanf(cdummy, "%f %f %f", &fRadius, &fuWF, &fwWF);
    OutputU.SetBinContent(NumBins, fuWF);
    OutputW.SetBinContent(NumBins, fwWF);
    NumBins++;
  }
  fclose(InFile);

  delete[] cdummy;
  delete[] Radius;
}

DLM_Histo<float> huWF;
DLM_Histo<float> hwWF;

static float Evaluate_d_u(double Radius) {
  return huWF.Eval(&Radius);
}
static float Evaluate_d_w(double Radius) {
  return hwWF.Eval(&Radius);
}

void TestDLMhisto() {
  ReadDeuteronWF("/home/sbhawani/Coalescence/WaveFunctions/WFChiEFT/deuwaves_n4lo500_rspace.d", huWF, hwWF);
}

