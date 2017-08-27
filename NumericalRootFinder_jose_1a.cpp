#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

using namespace std;

double iM = 0.0;

double myfunc(double x) {
   double sigma = x;

   return iM - exp(3.63*1.0e37*sigma)*4.24*1.0e-26*sqrt(sigma);
   // return iM - exp(3.63*1.0e37*sigma)*2.36*1.0e+25*sqrt(sigma);
   // return iM - exp(3.63*1.0e37*sigma)*4.17*1.0e+24*sqrt(sigma);
}

int NumericalRootFinder_jose()
{
 
   TCanvas *can = new TCanvas("can","can",800,600);

   const int nPoints = 1000;
   double minM = 10;
   double maxM = 1000000;

   double minSigma = 1.0e-50;
   double maxSigma = 1.0e-10;

   double xArr[nPoints];
   double yArr[nPoints];

   for (int i=0; i<nPoints; i++){
	   iM = minM+i*(maxM-minM)/nPoints;
	   xArr[i] = iM;

	   // Create the function and wrap it
	   TF1* f = new TF1("myfunc", "myfunc(x)", minM, maxM);
	   ROOT::Math::WrappedTF1 wf1(*f);
	 
	   // Create the Integrator
	   ROOT::Math::BrentRootFinder brf;
	 
	   // Set parameters of the method
	   brf.SetFunction( wf1, minSigma, maxSigma );
	   brf.SetNpx(100000);
	   brf.SetLogScan(true);
	   // brf.Solve(10000,0.1,0.1);
	   brf.Solve(10000);

	   cout << "[DEBUG]\tM: " << iM << "; sigma: " << brf.Root() << endl;
	   // double sigmaIn = 1.0e-40;
	   // cout << "[DEBUG]\tsigmaIn: " << sigmaIn << "; func: " << f->Eval(sigmaIn) << endl;
	 
	   yArr[i] = brf.Root()*1.0e+36;

   } 

   double sigmaIn = 1.0e-40;
   cout << "sigmaIn: " << sigmaIn <<  endl;
   for (int i=0; i<nPoints; i++){
	   iM = minM+i*(maxM-minM)/nPoints;
	   TF1* f = new TF1("myfunc", "myfunc(x)", minM, maxM);
	   // cout << "[DEBUG]\tiM: " << iM << "; func: " << f->Eval(sigmaIn) << endl;

   }

   TGraph* gr;
   gr = new TGraph(nPoints, xArr, yArr);
   gr->SetMarkerColor(kRed);
   gr->SetLineColor(kRed);
   gr->Draw("AP");

   // TLegend *leg = new TLegend(0.12,0.7,0.32,0.9);
   // leg->AddEntry(grMat[0],"#Omega_{m}","l");
   // leg->Draw("same");
   can->SaveAs("limit.png");

   return 0;
}

// int main(){
//
//         NumericalRootFinder2();
//         return 0;
// }
