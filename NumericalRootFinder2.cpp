#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

using namespace std;

double tUniv_years = 12.7-0.7;
double tUniv_seconds = tUniv_years*365*24*60*60*1.0e+09;
double H0_seconds = 2.26853e-18;
double tUnivTimesH0 = tUniv_seconds*H0_seconds;
double k = 0.0;

double myfunc(double x) { 
   double xMat = 1.0 + k - x;
   double xRad = x;

   double numeratorLog = xMat + 2.0*sqrt( (-k)*xRad );
   double denominatorLog = 2*(-k)+xMat+2*sqrt(xRad);
   double logCoef = xMat / ( 2 * pow((-k),3.0/2.0) );
   double leftover = (1.0 - sqrt(xRad) ) / (-k);
   return (logCoef * TMath::Log( numeratorLog / denominatorLog ) + leftover) - tUnivTimesH0; 
}

int NumericalRootFinder2()
{
 
   TCanvas *can = new TCanvas("can","can",800,600);

   const int nPoints = 1000;

   double curv[nPoints];
   double xRad[3][nPoints];
   double xMat[3][nPoints];
   double tUniverse[3] = {12.0,12.7,13.4};

   for (int kk=0; kk<3; kk++){
	   tUniv_years = tUniverse[kk];
	   double tUniv_seconds = tUniv_years*365*24*60*60*1.0e+09;
	   double H0_seconds = 2.26853e-18;
	   tUnivTimesH0 = tUniv_seconds*H0_seconds;
	   for (int i=0; i<nPoints; i++){
		   k = -(0.00001+i*0.99999/nPoints);

		   // Create the function and wrap it
		   TF1* f = new TF1("myfunc", "myfunc(x)", 0.0,1.0+k);
		   ROOT::Math::WrappedTF1 wf1(*f);
		 
		   // Create the Integrator
		   ROOT::Math::BrentRootFinder brf;
		 
		   // Set parameters of the method
		   brf.SetFunction( wf1, 0.0, 1.0+k );
		   brf.Solve();
		 
		   cout << "Omega0 = " << 1+k << "; xMat = " << brf.Root() << "; xRad = " << 1+k-brf.Root() << endl;
		   curv[i] = 1+k;
		   xMat[kk][i] = brf.Root();
		   xRad[kk][i] = 1+k-brf.Root();

	   } 
   }

   TGraph* grRad[3];
   TGraph* grMat[3];
   for (int kk=0; kk<3; kk++){
	   grRad[kk] = new TGraph(nPoints, curv, xRad[kk]);
	   grRad[kk]->SetMarkerColor(kRed);
	   grRad[kk]->SetLineColor(kRed);
	   if (kk==0){
		   grRad[kk]->Draw("AP");
		   grRad[kk]->GetXaxis()->SetTitle("1-#Omega_{0}");
		   grRad[kk]->SetTitle("");
	   }
	   else
		   grRad[kk]->Draw("P same");
	   grMat[kk] = new TGraph(nPoints, curv, xMat[kk]);
	   grMat[kk]->Draw("P same");
   }
   TLegend *leg = new TLegend(0.12,0.7,0.32,0.9);
   leg->AddEntry(grRad[0],"#Omega_{r}","l");
   leg->AddEntry(grMat[0],"#Omega_{m}","l");
   leg->Draw("same");
   can->SaveAs("plot.png");

   return 0;
}

// int main(){
//
//         NumericalRootFinder2();
//         return 0;
// }
