#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

using namespace std;

double k = -0.1;
double tUnivTimesH0 = 0.908562;

// error on Universe age:
// 0.500782

// double myfunc(double x ) {
//    return k*x*x-1;
// }
 
double myfunc(double x) { 
   // double xMat = 1.0 + k - x;
   // double xRad = x;
   double xMat = x;
   double xRad = 1.0 + k - x;
   double numeratorLog = xMat + 2.0*sqrt( (-k)*xRad );
   double denominatorLog = 2*(-k)+xMat+2*sqrt(xRad);
   double logCoef = xMat / ( 2 * pow((-k),3.0/2.0) );
   double leftover = (1.0 - sqrt(xRad) ) / (-k);
   return (logCoef * TMath::Log( numeratorLog / denominatorLog ) + leftover) - tUnivTimesH0; 
}

// double myfunc_deriv ( double x ) {
//    return 2*k*x;
// }

int NumericalRootFinder2()
{
 
   const int nPoints = 1;
   TCanvas *can = new TCanvas("can","can",800,600);
   TLegend *leg = new TLegend(0.3,0.65,0.9,0.9);
   TF1 *f[nPoints];
   for (int i=0; i<nPoints; i++){
	   k = 0.01+i*0.98/nPoints;

	   // Create the function and wrap it
	   f[i] = new TF1("myfunc", "myfunc(x)", 0.0,1.0+k);
	   ROOT::Math::WrappedTF1 wf1(*f[i]);
	 
	   // Create the Integrator
	   ROOT::Math::BrentRootFinder brf;
	 
	   // Set parameters of the method
	   brf.SetFunction( wf1, 0.0, 1.0+k );
	   brf.Solve();
	 
	   // cout << "k = " << k << "; root = " << brf.Root() << "; " ;
	   // cout << "sqrt(1/k) = " << TMath::Sqrt(1/k) << endl;
	    
	   // cout << "Omega0 = " << 1+k << "; fract = " << brf.Root()/(1+k) << "; root = " << brf.Root() << endl;
	   cout << "Omega0 = " << 1+k << "; xMat = " << brf.Root() << "; xRad = " << 1+k-brf.Root() << endl;

	   f[i]->SetLineColor(i+1);
	   const int nLinesToDraw = 10;
	   if (i==0){
		   f[i]->Draw("A");
		   leg->AddEntry(f[i],("k = " + to_string(k)).c_str());
	   }
	   else if (i<=nLinesToDraw){
		   f[i]->Draw("same");
		   leg->AddEntry(f[i],("k = " + to_string(k)).c_str());
		   if (i==nLinesToDraw || i==(nPoints-1)){
			   leg->Draw("same");
			   can->SaveAs("test2.png");
		   }
	   }
	   // if (i==0){
	   //         can->Print("test2.pdf[");
	   //         can->SaveAs("test2.png");
	   // }
	   // else if (i==(nPoints-1))
	   //         can->Print("test2.pdf]");
	   // else
	   //         can->Print("test2.pdf");
   } 

   return 0;
}

// int main(){
//
//         NumericalRootFinder2();
//         return 0;
// }
