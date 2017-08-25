#include "Math/Functor.h"
#include "Math/RootFinderAlgorithms.h"

double k = 1.0;

double myfunc(double x ) { 
   return k*x*x-1; 
}
 
double myfunc_deriv ( double x ) { 
   return 2*k*x; 
}
 
int NumericalRootFinder()
{

   for (int i=0; i<10; i++){
	   k = 0.1 + 0.1*i;
	   ROOT::Math::GradFunctor1D f1(&myfunc, &myfunc_deriv);
		 
	   // Create the Integrator
	   ROOT::Math::Roots::Newton rf;
	 
	   rf.SetFunction(f1, 2.0); 
	 
	   rf.Solve();
	 
	   cout << "k = " << k << "; root = " << rf.Root() << "; " ;
	   cout << "sqrt(1/k) = " << TMath::Sqrt(1/k) << endl;
   }
   TCanvas *can = new TCanvas("can","can",800,600);
   TF1 *myFunc = new TF1("myfunc", "myfunc(x)", -2.,2.);
   myFunc->Draw();
   can->SaveAs("cppTest.png");

   return 0;
}
