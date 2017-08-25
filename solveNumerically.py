from __future__ import print_function
import ROOT
from array import array
import sys
import math
import os

from ROOT import TCanvas, TGraph, TLegend, TF1
from ROOT import gROOT
from ROOT import 
###############################################################################
thetaPointsToDraw = [35, 90]

###############################################################################
def myFunc(x):
    return x[0]**2
###############################################################################
###############################################################################
#  def main(inXmlFile):
def main():
    c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 750, 800 )
    f = TF1("myFunc", myFunc, -2., 2.)
    f.Draw()
    c1.SaveAs("test.png")
###############################################################################
###############################################################################

if __name__ == "__main__":
        main()
	#  if len (sys.argv) == 1:
	#          print ('Please specify input xml-files!!!')
	#  for iFile in sys.argv[1:]:
	#          main(iFile)
