// example program to run the NGL shower
// usage:
//   ./example [-out output.dat] [-nev 1e6] [-drap 1.0] [-longrec]
//
// the program will create a file containing the dS/dt and dS/dlnkt
// histograms.


#include "EvolGrid.hh"
#include "CmdLine.hh"
#include <iostream>
#include <iomanip>
using namespace std;


int main(int argc, char ** argv) {
  CmdLine cmdline(argc,argv);
  int npts  = cmdline.value("-npoints",5000);
  double lmax  = cmdline.value("-lmax",7.0);
  double xmur = cmdline.value("-xmur", 1.0);
  double xQ = cmdline.value("-xQ", 1.0);
  EvolGrid eg(npts, xmur, xQ);
#ifdef NNET
  string fn = cmdline.value<string>("-model","nn_fits/lnkt_model_as0.118_ca3.0_cf1.3333_nf5.pt");
  cout << "# evolution scale grid initialized using: " << fn << endl;
  EvolGrid eg_nn(fn);
#endif
  for(int i=0; i<npts*100; ++i) {
    double t = double(i)*lmax/(npts*100);
    cout << t << " " << eg.ln_kt(t)
#ifdef NNET
	 << " " << eg_nn.ln_kt(t)
#endif
	 << endl;
  }
}
