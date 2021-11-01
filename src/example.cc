// example program to run the NGL shower
// usage:
//   ./example [-out output.dat] [-nev 1e6] [-drap 1.0] [-longrec]
//
// the program will create a file containing the dS/dt and dS/dlnkt
// histograms.


#include "Shower.hh"
#include "Observables.hh"
#include "CmdLine.hh"
#include <string>
using namespace std;


int main(int argc, char ** argv) {
  CmdLine cmdline(argc,argv);

  // run parameters
  int seed   = cmdline.value("-seed",0);
  int order  = cmdline.value("-order",0);
  double dy  = cmdline.value("-drap",1.0);
  double nev = cmdline.value("-nev",1e6);
  double xmur = cmdline.value("-xmur",1.0);
  double xQ = cmdline.value("-xQ",1.0);

  // observable
  double p = cmdline.value("-p",-1.0);
  int nbins = cmdline.value("-nbins",100);
  Slice slice(dy, p, nbins);

  // output
  string filename  = cmdline.value<string>("-out","output.dat");
  string header = "# "+cmdline.command_line()
    +"\n# file created on: "+cmdline.time_stamp();
#ifdef NNET
  string fn_evl_nn  = cmdline.value<string>("-evol","");
  Shower shower(slice, xmur, xQ, fn_evl_nn, order, header, seed);
#else
  Shower shower(slice, xmur, xQ, order, header, seed);
#endif
  shower.run(nev, filename);
}
