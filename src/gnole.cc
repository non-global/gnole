// example program to run non global evolution
// usage:
//   ./gnole [-out output.dat] [-nev 1e6] [-drap 1.0] [-seed 0] [-order 0] [-xmur 1.0] [-xQ 1.0] [-p -1] [-nbins 100]
//
// the program will create a file containing the dS/dlnET histogram.


#include "Shower.hh"
#include "Observables.hh"
#include "CmdLine.hh"
#include "Parameters.hh"
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

  // set shower IR cutoff (besides the Landau pole)
  double cutoff = cmdline.value("-lnktmax", 20);
  set_lnktmax(cutoff);

  // set collinear cutoff
  double etamax = cmdline.value("-etamax", 5);
  set_rapmax(etamax);

  // set strong coupling
  double alphas = cmdline.value("-as", 0.118);
  set_alphas_at_Q(alphas);

  // decide whether to expand NLL corrections
  bool nll_expanded = cmdline.present("-expand-nll");
  set_nll_expanded(nll_expanded);

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
