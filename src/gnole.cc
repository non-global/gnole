// example program to run non global evolution
// usage:
//   ./gnole [-out output.dat] [-nev 1e6] [-drap 1.0] [-seed 0] [-order 0] [-xmur 1.0] [-xQ 1.0] [-p -1] [-nbins 100]
//
// the program will create a file containing the dS/dlnET histogram.


#include "Shower.hh"
#include "Observables.hh"
#include "Parameters.hh"
#include "CmdLine.hh"
#include <string>
#include <sstream>
using namespace std;


int main(int argc, char ** argv) {
  CmdLine cmdline(argc,argv);

  ostringstream header;
  header << "# "+cmdline.command_line()
    +"\n# file created on: "+cmdline.time_stamp();

  // run parameters
  int seed   = cmdline.value("-seed",0);
  header << "# seed = " << seed << endl;
  int order  = cmdline.value("-order",0);
  header << "# order = " << order << endl;
  double dy  = cmdline.value("-drap",1.0).help("full rapidity width of the slice (centered on y=0)");
  header << "# dy = " << dy << endl;
  double nev = cmdline.value("-nev",1e6).help("number of events to run");
  header << "# nev_requested = " << nev << endl;
  double xmur = cmdline.value("-xmur",1.0);
  header << "# xmur = " << xmur << endl;
  double xQ = cmdline.value("-xQ",1.0);
  header << "# xQ = " << xQ << endl;

  // set shower IR cutoff (besides the Landau pole)
  double cutoff = cmdline.value("-lnktmax", 20);
  header << "# lnktmax = " << cutoff << endl;
  set_lnktmax(cutoff);

  // set collinear cutoff
  double etamax = cmdline.value("-etamax", 5);
  header << "# etamax = " << etamax << endl;
  set_rapmax(etamax);

  // set strong coupling
  double alphas = cmdline.value("-as", 0.118).help("Value of the coupling at mu=rts").argname("alphas(rts)");
  header << "# alphas(rts) = " << as << endl;
  set_alphas_at_Q(alphas);

  // decide whether to expand NLL corrections
  bool nll_expanded = cmdline.present("-expand-nll");
  set_nll_expanded(nll_expanded);

  // check whether the observable should be computed in SL approximation
  bool sl_observable = cmdline.present("-sl-obs");
  set_sl_observable(sl_observable);

  // observable
  double p = cmdline.value("-p",-1.0);
  header << "# p = " << p << endl;

  int nbins = cmdline.value("-nbins",100);
  header << "# nbins = " << nbins << endl;
  double obsmax = cmdline.value("-obsmax", cutoff);
  header << "# obsmax = " << obsmax << endl;
  Slice slice(dy, p, nbins, obsmax);

  // output
  string filename  = cmdline.value<string>("-out","output.dat");
  string header = "# "+cmdline.command_line()
    +"\n# file created on: "+cmdline.time_stamp();
#ifdef NNET
  string fn_evl_nn  = cmdline.value<string>("-evol","");
  Shower shower(slice, xmur, xQ, fn_evl_nn, order, header.str(), seed);
#else
  Shower shower(slice, xmur, xQ, order, header.str(), seed);
#endif
  // this makes sure there are no unused options
  // left and also triggers the code needed to produce
  // the output for -h
  cmdline.assert_all_options_used();
  

  shower.run(nev, filename);
}
