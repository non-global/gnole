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
#include <limits>
using namespace std;


int main(int argc, char ** argv) try {
  CmdLine cmdline(argc,argv);

  ostringstream header;
  header << "# "+cmdline.command_line()
    +"\n# file created on: "+cmdline.time_stamp();

  // run parameters
  int seed   = cmdline.value("-seed",0);
  header << "# seed = " << seed << endl;
  int order  = cmdline.value("-order",0).help("logarithmic accuracy: 0 = LL, 1 = expanded NLL");
  header << "# order = " << order << endl;
  header << "# evolution = " << (order == 1 ? "expanded NLL" : "LL") << endl;
  double dy  = cmdline.value("-drap",1.0).help("full rapidity width of the slice (centered on y=0)");
  header << "# dy = " << dy << endl;
  double nev = cmdline.value("-nev",1e6).help("number of events to run");
  header << "# nev_requested = " << nev << endl;
  if (!std::isfinite(nev) || nev < 1.0 || nev != std::floor(nev) ||
      nev > std::numeric_limits<int>::max())
    throw CmdLine::Error("-nev must be a positive integer fitting the existing int event counter");
  double xmur = cmdline.value("-xmur",1.0);
  header << "# xmur = " << xmur << endl;
  double xQ = cmdline.value("-xQ",1.0);
  header << "# xQ = " << xQ << endl;

  // set shower IR cutoff (besides the Landau pole)
  double cutoff = cmdline.value("-lnktmax", 20.0);
  header << "# lnktmax = " << cutoff << endl;
  set_lnktmax(cutoff);

  // set collinear cutoff
  double etamax = cmdline.value("-etamax", 5.0);
  header << "# etamax = " << etamax << endl;
  set_rapmax(etamax);

  // set strong coupling
  double alphas = cmdline.value("-as", 0.118).help("Value of the coupling at mu=rts").argname("alphas(rts)");
  header << "# alphas(rts) = " << alphas << endl;
  set_alphas_at_Q(alphas);

  // flag to run only the hard part of the calculation
  bool hard_only = cmdline.present("-hard-only");
  set_hard_only(hard_only);

  // flag to run only the single insertion part of the calculation
  bool single_insertion_only = cmdline.present("-single-insertion-only");
  set_single_insertion_only(single_insertion_only);

  // flag to run only the double insertion part of the calculation
  bool double_insertion_only = cmdline.present("-double-insertion-only");
  bool double_insertion_gauge_only = cmdline.present("-double-insertion-gauge-only")
    .help("double-insertion gauge current only, retaining the LL baseline; requires order 1");
  bool double_insertion_fermion_only = cmdline.present("-double-insertion-fermion-only")
    .help("double-insertion fermion current only, retaining the LL baseline; requires order 1");
  bool double_insertion_requested = double_insertion_only ||
    double_insertion_gauge_only || double_insertion_fermion_only;
  set_double_insertion_only(double_insertion_requested);
  set_double_insertion_gauge_only(double_insertion_gauge_only);
  set_double_insertion_fermion_only(double_insertion_fermion_only);

  // flag to perform NLO matching to QCD
  bool matching = cmdline.present("-matching");
  set_matching(matching);

  // check whether the observable should be computed in SL approximation
  bool obs_correction_only = cmdline.present("-obs-correction-only")
    .help("LL soft observable correction only; implies SL measurement and requires order 0, p=-1");
  bool sl_observable = cmdline.present("-sl-obs");
  if (obs_correction_only) sl_observable = true;
  set_sl_observable(sl_observable);

  VarianceOptions variance;
  string variance_mode = cmdline.value<string>("-variance-mode", "coupled")
    .help("variance reduction: legacy, merge, or coupled");
  if (variance_mode == "legacy") variance.mode = VarianceOptions::Legacy;
  else if (variance_mode == "merge") variance.mode = VarianceOptions::Merge;
  else if (variance_mode == "coupled") variance.mode = VarianceOptions::Coupled;
  else throw CmdLine::Error("-variance-mode must be legacy, merge, or coupled");
  variance.adaptive_insertion = cmdline.present("-adaptive-insertion");
  variance.insertion_power = cmdline.value("-insertion-power", 5.0)
    .help("first-insertion importance-map power (must exceed 2.5)");
  variance.second_power = cmdline.value("-second-power", 1.0)
    .help("second-insertion importance-map power (must exceed 0.5; default 1)");
  bool mc_errors = cmdline.present("-mc-errors");
  bool reference_kernels = cmdline.present("-reference-kernels");
  set_precise_kernels(!reference_kernels);
  variance.optimize = !reference_kernels;
  string replay_mode = cmdline.value<string>("-rng-replay", "tape")
    .help("correlated random-number replay: snapshots or tape");
  if (replay_mode != "snapshots" && replay_mode != "tape")
    throw CmdLine::Error("-rng-replay must be snapshots or tape");
  variance.replay_tape = replay_mode == "tape";
  variance.diagnostics = cmdline.present("-diagnostics");
  variance.output_precision = cmdline.value("-output-precision", 17)
    .help("significant output digits (1..17)");
  if (variance.output_precision < 1 || variance.output_precision > 17)
    throw CmdLine::Error("-output-precision must be in 1..17");
  header << "# precise_kernels = " << PRECISE_KERNELS << endl;
  header << "# rng_replay = " << replay_mode << endl;
  header << "# output_precision = " << variance.output_precision << endl;
  // Avoid introducing divergent proposal second moments into otherwise
  // bounded baseline/second-insertion contributions. Low powers remain useful
  // for algebraic unit probes but are not supported production samplers.
  if (!(variance.insertion_power > 2.5) || !std::isfinite(variance.insertion_power) ||
      !(variance.second_power > 0.5) || !std::isfinite(variance.second_power))
    throw CmdLine::Error("Require finite -insertion-power > 2.5 and -second-power > 0.5");
  header << "# variance_mode = " << variance_mode << endl;
  header << "# adaptive_insertion = " << variance.adaptive_insertion << endl;
  header << "# insertion_power = " << variance.insertion_power << endl;
  header << "# second_power = " << variance.second_power << endl;
  header << "# mc_errors = " << mc_errors << endl;

  // Reject combinations which would silently select a different calculation.
  if (obs_correction_only &&
      (order != 0 || matching || hard_only || single_insertion_only || double_insertion_requested))
    throw CmdLine::Error("-obs-correction-only requires -order 0, without matching or NLL component selectors");
  if (obs_correction_only && variance.adaptive_insertion)
    throw CmdLine::Error("-adaptive-insertion is not used in LL observable-correction mode");
  int selectors = int(hard_only) + int(single_insertion_only) + int(double_insertion_only) +
                  int(double_insertion_gauge_only) + int(double_insertion_fermion_only);
  if (selectors > 1 || (selectors && order != 1))
    throw CmdLine::Error("Use at most one component selector, with -order 1");
  if (matching && (order != 1 || single_insertion_only || double_insertion_requested))
    throw CmdLine::Error("Matching requires order one and the hard contribution");
  if (double_insertion_gauge_only)
    header << "# double_insertion_current = gauge (LL baseline retained)" << endl;
  if (double_insertion_fermion_only)
    header << "# double_insertion_current = fermion (LL baseline retained)" << endl;
  if (matching && sl_observable)
    header << "# matching with SL measurement: endpoint deposits are not physical fixed-order transverse-energy matching" << endl;
  if (order != 0 && order != 1) throw CmdLine::Error("-order must be 0 or 1");
  if (!(xmur > 0.0 && xQ > 0.0 && alphas > 0.0 && cutoff > 0.0 && etamax > 0.0) ||
      !std::isfinite(alphas2(xmur)) || !(alphas2(xmur) > 0.0) ||
      !(std::log(xQ) < 1.0/(2.0*alphas2(xmur)*b0)))
    throw CmdLine::Error("Invalid scale, coupling or cutoff");

  // observable
  double p = cmdline.value("-p",-1.0);
  header << "# p = " << p << endl;

  int nbins = cmdline.value("-nbins",100);
  header << "# nbins = " << nbins << endl;
  double obsmax = cmdline.value("-obsmax", cutoff);
  header << "# obsmax = " << obsmax << endl;
  if (nbins <= 0 || !(obsmax > 0.0) || !(dy > 0.0) || !(p == -1.0 || p > 0.0))
    throw CmdLine::Error("Invalid histogram, slice width or logarithm parameter");
  if (obs_correction_only && p != -1.0)
    throw CmdLine::Error("-obs-correction-only currently requires -p -1 (ordinary logarithms)");
  if (obs_correction_only) {
    header << "# observable_correction_only = 1" << endl;
    header << "# sl_observable = 1 (automatic for observable correction)" << endl;
  }
  Slice slice(dy, p, nbins, obsmax);
  if (obs_correction_only) slice.enable_obs_correction_only(true);
  slice.enable_mc_errors(mc_errors);
  slice.enable_event_score(variance.adaptive_insertion);

  // output
  string filename  = cmdline.value<string>("-out","output.dat");
  Shower shower(slice, xmur, xQ, order, header.str(), seed, variance);
  // this makes sure there are no unused options
  // left and also triggers the code needed to produce
  // the output for -h
  cmdline.assert_all_options_used();
  

  shower.run(static_cast<int>(nev), filename);
} catch (CmdLine::Error&) {
  return 1; // CmdLine already prints its diagnostic.
} catch (const std::exception& error) {
  std::cerr << error.what() << std::endl;
  return 1;
}
