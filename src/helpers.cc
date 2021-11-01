#include "helpers.hh"
#include <sstream>

unsigned long periodic_output_n_queries = 0;
unsigned long periodic_output_frequency = 1000;

bool periodic_output() {
  periodic_output_n_queries += 1;

  if (periodic_output_n_queries % periodic_output_frequency == 0) {
    if (periodic_output_n_queries / periodic_output_frequency >= 6) periodic_output_frequency *= 2;
    return true;
  } else {
    return false;
  }
  
}


// //----------------------------------------------------------------------
// string cumul_output(const NewHist & hist) {
//   ostringstream ostr;
//   double from_below = hist.underflow();
//   double total = hist.total_weight();
//   double norm = 1.0/total;
//   for (unsigned i = 0; true; i++) {
//     // a rough error estimate based on the smaller of the cumulative
//     // results from below / above
//     double err = sqrt(min(from_below, total - from_below));
//     ostr << hist.binlo(i) << " " << from_below * norm << " " << err * norm << endl;
//     if (i < hist.size()) {
//       from_below += hist[i];
//     } else {
//       break;
//     }
//   }
//   return ostr.str();
// }
