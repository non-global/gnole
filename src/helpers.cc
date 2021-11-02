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

