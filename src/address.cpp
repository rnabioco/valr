#include "pch.h"
#include <dplyr/main.h>

#include <tools/encoding.h>

using namespace Rcpp;

const char* address(SEXP x) {
  static char buffer[20];
  snprintf(buffer, 20, "%p", reinterpret_cast<void*>(x));
  return (const char*)buffer;
}

