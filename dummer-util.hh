// Author: Martin C. Frith 2025
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>
#include <iostream>

#define MACROS_SUCK(X) #X
#define STR(X) MACROS_SUCK(X)

bool isDash(const char *text) {
  return text[0] == '-' && text[1] == 0;
}

std::istream &fail(std::istream &s, const char *message) {
  std::cerr << message << "\n";
  s.setstate(std::ios::failbit);
  return s;
}

std::istream &openFile(std::ifstream &file, const char *name) {
  if (isDash(name)) return std::cin;
  file.open(name);
  if (!file) std::cerr << "can't open file: " << name << "\n";
  return file;
}

int badOpt() {
  std::cerr << "bad option value\n";
  return 1;
}
