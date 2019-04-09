#include <iostream>
#include <iomanip>
#include <cstring>

#include "IF97.h"
#include "unit_tests.h"

int main(int argc, char *argv[])
{
  std::cout << "Hello, IF97" << std::endl;

  if ((argc == 2) && (strcmp(argv[1], "-test_all") == 0))
  {
    Unit_Test_All();
    return 0;
  }
  else if ((argc == 3) && (strcmp(argv[1], "-unit_test") == 0))
  {
    std::string str(argv[2]);
    Unit_Test(str);
    return 0;
  }
  else
  {
    std::cerr << "I am afraid I do not understand the command line." << std::endl;
    return 1;
  }

  return 0;
}
