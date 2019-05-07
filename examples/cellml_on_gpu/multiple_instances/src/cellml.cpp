#include <iostream>
#include <cstdlib>
#include <unistd.h>

#include "opendihu.h"

int main(int argc, char *argv[])
{
  // 0D sub-cellular model
    // initialize everything, handle arguments and parse settings from input file
  DihuContext settings(argc, argv, false, true);
  
  char userInput;
  std::cout << "Do you wish to attach with (e.g.) gdb? (y/n)";
  std::cin >> userInput;
  if(userInput=='y' || userInput=='Y' || userInput=='j' || userInput=='J')
  {
    LOG(WARNING) << "Programme pausing for 1 min. Use the time to attach.";
    sleep(60);
    LOG(WARNING) << "Attachement time is over.";
  }
  
  TimeSteppingScheme::ExplicitEuler<
    CellmlAdapter<57>
  > equationDiscretized(settings);

  LOG(WARNING) << "reached *.run() in main";
  
  equationDiscretized.run();
  
  return EXIT_SUCCESS;
}
