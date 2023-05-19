
#include <libxml/xmlwriter.h>
#include <cassert>
#include <iostream>

int main()
{
  xmlTextWriterPtr writer;
  writer = xmlNewTextWriterFilename("xml_output.tmp", 0);
  assert(writer != NULL);
  std::cout << "xml writer success." << std::endl;

  return EXIT_SUCCESS;
}