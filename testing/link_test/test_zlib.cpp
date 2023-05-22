
#include <zlib.h>
#include <iostream>

int main()
{
  z_off_t len = 3000;
  uLong a = 1, b = 2;
  if (a == adler32_combine (a, b, len)) {;}

  return EXIT_SUCCESS;
}