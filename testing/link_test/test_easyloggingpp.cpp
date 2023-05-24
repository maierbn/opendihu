
#include <easylogging++.h>

INITIALIZE_EASYLOGGINGPP

int main(int argc, char* argv[]) {
   LOG(INFO) << "My first info log using default logger";
   return 0;
}