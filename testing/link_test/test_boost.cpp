
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>   
#include <boost/date_time.hpp>   
#include <boost/system/error_code.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

using namespace boost::system;
using namespace boost::program_options;
using namespace boost::filesystem;

void workerFunc()  
{  
    boost::posix_time::seconds workTime(3);          
    std::cout << "Worker: running" << std::endl;    
      
    // Pretend to do something useful... 
    boost::this_thread::sleep(workTime);          
    std::cout << "Worker: finished" << std::endl;  
}    

void fail(error_code &ec)
{
  ec = errc::make_error_code(errc::not_supported);
}

int main(int argc, char** argv)
{
    /* Code to test log module. */
    BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
    BOOST_LOG_TRIVIAL(debug) << "A debug severity message";
    BOOST_LOG_TRIVIAL(info) << "An informational severity message";
    BOOST_LOG_TRIVIAL(warning) << "A warning severity message";
    BOOST_LOG_TRIVIAL(error) << "An error severity message";
    BOOST_LOG_TRIVIAL(fatal) << "A fatal severity message";

    /* Code to test thread module. */
    std::cout << "main: startup" << std::endl;          
    boost::thread workerThread(workerFunc);      
    std::cout << "main: waiting for thread" << std::endl;          
    workerThread.join();
    std::cout << "main: done" << std::endl;

    /* Code to test system module. */
    error_code ec;
    fail(ec);
    boost::system::error_condition ecnd = ec.default_error_condition();
    std::cout << ecnd.value() << '\n';

    /* Code to test program_options module. */
    options_description desc{"Options"};
    variables_map vm;

    /* Code to test filesystem. */
    std::cout << argv[0] << " " << file_size(argv[0]) << '\n';

    return EXIT_SUCCESS;
}