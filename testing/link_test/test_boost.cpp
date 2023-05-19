
#include <boost/log/trivial.hpp>
#include <boost/thread.hpp>   
#include <boost/date_time.hpp>   
#include <iostream>

void workerFunc()  
{  
    boost::posix_time::seconds workTime(3);          
    std::cout << "Worker: running" << std::endl;    
      
    // Pretend to do something useful... 
    boost::this_thread::sleep(workTime);          
    std::cout << "Worker: finished" << std::endl;  
}    

int main(int, char*[])
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

    return EXIT_SUCCESS;
}