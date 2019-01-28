#include "control/dihu_context.h"

#include <Python.h>  // this has to be the first included header

#ifdef HAVE_MEGAMOL
#include "Console.h"     // contains megamol_main
#endif

#ifdef HAVE_ADIOS
#include <adios2.h>
#endif

#ifdef HAVE_MEGAMOL
std::shared_ptr<zmq::context_t> DihuContext::zmqContext_;  ///< the 0mq context
std::shared_ptr<zmq::socket_t> DihuContext::zmqSocket_;  ///< a socket that is connected to one megamol
#endif

void DihuContext::initializeAdios(int argc, char *argv[])
{
#ifdef HAVE_ADIOS
  LOG(DEBUG) << "initializeAdios";

  adios_ = std::make_shared<adios2::ADIOS>(MPI_COMM_WORLD);
  io_ = std::make_shared<adios2::IO>(adios_->DeclareIO("Output"));
  assert(io_);
#endif
}

#ifdef HAVE_MEGAMOL
void DihuContext::initializeZMQ()
{
  if (!zmqContext_)
  {
    int ioThreads = 1;
    // The io_threads argument specifies the size of the ØMQ thread pool to handle I/O operations. If your application is using only the inproc transport for messaging you may set this to zero, otherwise set it to at least one.
    zmqContext_ = std::make_shared<zmq::context_t>(ioThreads);
    zmqSocket_ = std::make_shared<zmq::socket_t>(*zmqContext_, ZMQ_PAIR);  // also possible: ZMQ_PAIR
    // types of sockets: https://linux.die.net/man/3/zmq_socket

    zmq::socket_t requesterSocket(*zmqContext_, ZMQ_REQ);
    requesterSocket.connect("tcp://localhost:33333");

    // send message to request port from megamol
    std::string messageStr = "Requesting Port!";
    zmq::message_t messageZmq(messageStr.length());
    memcpy(messageZmq.data(), messageStr.data(), messageStr.length());

    LOG(DEBUG) << "send message \"" << messageStr << "\" to MegaMol";
    requesterSocket.send(messageZmq);

    //
    LOG(DEBUG) << "recv message from MegaMol";
    zmq::message_t receivedMessage(1000);
    int replySize = requesterSocket.recv(&receivedMessage);

    if (replySize > 0)
    {
      std::string receivedMessageStr;
      receivedMessageStr = (char *)receivedMessage.data();

      LOG(DEBUG) << "received message: \"" << receivedMessageStr << "\", replySize: " << replySize;

      int port = atoi(receivedMessageStr.c_str());
      LOG(DEBUG) << "port: " << port;

      std::stringstream address;
      address << "tcp://localhost:" << port;
      zmqSocket_->connect(address.str().c_str());
      LOG(DEBUG) << "connect succeeded";

// test connection by sending a message
#if 0
      std::this_thread::sleep_for(std::chrono::milliseconds(2000));

      messageStr = "return mmListModules()";
      LOG(DEBUG) << "send message \"" << messageStr << "\" to MegaMol";

      zmq::message_t messageZmq(messageStr.length());
      memcpy(messageZmq.data(), messageStr.data(), messageStr.length());
      zmqSocket_->send(messageZmq);

      LOG(DEBUG) << "recv message from MegaMol";
      zmq::message_t receivedMessage(1000);
      int replySize = zmqSocket_->recv(&receivedMessage);

      if (replySize != 0)
      {
        std::string receivedMessageStr;
        receivedMessageStr = (char *)receivedMessage.data();

        LOG(DEBUG) << "received message: \"" << receivedMessageStr << "\", replySize: " << replySize;
      }
      else
      {
        LOG(DEBUG) << "replySize: " << replySize;
      }
#endif
    }
    else
    {
      LOG(ERROR) << "No port received from MegaMol";
    }
  }
}
#endif

void DihuContext::initializeMegaMol(int argc, char *argv[])
{
  initializeAdios(argc, argv);

  // extract MegaMol arguments from config
  if (pythonConfig_.hasKey("MegaMolArguments"))
  {
#ifdef HAVE_MEGAMOL

    // only start MegaMol and initialize ZMQ context if rank is 0
    if (ownRankNo() == 0)
    {

      std::string megamolArgumentsOption = pythonConfig_.getOptionString("MegaMolArguments", "");

      megamolArguments_.clear();
      // convert options from a string to a vector of strings
      size_t pos = 0;
      while (pos < megamolArgumentsOption.length())
      {
        size_t newPos = megamolArgumentsOption.find(" ", pos);
        if (newPos == std::string::npos)
        {
          megamolArguments_.push_back(megamolArgumentsOption.substr(pos));
          break;
        }

        megamolArguments_.push_back(megamolArgumentsOption.substr(pos, newPos-pos));
        pos = megamolArgumentsOption.find(" ", pos+1)+1;
      }

      // prepare arguments
      int megamolArgc = megamolArguments_.size()+1;
      megamolArgv_.resize(megamolArgc);
      static std::string commandName = "mmconsole";
      megamolArgv_[0] = (char *)commandName.c_str();

      for (int i = 0; i < megamolArguments_.size(); i++)
      {
        megamolArgv_[i+1] = (char *)megamolArguments_[i].c_str();
      }

      LOG(DEBUG) << "start MegaMol with arguments " << megamolArguments_;

      //start megamol main method
      megamolThread_ = std::make_shared<std::thread>(megamol_main, megamolArgc, megamolArgv_.data());

      initializeZMQ();
    }
#else
    LOG(ERROR) << "Not compiled with MegaMol, but \"MegaMolArguments\" are given.";
#endif
  }
}
