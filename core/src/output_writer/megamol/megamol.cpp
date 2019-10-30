#include "output_writer/megamol/megamol.h"

#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

#include "easylogging++.h"
#include <utility/python_utility.h>

#ifdef HAVE_MEGAMOL
#include <libzmq/zmq.hpp>
#endif
namespace OutputWriter
{

#ifdef HAVE_ADIOS
std::map<std::string,std::array<std::shared_ptr<MegaMol::adios_writer_t>, 2>> MegaMol::adiosWriters_;
std::map<std::string,std::shared_ptr<MegaMol::adios_writer_t>> MegaMol::adiosWriter_;
#endif

BoundingBox::BoundingBox():
  min(Vec3({0.0,0.0,0.0})),
  max(Vec3({0.0,0.0,0.0}))
{

}

#ifdef HAVE_ADIOS
MegaMol::MegaMol(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset) :
  Generic(context, settings, rankSubset), currentOpenWriterIndex_(0)
{
  combineNInstances_ = specificSettings_.getOptionInt("combineNInstances", 1);
  useFrontBackBuffer_ = specificSettings_.getOptionBool("useFrontBackBuffer", true);
}
#else

MegaMol::MegaMol(DihuContext context, PythonConfig settings, std::shared_ptr<Partition::RankSubset> rankSubset) :
Generic(context, settings, rankSubset)
{
}
#endif

#if defined(HAVE_MEGAMOL) && defined(HAVE_ADIOS)

void MegaMol::notifyMegaMol()
{
  if (DihuContext::ownRankNoCommWorld() == 0 && context_.zmqSocket())
  {
    std::stringstream message;
    //message << "return mmHelp()";

    std::stringstream adiosOutputFilename;
    adiosOutputFilename << lastFilename_ << ".bp";
    message << "mmSetParamValue(\"::dat::filename\", \"" << adiosOutputFilename.str() << "\")";

    std::string messageStr = message.str();
    zmq::message_t messageZmq(messageStr.length());

    memcpy(messageZmq.data(), messageStr.data(), messageStr.length());

    bool MegaMolReady = false;
    while(!MegaMolReady)
    {

      LOG(DEBUG) << "send message \"" << messageStr << "\" to MegaMol";
      assert(context_.zmqSocket());
      context_.zmqSocket()->send(messageZmq);

      LOG(DEBUG) << "recv message from MegaMol";
      zmq::message_t receivedMessage(1000);
      int replySize = context_.zmqSocket()->recv(&receivedMessage);

      if (replySize != 0)
      {
        std::string receivedMessageStr;
        receivedMessageStr = (char *)receivedMessage.data();

        LOG(DEBUG) << "received message: \"" << receivedMessageStr << "\", replySize: " << replySize;
        if (receivedMessageStr.find("READY") == std::string::npos)
        {
          LOG(DEBUG) << "MegaMol accepted request.";
          MegaMolReady = true;
        }
        else
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
      }
      else
      {
        LOG(DEBUG) << "replySize: " << replySize;
      }
    }
  }
}

#endif

} // namespace
