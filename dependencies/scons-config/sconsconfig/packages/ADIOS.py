import sys, os, multiprocessing
from Package import Package
import subprocess

class ADIOS(Package):

  def __init__(self, **kwargs):
    defaults = {
        'download_url': 'https://github.com/ornladios/ADIOS2/archive/v2.3.1.zip',
    }
    defaults.update(kwargs)
    super(ADIOS, self).__init__(**defaults)
    self.ext = '.cpp'
    #self.sub_dirs = [
    #    ('include/mysql', 'lib'),
    #    ('include/mysql', 'lib64'),
    #]
    #self.headers = ['mysql.h']
    #self.libs = ['adios2']
    #self.extra_libs = ['ADIOS', 'blas']
    #self.set_rpath = False    # do not use dynamic linkage
    self.check_text = r'''/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * helloBPWriter.cpp: Simple self-descriptive example of how to write a variable
 * to a BP File that lives in several MPI processes.
 *
 *  Created on: Feb 16, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <mpi.h>
#include <stdexcept> //std::invalid_argument std::exception
#include <string>
#include <vector>

#include <adios2.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /** Application variable */
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const std::size_t Nx = myFloats.size();

    try
    {
        /** ADIOS class factory of IO class objects, DebugON is recommended */
        adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");

        /** global array : name, { shape (total) }, { start (local) }, { count
         * (local) }, all are constant dimensions */
        adios2::Variable<float> bpFloats = bpIO.DefineVariable<float>(
            "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

        bpIO.DefineAttribute<std::string>("Single_String",
                                          "File generated with ADIOS2");

        std::vector<std::string> myStrings = {"one", "two", "three"};
        bpIO.DefineAttribute<std::string>("Array_of_Strings", myStrings.data(),
                                          myStrings.size());

        bpIO.DefineAttribute<double>("Attr_Double", 0.f);
        std::vector<double> myDoubles = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        bpIO.DefineAttribute<double>("Array_of_Doubles", myDoubles.data(),
                                     myDoubles.size());

        /** Engine derived class, spawned to start IO operations */
        adios2::Engine bpWriter =
            bpIO.Open("fileAttributes.bp", adios2::Mode::Write);

        /** Write variable for buffering */
        bpWriter.Put<float>(bpFloats, myFloats.data());

        /** Create bp file, engine becomes unreachable after this*/
        bpWriter.Close();

        adios2::IO bpReader = adios.DeclareIO("BPReader");

        adios2::Engine bpReaderEngine =
            bpReader.Open("fileAttributes.bp", adios2::Mode::Read);

        const auto attributesInfo = bpReader.AvailableAttributes();

        for (const auto &attributeInfoPair : attributesInfo)
        {
            std::cout << "Attribute: " << attributeInfoPair.first;
            for (const auto &attributePair : attributeInfoPair.second)
            {
                std::cout << "\tKey: " << attributePair.first
                          << "\tValue: " << attributePair.second << "\n";
            }
            std::cout << "\n";
        }

        bpReaderEngine.Close();
    }
    catch (std::invalid_argument &e)
    {
        std::cout << "Invalid argument exception, STOPPING PROGRAM from rank "
                  << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::ios_base::failure &e)
    {
        std::cout
            << "IO System base failure exception, STOPPING PROGRAM from rank "
            << rank << "\n";
        std::cout << e.what() << "\n";
    }
    catch (std::exception &e)
    {
        std::cout << "Exception, STOPPING PROGRAM from rank " << rank << "\n";
        std::cout << e.what() << "\n";
    }

    MPI_Finalize();

    return 0;
}
'''

    self.number_output_lines = 4626
      
    self.libs = ["adios2"]
    self.headers = ["adios2.h"]

  def check(self, ctx):
    env = ctx.env
    ctx.Message('Checking for ADIOS ...         ')
    self.check_options(env)

    # reference blas, cmake based, dynamic libraries
    self.set_build_handler([
      'mkdir -p ${PREFIX}',
      'cd ${SOURCE_DIR} && mkdir -p build && cd build && '+ctx.env["cmake"]+' -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_BUILD_TYPE=RELEASE -DADIOS2_USE_SST=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_BUILD_TESTING=OFF -DADIOS2_BUILD_EXAMPLE=OFF \
      -DCMAKE_C_COMPILER='+ctx.env["cc"]+' -DCMAKE_CXX_COMPILER='+ctx.env["CC"]+' ..',
      'cd ${SOURCE_DIR}/build && make all install'
    ])
    
    res = super(ADIOS, self).check(ctx)

    self.check_required(res[0], ctx)
    ctx.Result(res[0])
    return res[0]
