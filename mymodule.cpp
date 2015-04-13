/*
 *  mymodule.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// include necessary NEST headers
#include "config.h"
#include "network.h"
#include "model.h"
#include "dynamicloader.h"
#include "genericmodel.h"
#include "generic_connector.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "nestmodule.h"

// include headers with your own stuff
#include "mymodule.h"
#include "stdp_connection_ext.h"
#include "stdp_connection_alpha.h"
#include "stdp_connection_multi.h"
#include "iaf_psc_alpha_ext.h"
#include "iaf_psc_alpha_multi_ext.h"
#include "glif_psc_alpha_multi.h"
#include "iaf_freq_sensor.h"
#include "iaf_freq_sensor_v2.h"
#include "iaf_wsn_hermitian_1.h"
#include "iaf_wsn_hermitian_2.h"
#include "iaf_wsn_alpha.h"

// -- Interface to dynamic module loader ---------------------------------------

/*
 * The dynamic module loader must be able to find your module.
 * You make the module known to the loader by defining an instance of your
 * module class in global scope. This instance must have the name
 *
 * <modulename>_LTX_mod
 *
 * The dynamicloader can then load modulename and search for symbol "mod" in it.
 */

mynest::MyModule mymodule_LTX_mod;

// -- DynModule functions ------------------------------------------------------

mynest::MyModule::MyModule()
  {
#ifdef LINKED_MODULE
     // register this module at the dynamic loader
     // this is needed to allow for linking in this module at compile time
     // all registered modules will be initialized by the main app's dynamic loader
     nest::DynamicLoaderModule::registerLinkedModule(this);
#endif
   }

mynest::MyModule::~MyModule()
   {}

   const std::string mynest::MyModule::name(void) const
   {
     return std::string("My NEST Module"); // Return name of the module
   }

   const std::string mynest::MyModule::commandstring(void) const
   {
     // Instruct the interpreter to load mymodule-init.sli
     return std::string("(mymodule-init) run");
   }

   /* BeginDocumentation
   */

  //-------------------------------------------------------------------------------------

  void mynest::MyModule::init(SLIInterpreter *i, nest::Network*)
  {
    /* Register a neuron or device model.
       Give node type as template argument and the name as second argument.
       The first argument is always a reference to the network.
       Return value is a handle for later unregistration.
    */
    nest::register_model<iaf_psc_alpha_ext>(nest::NestModule::get_network(),
                                        "iaf_psc_alpha_ext");
    nest::register_model<iaf_psc_alpha_multi_ext>(nest::NestModule::get_network(),
                                        "iaf_psc_alpha_multi_ext");
    nest::register_model<glif_psc_alpha_multi>(nest::NestModule::get_network(),
                                        "glif_psc_alpha_multi");
    nest::register_model<iaf_freq_sensor>(nest::NestModule::get_network(),
                                        "iaf_freq_sensor");
    nest::register_model<iaf_freq_sensor_v2>(nest::NestModule::get_network(),
                                        "iaf_freq_sensor_v2");
    nest::register_model<iaf_wsn_hermitian_2>(nest::NestModule::get_network(),
                                        "wsn_hermitian_2");
    nest::register_model<iaf_wsn_hermitian_1>(nest::NestModule::get_network(),
                                        "wsn_hermitian_1");
    nest::register_model<iaf_wsn_alpha>(nest::NestModule::get_network(),
                                        "wsn_alpha");


    /* Register a synapse type.
       Give synapse type as template argument and the name as second argument.
       The first argument is always a reference to the network.
    */
    //nest::register_prototype_connection<DropOddSpikeConnection>(nest::NestModule::get_network(),
                                                       //"drop_odd_synapse");
    nest::register_prototype_connection<STDPConnectionExt>(nest::NestModule::get_network(),
        "stdp_synapse_ext");
    nest::register_prototype_connection<STDPConnectionAlpha>(nest::NestModule::get_network(),
            "stdp_synapse_alpha");
    nest::register_prototype_connection<STDPConnectionMulti>(nest::NestModule::get_network(),
            "stdp_synapse_multi");

    /* Register a SLI function.
       The first argument is the function name for SLI, the second a pointer to
       the function object. If you do not want to overload the function in SLI,
       you do not need to give the mangled name. If you give a mangled name, you
       should define a type trie in the mymodule-init.sli file.
    */
    //i->createcommand("StepPatternConnect_Vi_i_Vi_i_l",
                     //&stepPatternConnect_Vi_i_Vi_i_lFunction);

    /* Register a Topography connection kernel function
     *
     */
    nest::TopologyModule::register_parameter<LaplacianParameter>("laplacian");

  }  // MyModule::init()


