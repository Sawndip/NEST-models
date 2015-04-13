/*
 *  stdp_connection_multi.cpp
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
#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_connection_multi.h"
#include "event.h"
using namespace nest;

namespace mynest
{


  STDPConnectionMulti::STDPConnectionMulti() :
    ConnectionHetWD(),
    Aplus_(0.0016),
    Aneg_(0.0055),
    tplus_(11.0),
    tneg_(10.0),
    Wmax_(100.0),
    Esyn_(1.0),
    EmitSpk_(true)
  { }


  STDPConnectionMulti::STDPConnectionMulti(const STDPConnectionMulti &rhs) :
    ConnectionHetWD(rhs)
  {
    Aplus_  = rhs.Aplus_;
    Aneg_   = rhs.Aneg_;
    tplus_  = rhs.tplus_;
    tneg_   = rhs.tneg_;
    Wmax_   = rhs.Wmax_;
    Esyn_   = rhs.Esyn_;
    EmitSpk_= rhs.EmitSpk_;
  }

  void STDPConnectionMulti::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, "Aplus", Aplus_);
    def<double_t>(d, "Aneg", Aneg_);
    def<double_t>(d, "tplus", tplus_);
    def<double_t>(d, "tneg", tneg_);
    def<double_t>(d, "Esyn", Esyn_);
    def<double_t>(d, "Wmax", Wmax_);
    def<bool>    (d, "EmitSpk", EmitSpk_);
  }

  void STDPConnectionMulti::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    updateValue<double_t>(d, "Aplus", Aplus_);
    updateValue<double_t>(d, "Aneg", Aneg_);
    updateValue<double_t>(d, "tplus", tplus_);
    updateValue<double_t>(d, "tneg", tneg_);
    updateValue<double_t>(d, "Esyn", Esyn_);
    updateValue<double_t>(d, "Wmax", Wmax_);
    updateValue<bool>    (d, "EmitSpk", EmitSpk_);
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void STDPConnectionMulti::set_status(const DictionaryDatum & d, nest::index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    set_property<double_t>(d, "Aplus"   , p, Aplus_);
    set_property<double_t>(d, "Aneg"    , p, Aneg_);
    set_property<double_t>(d, "tplus"   , p, tplus_);
    set_property<double_t>(d, "tneg"    , p, tneg_);
    set_property<double_t>(d, "Esyn"    , p, Esyn_);
    set_property<double_t>(d, "Wmax"    , p, Wmax_);
    set_property<bool>    (d, "EmitSpk" , p, EmitSpk_);
  }

  void STDPConnectionMulti::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "Aplus"   );
    initialize_property_array(d, "Aneg"    );
    initialize_property_array(d, "tplus"   );
    initialize_property_array(d, "tneg"    );
    initialize_property_array(d, "Esyn"    );
    initialize_property_array(d, "Wmax"    );
    initialize_property_array(d, "EmitSpk" );
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPConnectionMulti::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property<double_t>(d, "Aplus", Aplus_);
    append_property<double_t>(d, "Aneg", Aneg_);
    append_property<double_t>(d, "tplus", tplus_);
    append_property<double_t>(d, "tneg", tneg_);
    append_property<double_t>(d, "Esyn", Esyn_);
    append_property<double_t>(d, "Wmax", Wmax_);
    append_property<bool>    (d, "EmitSpk", EmitSpk_);
  }

} // of namespace nest
