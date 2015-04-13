/*
 *  stdp_connection.cpp
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
#include "stdp_connection_alpha.h"
#include "event.h"
using namespace nest;

namespace mynest
{


  STDPConnectionAlpha::STDPConnectionAlpha() :
    ConnectionHetWD(),
    lambda_(0.01),
    amp_(1.2),
    shift_(-0.2),
    sigma_(1.7),
    center_(-2.0),
    Wmax_(100.0),
    Esyn_(1.0),
    EmitSpk_(true)
  { }


  STDPConnectionAlpha::STDPConnectionAlpha(const STDPConnectionAlpha &rhs) :
    ConnectionHetWD(rhs)
  {
    lambda_ = rhs.lambda_;
    amp_    = rhs.amp_;
    shift_  = rhs.shift_;
    sigma_  = rhs.sigma_;
    center_ = rhs.center_;
    Wmax_   = rhs.Wmax_;
    Esyn_   = rhs.Esyn_;
    EmitSpk_= rhs.EmitSpk_;
  }

  void STDPConnectionAlpha::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "amp", amp_);
    def<double_t>(d, "shift", shift_);
    def<double_t>(d, "sigma", sigma_);
    def<double_t>(d, "center", center_);
    def<double_t>(d, "Wmax", Wmax_);
    def<double_t>(d, "Esyn", Esyn_);
    def<bool>    (d, "EmitSpk", EmitSpk_);
  }

  void STDPConnectionAlpha::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    updateValue<double_t>(d, "lambda" , lambda_);
    updateValue<double_t>(d, "amp"    , amp_);
    updateValue<double_t>(d, "shift"  , shift_);
    updateValue<double_t>(d, "sigma"  , sigma_);
    updateValue<double_t>(d, "center" , center_);
    updateValue<double_t>(d, "Wmax"   , Wmax_);
    updateValue<double_t>(d, "Esyn"   , Esyn_);
    updateValue<bool>    (d, "EmitSpk", EmitSpk_);
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void STDPConnectionAlpha::set_status(const DictionaryDatum & d, nest::index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    set_property<double_t>(d, "lambda" , p, lambda_);
    set_property<double_t>(d, "amp"    , p, amp_);
    set_property<double_t>(d, "shift"  , p, shift_);
    set_property<double_t>(d, "sigma"  , p, sigma_);
    set_property<double_t>(d, "center" , p, center_);
    set_property<double_t>(d, "Wmax"   , p, Wmax_);
    set_property<double_t>(d, "Esyn"   , p, Esyn_);
    set_property<bool>    (d, "EmitSpk", p, EmitSpk_);
  }

  void STDPConnectionAlpha::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "lambda" );
    initialize_property_array(d, "amp"    );
    initialize_property_array(d, "shift"  );
    initialize_property_array(d, "sigma"  );
    initialize_property_array(d, "center" );
    initialize_property_array(d, "Wmax"   );
    initialize_property_array(d, "Esyn"   );
    initialize_property_array(d, "EmitSpk");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPConnectionAlpha::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property<double_t>(d, "lambda" , lambda_);
    append_property<double_t>(d, "amp"    , amp_);
    append_property<double_t>(d, "shift"  , shift_);
    append_property<double_t>(d, "sigma"  , sigma_);
    append_property<double_t>(d, "center" , center_);
    append_property<double_t>(d, "Wmax"   , Wmax_);
    append_property<double_t>(d, "Esyn"   , Esyn_);
    append_property<bool>    (d, "EmitSpk", EmitSpk_);
  }

} // of namespace nest
