/*
 *  stdp_connection_ext.cpp
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
#include "stdp_connection_ext.h"
#include "event.h"
using namespace nest;

namespace mynest
{

  STDPConnectionExt::STDPConnectionExt() :
    ConnectionHetWD(),
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0),
    Kplus_(0.0),
    Gpre_(0.0),
    Gpost_(0.0),
    LearnEn_(true),
    EmitSpk_(true),
    Esyn_(1.0)
  { }


  STDPConnectionExt::STDPConnectionExt(const STDPConnectionExt &rhs) :
    ConnectionHetWD(rhs)
  {
    tau_plus_ = rhs.tau_plus_;
    lambda_ = rhs.lambda_;
    alpha_ = rhs.alpha_;
    mu_plus_ = rhs.mu_plus_;
    mu_minus_ = rhs.mu_minus_;
    Wmax_ = rhs.Wmax_;
    Kplus_ = rhs.Kplus_;
    Gpre_ = rhs.Gpre_;
    Gpost_ = rhs.Gpost_;
    LearnEn_ = rhs.LearnEn_;
    EmitSpk_ = rhs.EmitSpk_;
    Esyn_ = rhs.Esyn_;
  }

  void STDPConnectionExt::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "alpha", alpha_);
    def<double_t>(d, "mu_plus", mu_plus_);
    def<double_t>(d, "mu_minus", mu_minus_);
    def<double_t>(d, "Wmax", Wmax_);
    def<double_t>(d, "Gpre", Gpre_);
    def<double_t>(d, "Gpost", Gpost_);
    def<bool>(d, "LearnEn", LearnEn_);
    def<bool>(d, "EmitSpk", EmitSpk_);
    def<double_t>(d, "Esyn", Esyn_);
  }

  void STDPConnectionExt::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "lambda", lambda_);
    updateValue<double_t>(d, "alpha", alpha_);
    updateValue<double_t>(d, "mu_plus", mu_plus_);
    updateValue<double_t>(d, "mu_minus", mu_minus_);
    updateValue<double_t>(d, "Wmax", Wmax_);
    updateValue<double_t>(d, "Gpre", Gpre_);
    updateValue<double_t>(d, "Gpost", Gpost_);
    updateValue<bool>(d, "LearnEn", LearnEn_);
    updateValue<bool>(d, "EmitSpk", EmitSpk_);
    updateValue<double_t>(d, "Esyn", Esyn_);
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void STDPConnectionExt::set_status(const DictionaryDatum & d, nest::index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "lambdas", p, lambda_);
    set_property<double_t>(d, "alphas", p, alpha_);
    set_property<double_t>(d, "mu_pluss", p, mu_plus_);
    set_property<double_t>(d, "mu_minuss", p, mu_minus_);
    set_property<double_t>(d, "Wmaxs", p, Wmax_);
    set_property<double_t>(d, "Gpres", p, Gpre_);
    set_property<double_t>(d, "Gposts", p, Gpost_);
    set_property<bool>(d, "LearnEns", p, LearnEn_);
    set_property<bool>(d, "EmitSpks", p, EmitSpk_);
    set_property<double_t>(d, "Esyns", p, Esyn_);
  }

  void STDPConnectionExt::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "tau_pluss"); 
    initialize_property_array(d, "lambdas"); 
    initialize_property_array(d, "alphas"); 
    initialize_property_array(d, "mu_pluss"); 
    initialize_property_array(d, "mu_minuss");
    initialize_property_array(d, "Wmaxs");
    initialize_property_array(d, "Gpres");
    initialize_property_array(d, "Gposts");
    initialize_property_array(d, "LearnEns");
    initialize_property_array(d, "EmitSpks");
    initialize_property_array(d, "Esyns");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPConnectionExt::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "tau_pluss", tau_plus_); 
    append_property<double_t>(d, "lambdas", lambda_); 
    append_property<double_t>(d, "alphas", alpha_); 
    append_property<double_t>(d, "mu_pluss", mu_plus_); 
    append_property<double_t>(d, "mu_minuss", mu_minus_);
    append_property<double_t>(d, "Wmaxs", Wmax_);
    append_property<double_t>(d, "Gpres", Gpre_);
    append_property<double_t>(d, "Gposts", Gpost_);
    append_property<bool>(d, "LearnEns", LearnEn_);
    append_property<bool>(d, "EmitSpks", EmitSpk_);
    append_property<double_t>(d, "Esyns", Esyn_);
  }

} // of namespace nest
