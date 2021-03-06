/*
 *  iaf_freq_sensor_v2.cpp
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

#include "exceptions.h"
#include "iaf_freq_sensor_v2.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

#include <limits>

nest::RecordablesMap<mynest::iaf_freq_sensor_v2> mynest::iaf_freq_sensor_v2::recordablesMap_;
using namespace nest;

namespace nest
{

  /*
   * Override the create() method with one call to RecordablesMap::insert_()
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<mynest::iaf_freq_sensor_v2>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &mynest::iaf_freq_sensor_v2::get_V_m_);
    insert_("V0"      , &mynest::iaf_freq_sensor_v2::get_V0_);
    insert_("V1"      , &mynest::iaf_freq_sensor_v2::get_V1_);
    //insert_("Var"     , &mynest::iaf_freq_sensor_v2::get_Var_);
    insert_("Syn"     , &mynest::iaf_freq_sensor_v2::get_Syn_);
    insert_("Ie"      , &mynest::iaf_freq_sensor_v2::get_Ie_);
    insert_("Currents", &mynest::iaf_freq_sensor_v2::get_Currents_);

  }
}

namespace mynest
{
  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */

  mynest::iaf_freq_sensor_v2::Parameters_::Parameters_()
    : Tau_       (  30.0   ),  // ms
      C_         (   1.0   ),  // pF
      TauR_      (   2.0   ),  // ms
      U0_        (   0.0   ),  // mV
      I_e_       (   0.0   ),  // pA
      V_reset_   ( -10.0   ),  // mV, rel to U0_
      Theta_     (  -1.0   ),  // mV, rel to U0_
      LowerBound_(-std::numeric_limits<double_t>::infinity()),
      Sigma_     (  30.0   ),
      D_Int_     (  0.0   )   // ms
      //Var_Alpha_ (  0.0   )
      //num_of_receptors_ ( 2 )
  {}

  mynest::iaf_freq_sensor_v2::State_::State_()
    : u_    (0.0),
      v0_   (0.0),
      v1_   (0.0),
      //var_  (1.0),
      s_    (0.0),
      Ie_   (0.0),
      currents_ (0.0),
      t_clk_   (-std::numeric_limits<double_t>::infinity()),
      r_    (0)
  {}

  /* ----------------------------------------------------------------
   * Parameter and state extractions and manipulation functions
   * ---------------------------------------------------------------- */

  void mynest::iaf_freq_sensor_v2::Parameters_::get(DictionaryDatum &d) const
  {
    def<double>(d, names::E_L, U0_);   // Resting potential
    def<double>(d, names::I_e, I_e_);
    def<double>(d, names::V_th, Theta_+U0_); // threshold value
    def<double>(d, names::V_reset, V_reset_+U0_);
    def<double>(d, names::V_min, LowerBound_+U0_);
    def<double>(d, names::C_m, C_);
    def<double>(d, names::tau_m, Tau_);
    def<double>(d, names::t_ref, TauR_);
    def<double>(d, "Sigma", Sigma_);
    def<double>(d, "D_Int", D_Int_);
    //def<double>(d, "VarRate", Var_Alpha_);
    //def<int>(d, "n_receptors", num_of_receptors_);
  }

  double mynest::iaf_freq_sensor_v2::Parameters_::set(const DictionaryDatum& d)
  {
    // if U0_ is changed, we need to adjust all variables defined relative to U0_
    const double ELold = U0_;
    updateValue<double>(d, names::E_L, U0_);
    const double delta_EL = U0_ - ELold;

    updateValue<double>(d, names::V_reset, V_reset_);
    updateValue<double>(d, names::V_th, Theta_);
    updateValue<double>(d, names::V_min, LowerBound_);

    updateValue<double>(d, names::I_e, I_e_);
    updateValue<double>(d, names::C_m, C_);
    updateValue<double>(d, names::tau_m, Tau_);
    updateValue<double>(d, names::t_ref, TauR_);
    updateValue<double>(d, "Sigma", Sigma_);
    updateValue<double>(d, "D_Int", D_Int_);
    //updateValue<double>(d, "VarRate", Var_Alpha_);


    if ( C_ <= 0.0 )
      throw BadProperty("Capacitance must be > 0.");

    if ( Tau_ <= 0.0 )
      throw BadProperty("Membrane time constant must be > 0.");

    if ( D_Int_ < 0.0 )
        throw BadProperty("Integration time must be >= 0.");

    if ( TauR_ < 0.0 )
    	throw BadProperty("The refractory time t_ref can't be negative.");

    //if ( Var_Alpha_ < 0.0 || Var_Alpha_ > 1.0 )
        //throw BadProperty("The VarRate should between 0.0 and 1.0");

    return delta_EL;
  }

  void mynest::iaf_freq_sensor_v2::State_::get(DictionaryDatum &d, const Parameters_& p) const
  {
    def<double>(d, names::V_m, u_); // Membrane potential
  }

  void mynest::iaf_freq_sensor_v2::State_::set(const DictionaryDatum& d, const Parameters_& p, double delta_EL)
  {
    updateValue<double>(d, names::V_m, u_);
  }

  mynest::iaf_freq_sensor_v2::Buffers_::Buffers_(iaf_freq_sensor_v2& n)
    : logger_(n)
  {}

  mynest::iaf_freq_sensor_v2::Buffers_::Buffers_(const Buffers_ &, iaf_freq_sensor_v2& n)
    : logger_(n)
  {}


  /* ----------------------------------------------------------------
   * Default and copy constructor for node
   * ---------------------------------------------------------------- */

  mynest::iaf_freq_sensor_v2::iaf_freq_sensor_v2()
    : Archiving_Node(),
      P_(),
      S_(),
      B_(*this)
  {
    recordablesMap_.create();
  }

  mynest::iaf_freq_sensor_v2::iaf_freq_sensor_v2(const iaf_freq_sensor_v2& n)
    : Archiving_Node(n),
      P_(n.P_),
      S_(n.S_),
      B_(n.B_, *this)
  {}

  /* ----------------------------------------------------------------
   * Node initialization functions
   * ---------------------------------------------------------------- */

  void mynest::iaf_freq_sensor_v2::init_state_(const Node& proto)
  {
    const iaf_freq_sensor_v2& pr = downcast<iaf_freq_sensor_v2>(proto);
    S_ = pr.S_;
  }

  void mynest::iaf_freq_sensor_v2::init_buffers_()
  {
    B_.spikes_.clear();       // includes resize
    B_.currents_.clear();        // includes resize

    B_.logger_.reset();

    Archiving_Node::clear_history();
  }

  void mynest::iaf_freq_sensor_v2::calibrate()
  {
    B_.logger_.init();  // ensures initialization in case mm connected after Simulate

    const double h = Time::get_resolution().get_ms();

    //P_.receptor_types_.resize(P_.num_of_receptors_);
    //B_.currents_.resize(P_.num_of_receptors_);
    //for (size_t i = 0; i < P_.num_of_receptors_; i++)
    //{
        //P_.receptor_types_[i] = i+1;
        //B_.currents_[i].resize();
    //}

    // these P are independent
    V_.P2_ = 2.0 / (std::sqrt(3.0 * P_.Sigma_) * std::pow(numerics::pi, 0.25) * P_.Sigma_);
    V_.P32_ = P_.Tau_ / P_.C_ * ( 1.0 -  std::exp(-h / P_.Tau_) );
    V_.P33_ = std::exp(-h / P_.Tau_);

    // TauR specifies the length of the absolute refractory period as
    // a double_t in ms. The grid based iaf_freq_sensor_v2 can only handle refractory
    // periods that are integer multiples of the computation step size (h).
    // To ensure consistency with the overall simulation scheme such conversion
    // should be carried out via objects of class nest::Time. The conversion
    // requires 2 steps:
    //     1. A time object is constructed defining representation of
    //        TauR in tics. This representation is then converted to computation time
    //        steps again by a strategy defined by class nest::Time.
    //     2. The refractory time in units of steps is read out get_steps(), a member
    //        function of class nest::Time.
    //
    // The definition of the refractory period of the iaf_freq_sensor_v2 is consistent
    // the one of iaf_freq_sensor_v2_ps.
    //
    // Choosing a TauR that is not an integer multiple of the computation time
    // step h will lead to accurate (up to the resolution h) and self-consistent
    // results. However, a neuron model capable of operating with real valued spike
    // time may exhibit a different effective refractory time.

    V_.RefractoryCounts_ = Time(Time::ms(P_.TauR_)).get_steps();
    assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
  }

  /* ----------------------------------------------------------------
   * Update and spike handling functions
   */

  inline
  double_t mynest::iaf_freq_sensor_v2::get_Im_(const double_t dt)
  {
      double_t tt = dt*dt / (P_.Sigma_* P_.Sigma_);
      return V_.P2_ * (1.0-tt) * std::exp(-tt/2.0) * S_.Ie_;
  }

  void mynest::iaf_freq_sensor_v2::update(Time const & origin, const long_t from, const long_t to)
  {
    assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
    assert(from < to);

    double t,dt;
    const double h = Time::get_resolution().get_ms();
    double Vm0;

    for ( long_t lag = from ; lag < to ; ++lag )
    {
      t = Time(Time::step(origin.get_steps()+lag+1)).get_ms();
      if(B_.spikes_.get_value(lag) > 0.1)
      {
          // Clock input, reset v0_, v1_, s_, var_, ti_
          S_.v1_ = std::abs(S_.v0_);
          S_.v0_ = 0.0;
          S_.u_ = P_.V_reset_;
          S_.s_ = 1.0;
          S_.t_clk_ = t;
      }
      Vm0 = S_.u_;

      if ( S_.r_ == 0 )
      {
        // neuron not refractory
        S_.u_ = V_.P32_ * S_.s_ * S_.v1_ + V_.P33_ * S_.u_;

        //Simpson's method for v integration
        dt = t - S_.t_clk_ - P_.D_Int_;
        S_.v0_ += (get_Im_(dt) + 4.0 * get_Im_(dt+h/2.0) + get_Im_(dt+h)) * h / 6.0;

        //S_.var_ = (1.0-P_.Var_Alpha_) * S_.var_ + P_.Var_Alpha_ * S_.Ie_ * S_.Ie_;

        S_.s_ = V_.P33_ * S_.s_;

        // lower bound of membrane potential
        S_.u_ = ( S_.u_ < P_.LowerBound_ ? P_.LowerBound_ : S_.u_);

      }
      else // neuron is absolute refractory
        --S_.r_;

      if ( Vm0 < P_.Theta_ && S_.u_ >= P_.Theta_)
      {
        S_.r_  = V_.RefractoryCounts_;
        S_.u_  = P_.V_reset_;
        S_.s_  = 0.0;
        S_.v1_ = 0.0;
        // A supra-threshold membrane potential should never be observable.
        // The reset at the time of threshold crossing enables accurate integration
        // independent of the computation step size, see [2,3] for details.

        set_spiketime(Time::step(origin.get_steps()+lag+1));
        SpikeEvent se;
        network()->send(*this, se, lag);
      }

      // set new input current
      S_.Ie_ = B_.currents_.get_value(lag);

      // log state data
      B_.logger_.record_data(origin.get_steps() + lag);
    }
  }

  //port mynest::iaf_freq_sensor_v2::connect_sender(SpikeEvent&, port receptor_type)
  //{
      //if (receptor_type <=0 || receptor_type > static_cast <port>(P_.num_of_receptors_))
          //throw IncompatibleReceptorType(receptor_type, get_name(), "SpikeEvent");
      //return receptor_type;
  //}

  void mynest::iaf_freq_sensor_v2::handle(SpikeEvent& e)
  {
    assert(e.get_delay() > 0);
    B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
                         e.get_weight() * e.get_multiplicity());
  }

  //port mynest::iaf_freq_sensor_v2::connect_sender(CurrentEvent&, port receptor_type)
  //{
      //if (receptor_type <=0 || receptor_type > static_cast <port>(P_.num_of_receptors_))
          //throw IncompatibleReceptorType(receptor_type, get_name(), "CurrentEvent");
      //return receptor_type;
  //}

  //void mynest::iaf_freq_sensor_v2::handle(CurrentEvent& e)
  //{
    //assert(e.get_delay() > 0);
    //for(size_t i=0;i<P_.num_of_receptors_;++i)
    //{
        //if (P_.receptor_types_[i] == e.get_rport()){
            //const double_t I = e.get_current();
            //const double_t w = e.get_weight();
            //B_.currents_[i].add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), w * I);
        //}
    //}

  //}

  void mynest::iaf_freq_sensor_v2::handle(CurrentEvent& e)
  {
    assert(e.get_delay() > 0);
    const double_t I = e.get_current();
    const double_t w = e.get_weight();
    B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), w * I);
  }

  void mynest::iaf_freq_sensor_v2::handle(DataLoggingRequest& e)
  {
    B_.logger_.handle(e);
  }

} // namespace
