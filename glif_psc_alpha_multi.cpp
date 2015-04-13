/*
 *  glif_psc_alpha_multi.cpp
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
#include "glif_psc_alpha_multi.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

#include <limits>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::glif_psc_alpha_multi> mynest::glif_psc_alpha_multi::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<mynest::glif_psc_alpha_multi>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &mynest::glif_psc_alpha_multi::get_V_m_);
    insert_(names::currents, &mynest::glif_psc_alpha_multi::get_current_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
namespace mynest
{
mynest::glif_psc_alpha_multi::Parameters_::Parameters_()
  : C_                   (  1.0    ),  // pF
    U0_                  (  0.0    ),  // mV
    I_e_                 (  0.0    ),  // pA
    V_reset_             (  0.0    ),  // mV, rel to U0_
    TauR_                (  2.0    ),  // ms
    Theta_               (  4.5    ),  // mV, rel to U0_
    LowerBound_          (-std::numeric_limits<double_t>::infinity()),
    g_L_                 (  0.3    ),
    i_L_                 (  0.0    ),
    num_of_ionchannels_  (   0     ),
    num_of_receptors_    (   0     ),
    has_connections_     ( false   )

{
  A_k_.clear();
  l_k_.clear();
  mu_k_.clear();
  g_k_.clear();
  E_k_.clear();
  tau_syn_r_.clear();  
  tau_syn_f_.clear();
}

mynest::glif_psc_alpha_multi::State_::State_()
  : y0_   (0.0),  
    y3_   (0.0),
    y4_   (0.0),
    r_    (0)
{
  y1_syn_.clear();
  y2_syn_.clear();
}


/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::glif_psc_alpha_multi::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::C_m,     C_);
  def<double>(d, names::E_L,     U0_);         // resting potential
  def<double>(d, names::I_e,     I_e_);
  def<double>(d, names::V_reset, V_reset_);
  def<double>(d, names::t_ref,   TauR_);
  def<double>(d, names::V_th,    Theta_);      // threshold value
  def<double>(d, names::V_min,   LowerBound_);
  def<int>   (d, "n_channels",   num_of_ionchannels_); 
  def<double>(d, "g_L",           g_L_);
  def<double>(d, "i_L",           i_L_);
  def<int>   (d, "n_synapses",   num_of_receptors_);
  def<bool>  (d, names::has_connections, has_connections_);

  ArrayDatum A_k_ad(A_k_);
  ArrayDatum l_k_ad(l_k_);
  ArrayDatum mu_k_ad(mu_k_);
  ArrayDatum g_k_ad(g_k_);
  ArrayDatum E_k_ad(E_k_);
  ArrayDatum tau_syn_r_ad(tau_syn_r_);
  ArrayDatum tau_syn_f_ad(tau_syn_f_);
  def<ArrayDatum>(d,"A_k", A_k_ad);
  def<ArrayDatum>(d,"l_k", l_k_ad);
  def<ArrayDatum>(d,"mu_k", mu_k_ad);
  def<ArrayDatum>(d,"g_k", g_k_ad);
  def<ArrayDatum>(d,"E_k", E_k_ad);
  def<ArrayDatum>(d,"tau_syn_r", tau_syn_r_ad);
  def<ArrayDatum>(d,"tau_syn_f", tau_syn_f_ad);

}

double mynest::glif_psc_alpha_multi::Parameters_::set(const DictionaryDatum& d)
{
   //if U0_ is changed, we need to adjust all variables defined relative to U0_
  const double ELold = U0_;
  //updateValue<double>(d, names::E_L, U0_);
  //const double delta_EL = U0_ - ELold;

  //if(updateValue<double>(d, names::V_reset, V_reset_))
    //V_reset_ -= U0_;
  //else
    //V_reset_ -= delta_EL;

  //if (updateValue<double>(d, names::V_th, Theta_)) 
    //Theta_ -= U0_;
  //else
    //Theta_ -= delta_EL;

  //if (updateValue<double>(d, names::V_min, LowerBound_)) 
    //LowerBound_ -= U0_;
  //else
    //LowerBound_ -= delta_EL;

    
  updateValue<double>(d, names::C_m,   C_);
  bool renew_iL=updateValue<double>(d, names::E_L, U0_);
  updateValue<double>(d, names::I_e,   I_e_);
  updateValue<double>(d, names::V_reset,  V_reset_);
  updateValue<double>(d, names::t_ref, TauR_);
  updateValue<double>(d, names::V_th, Theta_);
  updateValue<double>(d, names::V_min, LowerBound_);
  if(updateValue<double>(d, "g_L", g_L_))
      renew_iL=true;

  const double delta_EL = U0_ - ELold;

  if ( C_ <= 0 )
    throw BadProperty("Capacitance must be > 0.");

  std::vector<double> A_k_tmp;
  std::vector<double> l_k_tmp;
  std::vector<double> mu_k_tmp;
  std::vector<double> g_k_tmp;
  std::vector<double> E_k_tmp;

  bool renew_ion=false;
  if(updateValue<std::vector<double> >(d, "A_k", A_k_tmp))
  {
      renew_ion=true;
      A_k_ = A_k_tmp;
  }
  if(updateValue<std::vector<double> >(d, "l_k", l_k_tmp))
  {
      renew_ion=true;
      l_k_ = l_k_tmp;
  }
  if(updateValue<std::vector<double> >(d, "mu_k", mu_k_tmp))
  {
      renew_ion=true;
      mu_k_ = mu_k_tmp;
  }
  if(updateValue<std::vector<double> >(d, "g_k", g_k_tmp))
  {
      renew_ion = true;
      renew_iL = true;
      g_k_ = g_k_tmp;
  }
  if(updateValue<std::vector<double> >(d, "E_k", E_k_tmp))
  {
      renew_ion=true;
      E_k_ = E_k_tmp;
  }

  if(renew_ion){
      num_of_ionchannels_ = A_k_.size();
      if (    l_k_.size() != num_of_ionchannels_ ||
              mu_k_.size() != num_of_ionchannels_ ||
              g_k_.size() != num_of_ionchannels_ ||
              E_k_.size() != num_of_ionchannels_ )
          throw BadProperty("The parameters for ion channels should be arrays with the same size.");
  }
  
  if(renew_iL)
  {
      double_t gall,iall;
      gall=g_L_;
      iall=0.0;
      for (size_t i = 0; i < num_of_ionchannels_; ++i)
      {
          gall+=g_k_[i];
          iall+=g_k_[i]*E_k_[i];
      }
      i_L_=U0_*gall-iall;
  }
  
  bool renew_tau_syn=false;
  std::vector<double> tau_tmp;
  if (updateValue<std::vector<double> >(d, "tau_syn_r", tau_tmp))
  {
    renew_tau_syn=true;
    for (size_t i = 0; i < tau_tmp.size(); ++i)
    {
      if (tau_tmp.size() < tau_syn_r_.size() && has_connections_ == true)
        throw BadProperty("The neuron has connections, therefore the number of ports cannot be reduced.");
      if (tau_tmp[i] <= 0)
        throw BadProperty("All synaptic time constants must be > 0.");
      //if (tau_tmp[i] == Tau_)
        //throw BadProperty("Membrane and synapse time constant(s) must differ. See note in documentation.");
    }

    tau_syn_r_ = tau_tmp;
  }

  if (updateValue<std::vector<double> >(d, "tau_syn_f", tau_tmp))
  {
    renew_tau_syn=true;
    for (size_t i = 0; i < tau_tmp.size(); ++i)
    {
      if (tau_tmp.size() < tau_syn_f_.size() && has_connections_ == true)
        throw BadProperty("The neuron has connections, therefore the number of ports cannot be reduced.");
      if (tau_tmp[i] <= 0)
        throw BadProperty("All synaptic time constants must be > 0.");
      //if (tau_tmp[i] == Tau_)
        //throw BadProperty("Membrane and synapse time constant(s) must differ. See note in documentation.");
    }

    tau_syn_f_ = tau_tmp;
    //num_of_receptors_ = tau_syn_.size();
  }

  if(renew_tau_syn)
  {
    num_of_receptors_ = tau_syn_r_.size();
    if(tau_syn_f_.size() != num_of_receptors_)
    {
        throw BadProperty("Synaptic time constants should be arrays of the same size");
    }
  }

  if ( TauR_ < 0. )
  	throw BadProperty("The refractory time t_ref can't be negative.");

  if ( V_reset_ >= Theta_ )
    throw BadProperty("Reset potential must be smaller than threshold.");

  return delta_EL;
}

void mynest::glif_psc_alpha_multi::State_::get(DictionaryDatum& d, const Parameters_& p) const
{
  def<double>(d, names::V_m, y3_ ); // Membrane potential
}

void mynest::glif_psc_alpha_multi::State_::set(const DictionaryDatum& d, const Parameters_& p, const double delta_EL)
{
  updateValue<double>(d, names::V_m, y3_);
}

mynest::glif_psc_alpha_multi::Buffers_::Buffers_(glif_psc_alpha_multi& n)
  : logger_(n)
{}

mynest::glif_psc_alpha_multi::Buffers_::Buffers_(const Buffers_ &, glif_psc_alpha_multi& n)
  : logger_(n)
{}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

mynest::glif_psc_alpha_multi::glif_psc_alpha_multi()
  : Archiving_Node(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
}

mynest::glif_psc_alpha_multi::glif_psc_alpha_multi(const glif_psc_alpha_multi& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::glif_psc_alpha_multi::init_state_(const Node& proto)
{
  const glif_psc_alpha_multi& pr = downcast<glif_psc_alpha_multi>(proto);
  S_ = pr.S_;
}

void mynest::glif_psc_alpha_multi::init_buffers_()
{
  B_.spikes_.clear();          // includes resize
  B_.currents_.clear();        // includes resize

  B_.logger_.reset();

  Archiving_Node::clear_history();
}

void mynest::glif_psc_alpha_multi::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate

  const double h = Time::get_resolution().get_ms();

  P_.receptor_types_.resize(P_.num_of_receptors_);
  for (size_t i=0; i < P_.num_of_receptors_; i++)
  {
    P_.receptor_types_[i] = i+1;
  }

  V_.P11_syn_.resize(P_.num_of_receptors_);
  V_.P21_syn_.resize(P_.num_of_receptors_);
  V_.P22_syn_.resize(P_.num_of_receptors_);
  //V_.P31_syn_.resize(P_.num_of_receptors_);
  //V_.P32_syn_.resize(P_.num_of_receptors_);
  V_.P40_.resize(P_.num_of_ionchannels_);
  V_.P44_.resize(P_.num_of_ionchannels_);
  V_.Y40_.resize(P_.num_of_ionchannels_);
  
  S_.y1_syn_.resize(P_.num_of_receptors_);
  S_.y2_syn_.resize(P_.num_of_receptors_);
  S_.y4_.resize(P_.num_of_ionchannels_);
  
  for (size_t i=0; i< P_.num_of_ionchannels_; ++i)
  {
      V_.Y40_[i]=std::exp(P_.mu_k_[i]/P_.l_k_[i]);
      S_.y4_[i]=0.0;
      V_.P40_[i]=P_.A_k_[i]/P_.l_k_[i];
      V_.P44_[i]=std::exp(-h/P_.l_k_[i]);
  }
  V_.minus_h_Cm_ = -h / P_.C_;
  //V_.PSCInitialValues_.resize(P_.num_of_receptors_);

  B_.spikes_.resize(P_.num_of_receptors_);

  //V_.P33_ = std::exp(-h/P_.Tau_);
  //V_.P30_ = 1.0/P_.C_*(1.0-V_.P33_);


  for (size_t i=0; i < P_.num_of_receptors_; i++)
  {
    V_.P11_syn_[i] = std::exp(-h/P_.tau_syn_r_[i]);
    V_.P22_syn_[i] = std::exp(-h/P_.tau_syn_f_[i]);
    V_.P21_syn_[i] = 1.0-V_.P22_syn_[i];

    //V_.P11_syn_[i] = V_.P22_syn_[i] =std::exp(-h/P_.tau_syn_[i]);
    //V_.P21_syn_[i] = h*V_.P11_syn_[i];
    //V_.P31_syn_[i] = 1/P_.C_ * ((V_.P11_syn_[i]-V_.P33_)/(-1/P_.tau_syn_[i]- -1/P_.Tau_)- h*V_.P11_syn_[i])
      ///(-1/P_.Tau_ - -1/P_.tau_syn_[i]);
    //V_.P32_syn_[i] = 1/P_.C_*(V_.P33_-V_.P11_syn_[i])/(-1/P_.Tau_ - -1/P_.tau_syn_[i]);

    //V_.PSCInitialValues_[i] = 1.0 * numerics::e/P_.tau_syn_[i];
    B_.spikes_[i].resize();
  }

  
  Time r=Time::ms(P_.TauR_);
  V_.RefractoryCounts_=r.get_steps();
  
  if ( V_.RefractoryCounts_ < 1 )
    throw BadProperty("Absolute refractory time must be at least one time step.");
}

void mynest::glif_psc_alpha_multi::update(Time const& origin, const long_t from, const long_t to)
{
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    if ( S_.r_ == 0 )
    {
      // neuron not refractory
      //S_.y3_ = V_.P30_*(S_.y0_ + P_.I_e_) + V_.P33_*S_.y3_;
      //S_.current_= 0.0;
      S_.current_ = S_.y0_ + P_.I_e_;
      for (size_t i=0; i < P_.num_of_receptors_; i++){
        //S_.y3_ += V_.P31_syn_[i]*S_.y1_syn_[i] + V_.P32_syn_[i]*S_.y2_syn_[i];
        //S_.y3_ += V_.P30_*S_.y2_syn_[i];
        S_.current_ += S_.y2_syn_[i];
      }

    }
    else
    { // neuron is absolute refractory
      --S_.r_;
      S_.current_=0.0;
    }

    double gall=P_.g_L_;
    double iall=P_.i_L_;
    double tmp1,tmp2;
    for (size_t i=0; i < P_.num_of_ionchannels_; i++)
    {
        tmp1=1+S_.y4_[i];
        tmp2= V_.P40_[i]*S_.y4_[i]/(tmp1*tmp1)+P_.g_k_[i];
        gall+=tmp2;
        iall+=tmp2*P_.E_k_[i];
        S_.y4_[i] = S_.y4_[i]*V_.P44_[i];
    }
    
    double pt=std::exp(gall*V_.minus_h_Cm_);
    S_.y3_ = pt * S_.y3_ + (iall+S_.current_)/gall * (1.0-pt);

    // lower bound of membrane potential
    S_.y3_ = ( S_.y3_<P_.LowerBound_ ? P_.LowerBound_ : S_.y3_); 

    for (size_t i=0; i < P_.num_of_receptors_; i++)
    {      
      // alpha shape PSCs
      S_.y2_syn_[i] = V_.P21_syn_[i] * S_.y1_syn_[i] + V_.P22_syn_[i] * S_.y2_syn_[i];
      S_.y1_syn_[i] *= V_.P11_syn_[i];

      // collect spikes
      //S_.y1_syn_[i] += V_.PSCInitialValues_[i] * B_.spikes_[i].get_value(lag);   
      S_.y1_syn_[i] += B_.spikes_[i].get_value(lag);   
    }

    if (S_.y3_ >= P_.Theta_ && S_.r_==0)  // threshold crossing
    {
      S_.r_ = V_.RefractoryCounts_;
      //S_.y3_= P_.V_reset_; 
      for(size_t i=0; i<P_.num_of_ionchannels_; ++i)
          S_.y4_[i]=V_.Y40_[i];

      // A supra-threshold membrane potential should never be observable.
      // The reset at the time of threshold crossing enables accurate integration
      // independent of the computation step size, see [2,3] for details.

      set_spiketime(Time::step(origin.get_steps()+lag+1));
      SpikeEvent se;
      network()->send(*this, se, lag);
    }

    // set new input current
    S_.y0_ = B_.currents_.get_value(lag);

    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);
  }  
}

port mynest::glif_psc_alpha_multi::connect_sender(SpikeEvent&, port receptor_type)
{
  if (receptor_type <= 0 || receptor_type > static_cast <port>(P_.num_of_receptors_))
    throw IncompatibleReceptorType(receptor_type, get_name(), "SpikeEvent");

  P_.has_connections_ = true;
  return receptor_type;
}

void mynest::glif_psc_alpha_multi::handle(SpikeEvent& e)
{
  assert(e.get_delay() > 0);

  for (size_t i=0; i < P_.num_of_receptors_; ++i)
  {
    if (P_.receptor_types_[i] == e.get_rport()){
      B_.spikes_[i].add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			   e.get_weight() * e.get_multiplicity());
    }
  }
}

void mynest::glif_psc_alpha_multi::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t I = e.get_current();
  const double_t w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), w * I);
}

void mynest::glif_psc_alpha_multi::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

} // namespace
