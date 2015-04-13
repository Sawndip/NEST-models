/*
 *  iaf_wsn_hermitian_2.h
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

#ifndef IAF_WSN_HERMITIAN_2
#define IAF_WSN_HERMITIAN_2

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

/* BeginDocumentation
Name: iaf_wsn_hermitian_2 - Leaky integrate-and-fire neuron model.

Description:

  iaf_wsn_hermitian_2 is an implementation of a leaky integrate-and-fire model
  with alpha-function shaped synaptic currents. Thus, synaptic currents
  and the resulting post-synaptic potentials have a finite rise time. 

  The threshold crossing is followed by an absolute refractory period
  during which the membrane potential is clamped to the resting potential.

  The linear subthresold dynamics is integrated by the Exact
  Integration scheme [1]. The neuron dynamics is solved on the time
  grid given by the computation step size. Incoming as well as emitted
  spikes are forced to that grid.  

  An additional state variable and the corresponding differential
  equation represents a piecewise constant external current.

  The general framework for the consistent formulation of systems with
  neuron like dynamics interacting by point events is described in
  [1].  A flow chart can be found in [2].

  Critical tests for the formulation of the neuron model are the
  comparisons of simulation results for different computation step
  sizes. sli/testsuite/nest contains a number of such tests.
  
  The iaf_wsn_hermitian_2 is the standard model used to check the consistency
  of the nest simulation kernel because it is at the same time complex
  enough to exhibit non-trivial dynamics and simple enough compute
  relevant measures analytically.

Remarks:

  The present implementation uses individual variables for the
  components of the state vector and the non-zero matrix elements of
  the propagator.  Because the propagator is a lower triangular matrix
  no full matrix multiplication needs to be carried out and the
  computation can be done "in place" i.e. no temporary state vector
  object is required.

  The template support of recent C++ compilers enables a more succinct
  formulation without loss of runtime performance already at minimal
  optimization levels. A future version of iaf_wsn_hermitian_2 will probably
  address the problem of efficient usage of appropriate vector and
  matrix objects.


Parameters: 

  The following parameters can be set in the status dictionary.

  V_m        double - Membrane potential in mV 
  E_L        double - Resting membrane potential in mV. 
  C_m        double - Capacity of the membrane in pF
  tau_m      double - Membrane time constant in ms.
  t_ref      double - Duration of refractory period in ms. 
  V_th       double - Spike threshold in mV.
  V_reset    double - Reset potential of the membrane in mV.
  tau_syn_ex_r double - Rise time of the excitatory synaptic alpha function in ms.
  tau_syn_ex_f double - Falling time of the excitatory synaptic alpha function in ms.
  tau_syn_in_r double - Rise time of the inhibitory synaptic alpha function in ms.
  tau_syn_in_f double - Falling time of the inhibitory synaptic alpha function in ms.
  I_e        double - Constant external input current in pA.
  V_min      double - Absolute lower value for the membrane potential.
 
Note:
  tau_m != tau_syn_{ex,in} is required by the current implementation to avoid a
  degenerate case of the ODE describing the model [1]. For very similar values,
  numerics will be unstable.

References:
  [1] Rotter S & Diesmann M (1999) Exact simulation of time-invariant linear
      systems with applications to neuronal modeling. Biologial Cybernetics
      81:381-402.
  [2] Diesmann M, Gewaltig M-O, Rotter S, & Aertsen A (2001) State space 
      analysis of synchronous spiking in cortical neural networks. 
      Neurocomputing 38-40:565-571.
  [3] Morrison A, Straube S, Plesser H E, & Diesmann M (2006) Exact subthreshold 
      integration with continuous spike times in discrete time neural network 
      simulations. Neural Computation, in press

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest
FirstVersion: September 1999
Author:  Diesmann, Gewaltig
SeeAlso: iaf_psc_delta, iaf_psc_exp, iaf_cond_exp
*/
using namespace nest;
namespace mynest
{
  class Network;

  /**
   * Leaky integrate-and-fire neuron with alpha-shaped PSCs.
   */
  class iaf_wsn_hermitian_2 : public Archiving_Node
  {
    
  public:
    
    iaf_wsn_hermitian_2();
    iaf_wsn_hermitian_2(const iaf_wsn_hermitian_2&);

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
     */

    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection&, port);
    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &); 
    
    port connect_sender(SpikeEvent&, port);
    port connect_sender(CurrentEvent&, port);
    port connect_sender(DataLoggingRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  private:

    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();

    void update(Time const &, const long_t, const long_t);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<iaf_wsn_hermitian_2>;
    friend class UniversalDataLogger<iaf_wsn_hermitian_2>;

    // ---------------------------------------------------------------- 

    struct Parameters_ {
  
      /** Membrane time constant in ms. */
      double_t Tau_; 

      /** Membrane capacitance in pF. */
      double_t C_;
    
      /** Refractory period in ms. */
      double_t TauR_;

      /** Resting potential in mV. */
      double_t U0_;

      /** External current in pA */
      double_t I_e_;

      /** Reset value of the membrane potential */
      double_t V_reset_;

      /** Threshold, RELATIVE TO RESTING POTENTIAL(!).
          I.e. the real threshold is (U0_+Theta_). */
      double_t Theta_;

      /** Rate of influence of Ie to the V amplitude **/
      double_t K_Ie_;

      /** Lower bound, RELATIVE TO RESTING POTENTIAL(!).
          I.e. the real lower bound is (LowerBound_+U0_). */
      double_t LowerBound_;

      /** all sigmas **/
      std::vector<double_t> Sigmas_;

      /** Wavelet convolution delay **/
      double_t D_Int_;

      /** Number of continuous wavelet scales **/
      size_t N_Sigmas_;

      /** Synapse number, constant **/
      //size_t num_of_receptors_;
      //std::vector<long> receptor_types_;

      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

      /** Set values from dictionary.
       * @returns Change in reversal potential E_L, to be passed to State_::set()
       */
      double set(const DictionaryDatum&);

    };
    
    // ---------------------------------------------------------------- 

    struct State_ {
      double_t u_;  //membrane voltage
      std::vector<double_t> v_;  //real-time v for all scales
      double_t vb_; //bufferred total v
      double_t Imean_; //running mean of v
      double_t Ivar_; //running variance of v
      double_t Vth_Boost_; //buffered standard variance of v
      double_t s_;  //S_enc synaptic current
      double_t Ie_; //Constant current
      //double_t Im_; //Total current
      double_t t_clk_; //last clock time
      //double_t dt_; //last clock time

      int_t    r_;  //!< Number of refractory steps remaining

      State_();  //!< Default initialization
      
      void get(DictionaryDatum&, const Parameters_&) const;

      /** Set values from dictionary.
       * @param dictionary to take data from
       * @param current parameters
       * @param Change in reversal potential E_L specified by this dict
       */
      void set(const DictionaryDatum&, const Parameters_&, double);

    };

    // ---------------------------------------------------------------- 

    struct Buffers_ {

      Buffers_(iaf_wsn_hermitian_2&);
      Buffers_(const Buffers_&, iaf_wsn_hermitian_2&);

      /** buffers and summs up incoming spikes/currents */
      RingBuffer spikes_;
      //std::vector<RingBuffer> currents_;
      RingBuffer currents_;

      //! Logger for all analog data
      UniversalDataLogger<iaf_wsn_hermitian_2> logger_;

    };
    
    // ---------------------------------------------------------------- 

    struct Variables_ {

      /** Amplitude of the synaptic current.
	  This value is chosen such that a post-synaptic potential with
	  weight one has an amplitude of 1 mV.
       */
      int_t    RefractoryCounts_;
    
      std::vector<double_t> P2_;
      double_t P32_;
      double_t P33_;

    };

    //ODE
    double_t get_Im_(const double_t, const size_t);

    // Access functions for UniversalDataLogger -------------------------------

    //! Read out the real membrane potential
    double_t get_V_m_() const { return S_.u_;}
    double_t get_V_() const {return S_.v_[0];}
    double_t get_VB_() const {return S_.vb_;}
    double_t get_Syn_() const {return S_.s_;}
    double_t get_Ie_() const {return S_.Ie_;}
    double_t get_Vth_Boost_() const {return S_.Vth_Boost_;}
    //double_t get_I_Var_() const {return S_.Ivar_;}
    //double_t get_I_Mean_() const {return S_.Imean_;}
    //double_t get_dt_() const {return S_.dt_;}



    // Data members ----------------------------------------------------------- 
    
    /**
     * @defgroup iaf_wsn_hermitian_2_data
     * Instances of private data structures for the different types
     * of data pertaining to the model.
     * @note The order of definitions is important for speed.
     * @{
     */   
    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;
    /** @} */
    
    //! Mapping of recordables names to access functions
    static RecordablesMap<iaf_wsn_hermitian_2> recordablesMap_;
  };

  inline
  port iaf_wsn_hermitian_2::check_connection(Connection& c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }
    
  inline
  port iaf_wsn_hermitian_2::connect_sender(SpikeEvent&, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
   
  inline
  port iaf_wsn_hermitian_2::connect_sender(CurrentEvent&, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
   
  inline
  port iaf_wsn_hermitian_2::connect_sender(DataLoggingRequest& dlr, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }
  
  inline
  void iaf_wsn_hermitian_2::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d, P_);
    Archiving_Node::get_status(d);
  
    (*d)[names::recordables] = recordablesMap_.get_list();
  }
  
  inline
  void iaf_wsn_hermitian_2::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;            // temporary copy in case of errors
    const double delta_EL = ptmp.set(d);         // throws if BadProperty
    State_      stmp = S_;            // temporary copy in case of errors
    stmp.set(d, ptmp, delta_EL);                 // throws if BadProperty
  
    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    Archiving_Node::set_status(d);
  
    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }

} // namespace

#endif /* #ifndef IAF_PSC_ALPHA_H_EXT */
