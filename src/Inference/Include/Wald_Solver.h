#ifndef __WALD_SOLVER_H__
#define __WALD_SOLVER_H__

// HEADERS
#include "../../FdaPDE.h"
#include "../../Regression/Include/Mixed_FE_Regression.h"
#include "../../Regression/Include/Regression_Data.h"
#include "../../Lambda_Optimization/Include/Optimization_Data.h"
#include "../../Lambda_Optimization/Include/Solution_Builders.h"
#include "Inference_Data.h"
#include "Inference_Carrier.h"
#include "Inverter.h" // To be implemented

// *** Wald_Solver Class ***
//! Hypothesis testing and confidence intervals using Wald implementation
/*!
 This class performes hypothesis testing and/or computes confidence intervals using a Wald-type approach. It contains a reference to an inverter, that manages to compute the invertion of matrixNoCov in an exact or non-exact way; It contains a reference to an Inference_Carrier object that wraps all the information needed to make inference. There is only one public method that calls the proper private methods to compute what is requested by the user.
*/
template<typename InputHandler>
class Wald_Solver{
 private:
  // inverter object that computes the inverse of matrixNoCov in exact/non-exact way.
  Inverse_Base & inverter;
  // inference carrier that contains all the information needed for inference. 
  const InferenceCarrier<InputHandler> & inf_car;
  // smoothing matrix. 
  MatrixXr S;
  // transpose of the smoothing matrix.
  MatrixXr S_t;   
    // boolean to know whether S has been computed
  bool is_S_computed = false;
  // Variance-Covariance matrix of the beta parameters
  MatrixXr V;
  // boolean to know whether V has been computed
  bool is_V_computed = false;
  // method used to compute S
  void compute_S(void);
  // method used to compute V
  void compute_V(void);
  // method to compute the pvalue of the test (they can be more than one if one-at-the-time)
  VectorXr compute_pvalue(void);
  // method to compute the confidence intervals (they can be simultaneous, one-at-the-time or bonferroni)
  MatrixXv compute_CI(void); // To be implemented

  public:
  // CONSTUCTOR
  Wald_Solver()=delete;	//The default constructor is deleted
  Wald_Solver(Inverse_Base & inverter_, const Inference_Carrier<InputHandler> & inf_car_):inverter(inverter_), inf_car(inf_car_){}; 

  // GETTERS
  inline const MatrixXr * getSp (void) const {return &this->S;}      //!< Getter of Sp \return Sp
  inline const MatrixXr * getS_tp (void) const {return &this->S_t;}  //!< Getter of Sp_tp \return Sp_tp
  inline const MatrixXr * getVp (void) const {return &this->V;}      //!< Getter of Vp \ return Vp

  // public method that calls the requested functions according to test_type and interval_type
  MatrixXv compute_inference_output (void);
	
};



#endif
