/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta, Sylvain Kuppel
 *******************************************************************************/
/*
 * Infilt_GreenAmpt.cpp
 *
 *  Created on: Oct 12, 2009
 *      Author: Marco Maneta
 */

#include "Basin.h"

void Basin::Infilt_GreenAmpt(Control &ctrl, double &f, double &F, double &theta,
		double &theta2, double &theta3, double &pond, double &gw,
			     double dt, int r, int c, int step) //time step
{

  double fF, dfF = 0;
  double ef_poros1 = _porosityL1->matrix[r][c];
  double ef_poros2 = _porosityL2->matrix[r][c];
  double ef_poros3 = _porosityL3->matrix[r][c];
  double thetar = _theta_r->matrix[r][c];

  double KvKh = _KvKs->matrix[r][c];
  double Ks1 = _KsatL1->matrix[r][c] * KvKh;
  double psi = _psi_ae->matrix[r][c];

  double DT = dt; //secs
  double tp = 0; //ponding time in secs
  int k = 0;

  double fc3 = _fieldcapL3->matrix[r][c];

    if (step == 0){
      _FluxL1toL2->matrix[r][c] = 0;
      _FluxL2toL3->matrix[r][c] = 0;
    } else {
	  _FluxL1toL2_r->matrix[r][c] = 0;
      _FluxL2toL3_r->matrix[r][c] = 0;
    }
  gw = 0;

  double inp = pond / dt; //inp is potential water input in ms-1

  if (inp < RNDOFFERR){
    if (step == 0){
      _FluxSrftoL1->matrix[r][c] = 0;
    } else {
      _FluxSrftoL1_r->matrix[r][c] = 0;
	}
    return;
  }
	
  double S = (theta - thetar) / (ef_poros1 - thetar);
  double dtheta = (1 - S) * ef_poros1;

  if (dtheta < RNDOFFERR) { //if there is no room for more water
    if (step == 0){
      _FluxSrftoL1->matrix[r][c] = 0;
    } else {
      _FluxSrftoL1_r->matrix[r][c] = 0;
	}
    f = 0;
    _infilt_cap->matrix[r][c] = f;
    return;
  }

  double depth = _depth_layer1->matrix[r][c]; //_soildepth->matrix[r][c];
  double depth2 = _depth_layer2->matrix[r][c];
  double depth3 = this->_soildepth->matrix[r][c] - depth - depth2;

  double F1, Fp;
  double psidtheta = fabs(psi) * dtheta;

  double deltaF = 0;
  double i = 0;

  F = theta * depth; // / ef_poros;
  f = Ks1 * ((psidtheta / F) + 1); //infiltration capacity at time t

  if (inp <= Ks1) {
    F += inp * DT;
    theta = F / depth; //* ef_poros / depth; //cf. line 27 of this file
    pond = 0;

    // Tracking
    if (step == 0){
      _FluxSrftoL1->matrix[r][c] = inp * DT;
    } else {
      _FluxSrftoL1_r->matrix[r][c] = inp * DT;
	}

    if (theta > ef_poros1) {
      // Tracking
      //if(ctrl.sw_trck)
      if (step == 0){
	    _FluxL1toL2->matrix[r][c] = (theta - ef_poros1) * depth;
      } else {
        _FluxL1toL2_r->matrix[r][c] = (theta - ef_poros1) * depth;
	  }

      theta2 += (theta - ef_poros1) * depth / depth2;
      theta = ef_poros1;
    }
    if (theta2 > ef_poros2) {
      _FluxRecharge->matrix[r][c] += (theta2 - ef_poros2) * depth2;
      // Tracking
      //if(ctrl.sw_trck)
	  if (step == 0){
	    _FluxL2toL3->matrix[r][c] = (theta2 - ef_poros2) * depth2;
      } else {
        _FluxL2toL3_r->matrix[r][c] = (theta2 - ef_poros2) * depth2;
	  }      
      theta3 += (theta2 - ef_poros2) * depth2 / depth3;
      theta2 = ef_poros2;
    }

 
    if (theta3 > ef_poros3) {
      pond += (theta3 - ef_poros3) * depth3;
      // Subtract the "excess" infiltration" (as seen by L3)
      // it is OK because here SrftoL1>=L1toL2>=L2toL3>=L3toGW
	  if (step == 0){
        _FluxSrftoL1->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
      } else {
        _FluxSrftoL1_r->matrix[r][c]  -= (theta3 - ef_poros3) * depth3;
	  }
      _FluxRecharge->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
      // Tracking
	  if (step == 0){
        _FluxL1toL2->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
		_FluxL2toL3->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
      } else {
        _FluxL1toL2_r->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
	    _FluxL2toL3_r->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
	  }
      theta3 = ef_poros3;
    }


    gw = max<double>(0,(theta3 - fc3) * depth3);

    return;

  } else if (inp > Ks1) {
    i = min<double>(inp, f);
    tp = (Ks1 * psidtheta) / (i * (i - Ks1));

    if (tp > DT) { //if time to ponding does not happen within the time step everything infiltrates
      deltaF = i * DT;
    } else {
      //if(inp >= f) //if input intensity is larger than infiltration capacity at time 0, ponding occurs the entire time step
      //tp = 0;

      Fp = i * tp;

      //S = (theta + (Fp / depth)  - thetar)/(ef_poros - thetar); //updates relative soil moisture with moisture infiltrated before ponding
      //psidtheta = fabs(psi) * (1 - S) * ef_poros;//to recalculate infiltration after ponding using N-R with correct suction forces

      F1 = Ks1 * DT; //initial guess
      k = 0;

      do {
	deltaF = F1;
	fF = deltaF - Fp - Ks1 * (DT - tp)
	  - psidtheta * log((psidtheta + deltaF) / (psidtheta + Fp));
	dfF = deltaF / (psidtheta + deltaF);
	F1 -= fF / dfF;
	k++;
      } while (fabs(deltaF - F1) > 0.00001 && k < MAX_ITER);

      if (k >= MAX_ITER)
	cout << "WARNING: Max no iterations reached for G&A solution "
	     << endl;
      //deltaF += Fp;
      //		cout << "GA case 1 or 2 ocurring " << Fp << " " << deltaF << endl;
    }

  }
  F += deltaF;
  f = Ks1 * ((psidtheta / F) + 1);
  theta = F / depth;			//* ef_poros / depth; //cf. line 27 of this file
  pond -= deltaF;

  // Tracking
  if (step == 0){
    _FluxSrftoL1->matrix[r][c] = deltaF;
  } else {
    _FluxSrftoL1_r->matrix[r][c]  = deltaF;
  }
  _infilt_cap->matrix[r][c] = f;
  if (theta > ef_poros1) {
    theta2 += (theta - ef_poros1) * depth / depth2;
    // Tracking
    //if(ctrl.sw_trck)
    if (step == 0){
      _FluxL1toL2->matrix[r][c] = (theta - ef_poros1) * depth;
    } else {
      _FluxL1toL2_r->matrix[r][c]= (theta - ef_poros1) * depth;
	}
    theta = ef_poros1;
  }
  if (theta2 > ef_poros2) {
    _FluxRecharge->matrix[r][c] = (theta2 - ef_poros2) * depth2;
    //if(ctrl.sw_trck)
	if (step == 0){
	  _FluxL2toL3->matrix[r][c] = (theta2 - ef_poros2) * depth2;
    } else {
      _FluxL2toL3_r->matrix[r][c] = (theta2 - ef_poros2) * depth2;
	} 

    theta3 += (theta2 - ef_poros2) * depth2 / depth3;
    theta2 = ef_poros2;
  }


  if (theta3 > ef_poros3) {
    pond += (theta3 - ef_poros3) * depth3;
    // Subtract the "excess" infiltration" (as seen by L3)
    // it is OK because here SrftoL1>=L1toL2>=L2toL3
	if (step == 0){
      _FluxSrftoL1->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
    } else {
      _FluxSrftoL1_r->matrix[r][c]  -= (theta3 - ef_poros3) * depth3;
	}
    _FluxRecharge->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
    // Tracking
    //if(ctrl.sw_trck){
	if (step == 0){
      _FluxL1toL2->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
      _FluxL2toL3->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
    } else {
      _FluxL1toL2_r->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
	  _FluxL2toL3_r->matrix[r][c] -= (theta3 - ef_poros3) * depth3;
	}
      
    //}
    theta3 = ef_poros3;
  }


  gw = max<double>(0,(theta3 - fc3) * depth3);

}
