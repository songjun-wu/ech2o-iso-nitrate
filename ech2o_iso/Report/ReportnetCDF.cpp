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
 *    Marco Maneta
 *******************************************************************************/
/*
 * ReportnetCDF.cpp
 *
 *  Created on: Feb 10, 2021
 *      Author: Xiaoqiang Yang
 *      This script includes two functions creating (ts = 1) and updating the nc files
 *      for water fluxes (Water_Fluxes_States.nc), deutrium tracer (D_Fluxes_States.nc),
 *      oxygen-18 tracer (O18_Fluxes_States.nc), water age (Age_Fluxes_States.nc) and
 *      vegtation dynamics (VegDyn_Fluxes_States.nc).
 *  Technicial notes:
 *      netCDF C++4 library should be installed (e.g., in UFZ eve-cluster: module load netCDF-C++4) 
 *      and configured in the compiling folder (e.g., add "-lnetcdf_c++4" in ./Release-Linux/objects.mk)
 */
#include <fstream>
#include <netcdf>
#include "Report.h"
#include "Basin.h"
#include "Sativa.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;



// Names of things. 
#define LAT_NAME "latitude"
#define LON_NAME "longitude"
#define REC_NAME "time"

string  UNITS = "units";
string  DEGREES_EAST =  "degrees_east";
string  DEGREES_NORTH = "degrees_north";

// For the units attributes. 
string LAT_UNITS = "degrees_north";
string LON_UNITS = "degrees_east";

// Return this code to the OS in case of failure.
#define NC_ERR 2

int Report::CreatOutputNC(string filepath, string outTP){
   
   REAL8 Rsol;
   UINT4  NLAT, NLON;
   REAL8 START_LAT,START_LON;
   REAL8 NODATA;
   std::vector<std::string> RprtNamesW,RprtNamesTD,RprtNamesTO,RprtNamesTA,RprtNamesVG;
   std::vector<std::string> outvars;
   string filename;
   
  //checking output variables
    if (outTP == "W"){
	  if(oControl->Rep_Long_Rad_Down)
		RprtNamesW.push_back("Ldown");
	  if(oControl->Rep_Short_Rad_Down)
        RprtNamesW.push_back("Sdown");
	  if(oControl->Rep_Precip)
        RprtNamesW.push_back("Pp");
	  if(oControl->Rep_Rel_Humidity)
		RprtNamesW.push_back("RH");
	  if(oControl->Rep_Wind_Speed)
        RprtNamesW.push_back("WndSp");
	  if(oControl->Rep_AvgAir_Temperature)
		RprtNamesW.push_back("Tp");
	  if(oControl->Rep_MinAir_Temperature)
		RprtNamesW.push_back("TpMin");
	  if(oControl->Rep_MaxAir_Temperature)
		RprtNamesW.push_back("TpMax");

	  if (oControl->Rep_Canopy_Water_Stor_sum)
		RprtNamesW.push_back("Cs");
	  if (oControl->Rep_SWE)
		RprtNamesW.push_back("SWE");
	  if (oControl->Rep_Infilt_Cap)
		RprtNamesW.push_back("IfCap");
	  if (oControl->Rep_Streamflow)
		RprtNamesW.push_back("Q");
	  if (oControl->Rep_Saturation_Area)
		RprtNamesW.push_back("SatArea_");
	  if (oControl->Rep_Ponding)
		RprtNamesW.push_back("Ponding_");
	  if (oControl->Rep_Soil_Water_Content_Average)
		RprtNamesW.push_back("SWCav");
	  if (oControl->Rep_Soil_Water_Content_12)
		RprtNamesW.push_back("SWCup");
	  if (oControl->Rep_Soil_Water_Content_L1)
		RprtNamesW.push_back("SWC1_");
	  if (oControl->Rep_Soil_Water_Content_L2)
		RprtNamesW.push_back("SWC2_");
	  if (oControl->Rep_Soil_Water_Content_L3)
		RprtNamesW.push_back("SWC3_");
	  if (oControl->Rep_WaterTableDepth)
		RprtNamesW.push_back("WTD_");
	  if (oControl->Rep_Field_Capacity_L1)
		RprtNamesW.push_back("FCap1");
	  if (oControl->Rep_Field_Capacity_L2)
		RprtNamesW.push_back("FCap2");
	  if (oControl->Rep_Field_Capacity_L3)
		RprtNamesW.push_back("FCap3");
	  if (oControl->Rep_Soil_Sat_Deficit)
		RprtNamesW.push_back("SatDef");
	  if (oControl->Rep_GWater)
		RprtNamesW.push_back("GW");
	  if (oControl->Rep_Soil_Net_Rad)
		RprtNamesW.push_back("NRs");
	  if (oControl->Rep_Net_Rad_sum)
		RprtNamesW.push_back("NRtot");
	  if (oControl->Rep_Soil_LE)
		RprtNamesW.push_back("LEs");
	  if (oControl->Rep_Sens_Heat)
		RprtNamesW.push_back("SensH");
	  if (oControl->Rep_Grnd_Heat)
		RprtNamesW.push_back("GrndH");
	  if (oControl->Rep_Snow_Heat)
		RprtNamesW.push_back("SnowH");
	  if (oControl->Rep_Soil_Temperature)
		RprtNamesW.push_back("Ts");
	  if (oControl->Rep_Skin_Temperature)
		RprtNamesW.push_back("Tskin");
	  if (oControl->Rep_Skin_Temperature)
		RprtNamesW.push_back("Twater");
	  if (oControl->Rep_Skin_Temperature)
		RprtNamesW.push_back("EvapC");
	  if (oControl->Rep_Total_ET)
		RprtNamesW.push_back("Evap");
	  if (oControl->Rep_Transpiration_sum)
		RprtNamesW.push_back("EvapT");
	  if (oControl->Rep_Einterception_sum)
		RprtNamesW.push_back("EvapI");
	  if (oControl->Rep_Esoil_sum)
		RprtNamesW.push_back("EvapS");

	  // Some vertical- lateral fluxes
	  if (oControl->Rep_GWtoChn)
		RprtNamesW.push_back("GWChn");
	  if (oControl->Rep_SrftoChn)
		RprtNamesW.push_back("SrfChn");
	  if (oControl->Rep_Infilt)
		RprtNamesW.push_back("Inf");
	  if (oControl->Rep_Exfilt)
		RprtNamesW.push_back("RSrf");
	  if (oControl->Rep_Recharge)
		RprtNamesW.push_back("Rchg");
	  if (oControl->Rep_LattoSrf)
		RprtNamesW.push_back("LSrfi");
	  if (oControl->Rep_LattoChn)
		RprtNamesW.push_back("LChni");
	  if (oControl->Rep_LattoGW)
		RprtNamesW.push_back("LGWi");
	  if (oControl->Rep_SrftoLat)
		RprtNamesW.push_back("LSrfo");
	  if (oControl->Rep_GWtoLat)
		RprtNamesW.push_back("LGWo");
	  if (oControl->Rep_BedRock_leakage)
		RprtNamesW.push_back("LeakF");	

	  if(oControl->sw_trck and oControl->sw_TPD){
		if (oControl->Rep_Moist_MW1)
		  RprtNamesW.push_back("SWCMW1");
		if (oControl->Rep_Moist_MW2)
		  RprtNamesW.push_back("SWCMW2");
		if (oControl->Rep_Frac_MW1)
		  RprtNamesW.push_back("fMW1");
		if (oControl->Rep_Frac_MW2)
		  RprtNamesW.push_back("fMW2");
		if (oControl->Rep_Frac_MW12)
		  RprtNamesW.push_back("fMW12");
	  }
	  //extra GW yangx 2020-05
	  if(oControl->sw_extraGW){
		if (oControl->Rep_ExtraGWater)
		  RprtNamesW.push_back("ExtraGW");
		if (oControl->Rep_ExtraGWtoChn)
		  RprtNamesW.push_back("ExGWtC");
		if (oControl->Rep_ExtraGWtoChnacc)
		  RprtNamesW.push_back("ExGWtCA");
		if (oControl->Rep_LattoExtraGW)
		  RprtNamesW.push_back("LtExGW");
		if (oControl->Rep_LattoExtraGWacc)
		  RprtNamesW.push_back("LtExGWA");
		if (oControl->Rep_ExtraGWtoLat)
		  RprtNamesW.push_back("ExGWtL");
		if (oControl->Rep_ExtraGWtoLatacc)
		  RprtNamesW.push_back("ExGWtLA");
	  }
    }
	  // Tracking maps
	if(outTP == "TD" && oControl->sw_trck && oControl->sw_2H){
		if (oControl->Rep_d2Hcanopy_sum)
		  RprtNamesTD.push_back("dHcnp");
		if (oControl->Rep_d2Hprecip)
		  RprtNamesTD.push_back("dHpcp");
		if (oControl->Rep_d2Hsnowpack)
		  RprtNamesTD.push_back("dHsnw");
		if (oControl->Rep_d2Hsurface)
		  RprtNamesTD.push_back("dHsrf");
		if (oControl->Rep_d2HSrftoChn)
		  RprtNamesTD.push_back("dHSrfQ");
		if (oControl->Rep_d2Hchan)
		  RprtNamesTD.push_back("dHchn");
		if (oControl->Rep_d2Hsoil1)
		  RprtNamesTD.push_back("dHsL1");
		if (oControl->Rep_d2Hsoil2)
		  RprtNamesTD.push_back("dHsL2");
		if (oControl->Rep_d2HsoilUp)
		  RprtNamesTD.push_back("dHsUp");
		if (oControl->Rep_d2Hsoil3)
		  RprtNamesTD.push_back("dHsL3");
		if (oControl->Rep_d2HsoilAv)
		  RprtNamesTD.push_back("dHsAv");
		if (oControl->Rep_d2Hgroundwater)
		  RprtNamesTD.push_back("dHgw");
		if (oControl->Rep_d2HevapS_sum)
		  RprtNamesTD.push_back("dHeS");
		if (oControl->Rep_d2HevapI_sum)
		  RprtNamesTD.push_back("dHeI");
		if (oControl->Rep_d2HevapT_sum)
		  RprtNamesTD.push_back("dHeT");
		if (oControl->Rep_d2Hleakage)
		  RprtNamesTD.push_back("dHLeak");		  
		if(oControl->sw_TPD){
		  if (oControl->Rep_d2H_MW1)
		    RprtNamesTD.push_back("dHMW1");
		  if (oControl->Rep_d2H_MW2)
		    RprtNamesTD.push_back("dHMW2");
		  if (oControl->Rep_d2H_TB1)
		    RprtNamesTD.push_back("dHTB1");
		  if (oControl->Rep_d2H_TB2)
		    RprtNamesTD.push_back("dHTB2");
		}
		//extra GW yangx 2020-05
		if(oControl->sw_extraGW){
		  if (oControl->Rep_d2HExtraGWtoLat)
		    RprtNamesTD.push_back("dHExGWO");
		}
	}
	  
    if(outTP == "TO" && oControl->sw_trck && oControl->sw_18O){
		if (oControl->Rep_d18Ocanopy_sum)
		  RprtNamesTO.push_back("dOcnp");
		if (oControl->Rep_d18Oprecip)
		  RprtNamesTO.push_back("dOpcp");
		if (oControl->Rep_d18Osnowpack)
		  RprtNamesTO.push_back("dOsnw");
		if (oControl->Rep_d18Osurface)
		  RprtNamesTO.push_back("dOsrf");
		if (oControl->Rep_d18OSrftoChn)
		  RprtNamesTO.push_back("dOSrfQ");
		if (oControl->Rep_d18Ochan)
		  RprtNamesTO.push_back("dOchn");
		if (oControl->Rep_d18Osoil1)
		  RprtNamesTO.push_back("dOsL1");
		if (oControl->Rep_d18Osoil2)
		  RprtNamesTO.push_back("dOsL2");
		if (oControl->Rep_d18OsoilUp)
		  RprtNamesTO.push_back("dOsUp");
		if (oControl->Rep_d18Osoil3)
		  RprtNamesTO.push_back("dOsL3");
		if (oControl->Rep_d18OsoilAv)
		  RprtNamesTO.push_back("dOsAv");
		if (oControl->Rep_d18Ogroundwater)
		  RprtNamesTO.push_back("dOgw");
		if (oControl->Rep_d18OevapS_sum)
		  RprtNamesTO.push_back("dOeS");
		if (oControl->Rep_d18OevapI_sum)
		  RprtNamesTO.push_back("dOeI");
		if (oControl->Rep_d18OevapT_sum)
		  RprtNamesTO.push_back("dOeT");
		if(oControl->sw_TPD){
		  if (oControl->Rep_d18O_MW1)
		RprtNamesTO.push_back("dOMW1");
		  if (oControl->Rep_d18O_MW2)
		RprtNamesTO.push_back("dOMW2");
		  if (oControl->Rep_d18O_TB1)
		RprtNamesTO.push_back("dOTB1");
		  if (oControl->Rep_d18O_TB2)
		RprtNamesTO.push_back("dOTB2");
		}
		//extra GW yangx 2020-05
		if(oControl->sw_extraGW){
		  if (oControl->Rep_d18Oleakage)
		    RprtNamesTO.push_back("dOLeak");
		  if (oControl->Rep_d18OExtraGWtoLat)
		    RprtNamesTO.push_back("dOExGWO");
		}
	}
 
    if(outTP == "TA" && oControl->sw_trck && oControl->sw_Age){
		if (oControl->Rep_Agecanopy_sum)
		 RprtNamesTA.push_back("Agecnp");
		if (oControl->Rep_Agesnowpack)
		  RprtNamesTA.push_back("Agesnw");
		if (oControl->Rep_Agesurface)
		  RprtNamesTA.push_back("Agesrf");
		if (oControl->Rep_Agechan)
		  RprtNamesTA.push_back("Agechn");
		if (oControl->Rep_Agesoil1)
		  RprtNamesTA.push_back("AgesL1");
		if (oControl->Rep_Agesoil2)
		  RprtNamesTA.push_back("AgesL2");
		if (oControl->Rep_AgesoilUp)
		  RprtNamesTA.push_back("AgesUp");
		if (oControl->Rep_Agesoil3)
		  RprtNamesTA.push_back("AgesL3");
		if (oControl->Rep_AgesoilAv)
		  RprtNamesTA.push_back("AgesAv");
		if (oControl->Rep_Agegroundwater)
		  RprtNamesTA.push_back("Agegw");
		if (oControl->Rep_AgeevapS_sum)
		  RprtNamesTA.push_back("AgeeS");
		if (oControl->Rep_AgeevapI_sum)
		  RprtNamesTA.push_back("AgeeI");
		if (oControl->Rep_AgeevapT_sum)
		  RprtNamesTA.push_back("AgeeT");
		if (oControl->Rep_AgeGWtoChn)
		  RprtNamesTA.push_back("AgeGWQ");
		if (oControl->Rep_AgeSrftoChn)
		  RprtNamesTA.push_back("AgeSrfQ");
		if (oControl->Rep_AgeRecharge)
		  RprtNamesTA.push_back("AgeR");

		if(oControl->sw_TPD){
		  if (oControl->Rep_Age_MW1)
		RprtNamesTA.push_back("AgeMW1");
		  if (oControl->Rep_Age_MW2)
		RprtNamesTA.push_back("AgeMW2");
		  if (oControl->Rep_Age_MWUp)
		RprtNamesTA.push_back("AgeMW12");
		  if (oControl->Rep_Age_TB1)
		RprtNamesTA.push_back("AgeTB1");
		  if (oControl->Rep_Age_TB2)
		RprtNamesTA.push_back("AgeTB2");
		  if (oControl->Rep_Age_TBUp)
		RprtNamesTA.push_back("AgeTB12");
		}
		//extra GW yangx 2020-05
		if(oControl->sw_extraGW){
		  if (oControl->Rep_Ageleakage)
		RprtNamesTA.push_back("AgLeak");
		  if (oControl->Rep_AgeExtraGWtoLat)
		RprtNamesTA.push_back("AgExGWO");
		  if (oControl->Rep_AgeExtraGWtoChn)
		RprtNamesTA.push_back("AgExOC");
		}
    }
    if (outTP == "VG"){  
	  for (int i = 0; i < oControl->NumSpecs; i++) {
		stringstream name;

		if (oControl->Rep_Veget_frac) {
		  name << "p" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Stem_Density) {
		  name << "ntr" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_RootFrac1Species) {
		  name << "L1Ro" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}
		if (oControl->Rep_RootFrac2Species) {
		  name << "L2Ro" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Leaf_Area_Index) {
		  name << "lai" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Stand_Age) {
		  name << "age" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Canopy_Conductance) {
		  name << "gc" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_GPP) {
		  name << "gpp" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_NPP) {
		  name << "npp" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Basal_Area) {
		  name << "bas" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Tree_Height) {
		  name << "hgt" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Root_Mass) {
		  name << "root" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Canopy_Temp) {
		  name << "Tc" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Canopy_NetR) {
		  name << "NRc" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Canopy_LE_E) {
		  name << "LEEi" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}
		if (oControl->Rep_Canopy_LE_T) {
		  name << "LETr" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Canopy_Sens_Heat) {
		  name << "Hc" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_Canopy_Water_Stor) {
		  name << "Cs" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		if (oControl->Rep_ETspecies) {
		  name << "ETc" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}
		if (oControl->Rep_Transpiration) {
		  name << "Tr" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}
		if (oControl->Rep_Einterception) {
		  name << "Ei" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}
		if (oControl->Rep_Esoil) {
		  name << "Es" << i << "_";
		  RprtNamesVG.push_back(name.str());
		  name.str("");
		}

		// Tracking
		if(oControl->sw_trck){
		  if (oControl->sw_2H && oControl->Rep_d2Hcanopy) {
		name << "dHCs" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }
		  
		  if (oControl->sw_18O && oControl->Rep_d18Ocanopy) {
		name << "dOCs" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }
		  
		  if (oControl->sw_Age && oControl->Rep_Agecanopy) {
		name << "AgeCs" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }
		  if (oControl->sw_2H && oControl->Rep_d2HevapI) {
		name << "dHeI" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }

		  if (oControl->sw_18O && oControl->Rep_d18OevapI) {
		name << "dOeI" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }

		  if (oControl->sw_Age && oControl->Rep_AgeevapI) {
		name << "AgeeI" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }
		  if (oControl->sw_2H && oControl->Rep_d2HevapT) {
		name << "dHeT" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }

		  if (oControl->sw_18O && oControl->Rep_d18OevapT) {
		name << "dOeT" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }

		  if (oControl->sw_Age && oControl->Rep_AgeevapT) {
		name << "AgeeT" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }
		  if (oControl->sw_2H && oControl->Rep_d2HevapS) {
		name << "dHeS" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }

		  if (oControl->sw_18O && oControl->Rep_d18OevapS) {
		name << "dOeS" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }

		  if (oControl->sw_Age && oControl->Rep_AgeevapS) {
		name << "AgeeS" << i << "_";
		RprtNamesVG.push_back(name.str());
		name.str("");
		  }
		}
      }
	}
	//CREATE OUTPUTFILES
    try
    {      
      //only create nc files with output variables
      if(outTP == "W") {  
	    if (!RprtNamesW.empty()){
		  // Create the file.
          filename =  filepath + "Water_Fluxes_States.nc";
		  outvars = RprtNamesW;
        }
	  }else if(outTP == "TD"){
	    if (!RprtNamesTD.empty()){
	      // Create the file.
		  filename =  filepath + "D_Fluxes_States.nc";
		  outvars = RprtNamesTD;
        }
	  }else if(outTP == "TO"){
	    if (!RprtNamesTO.empty()){
	     // Create the file.
          filename =  filepath + "O18_Fluxes_States.nc";
		  outvars = RprtNamesTO;
        }
	  }else if(outTP == "TA"){
	    if (!RprtNamesTA.empty()){
	      // Create the file.
          filename =  filepath + "Age_Fluxes_States.nc";
          outvars = RprtNamesTA;
	    }
	  }else if(outTP == "VG"){
	    if (!RprtNamesVG.empty()){
	      // Create the file.
          filename =  filepath + "VegDyn_Fluxes_States.nc";
		  outvars = RprtNamesVG;
        }	
      }

	  if (!outvars.empty()){
		//get basic information
		Rsol = oBasin->getCellSize(); //m
		NLAT = oBasin->getNumRows();
		NLON = oBasin->getNumCols();
		START_LAT = oBasin->getSouthCoord();
		START_LON = oBasin->getWestCoord();
		NODATA = oBasin->getModelnodata();
		  
		//declare variable to store data
		REAL8 lats[NLAT],lons[NLON];
		REAL8 varout[NLAT][NLON];

		//prepare the lat and lon values at the central of each grid
		for (UINT4 lat = 0; lat < NLAT; lat++)
	      lats[lat] = START_LAT + Rsol * 0.5 + Rsol * lat;
		for (UINT4 lon = 0; lon < NLON; lon++)
		  lons[lon] = START_LON + Rsol * 0.5 + Rsol * lon;
		//put initial NLL values
		for(unsigned i = 0; i < NLAT; i++){
          for(unsigned j = 0; j < NLON; j++){
            varout[i][j] = NODATA;
	      }
        }
		//create the file		   
	  	NcFile file0(filename, NcFile::replace);
		  
		// Define the dimensions. NetCDF will hand back an ncDim object for
		// each.
		NcDim latDim = file0.addDim(LAT_NAME, NLAT);
		NcDim lonDim = file0.addDim(LON_NAME, NLON);
		NcDim recDim = file0.addDim(REC_NAME);  //adds an unlimited dimension
		   
		  // Define the coordinate variables.
		NcVar latVar = file0.addVar(LAT_NAME, ncDouble, latDim);
		NcVar lonVar = file0.addVar(LON_NAME, ncDouble, lonDim);
		   
		  // Define units attributes for coordinate vars. This attaches a
		  // text attribute to each of the coordinate variables, containing
		  // the units.
		latVar.putAtt(UNITS, DEGREES_NORTH);
		lonVar.putAtt(UNITS, DEGREES_EAST);
		  // Write the coordinate variable data to the file.
		latVar.putVar(lats);
		lonVar.putVar(lons);
		  
		  // Define the netCDF variables should be reported
		  // data.
		vector<NcDim> dimVector;
		dimVector.push_back(recDim);
		dimVector.push_back(latDim);
		dimVector.push_back(lonDim);

		vector<size_t> startp,countp;
		startp.push_back(0);
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(1);
		countp.push_back(NLAT);
		countp.push_back(NLON);
		//assign initial record with nodata value
		for (auto item : outvars) { 
	      NcVar rec = file0.addVar(item.c_str(), ncDouble, dimVector);
		  rec.putVar(startp,countp,varout);
		}
		  
		// Global attributes
		file0.putAtt("title", "Map of modelled information");
		file0.putAtt("institution", "Landscape ecohydrology, Leibniz - IGB, Berlin");
		file0.putAtt("source", "Model outputs of the EcH2O_iso model");
        file0.putAtt("history", "The same variables as the original PCRmaps");
		file0.putAtt("comment", "The first record of each variable is initialized as NULL");

		//created nc file
		cout << "*** SUCCESS Creating file " << filename << "!" << endl;
	  }else{
	      cout << "*** No state variable or fluxe output should be wirtten as nc file!" << endl;
      }
      // The file is automatically closed by the destructor. This frees
      // up any internal netCDF resources associated with the file, and
      // flushes any buffers.
      return 0;
    }
    catch(NcException& e)
    {
      e.what(); 
      return NC_ERR;
    }
}


int Report::CreatALLOutputNC(string filepath, string outTP){
   
   REAL8 Rsol;
   UINT4  NLAT, NLON;
   REAL8 START_LAT,START_LON;
   REAL8 NODATA;
   std::vector<std::string> RprtNamesW,RprtNamesTD,RprtNamesTO,RprtNamesTA,RprtNamesVG;
   std::vector<std::string> outvars;
   string filename;
   
  //checking output variables
        RprtNamesW.push_back("Pp");
		RprtNamesW.push_back("Tp");
		RprtNamesW.push_back("Cs");
		RprtNamesW.push_back("SWE");
		RprtNamesW.push_back("Ponding_");
		RprtNamesW.push_back("EvapI");		
		RprtNamesW.push_back("ThrghRn");//NEW	
		RprtNamesW.push_back("ThrghSn");//NEW		
		RprtNamesW.push_back("SnMelt");//NEW
		RprtNamesW.push_back("Inf");
		RprtNamesW.push_back("L1toL2");//NEW
		RprtNamesW.push_back("L2toL3");//NEW
	 //reinfiltration	
        if (oControl->sw_reinfilt){
		  RprtNamesW.push_back("InfR");//NEW
		  RprtNamesW.push_back("L1toL2R");//NEW
		  RprtNamesW.push_back("L2toL3R");//NEW		
		}
		RprtNamesW.push_back("RSrf");
		RprtNamesW.push_back("L2toL1");
		RprtNamesW.push_back("L3toL2");
		RprtNamesW.push_back("EvapS");
		RprtNamesW.push_back("EvapT1");//NEW	
		RprtNamesW.push_back("EvapT2");//NEW	
		RprtNamesW.push_back("EvapT3");//NEW	
		RprtNamesW.push_back("SWD1_");//NEW
		RprtNamesW.push_back("SWD2_");//NEW
		RprtNamesW.push_back("SWD3_");//NEW
		RprtNamesW.push_back("satSM1");//NEW
		RprtNamesW.push_back("satSM2");//NEW
		RprtNamesW.push_back("satSM3");//NEW
		RprtNamesW.push_back("Wiltpnt");//NEW
		RprtNamesW.push_back("Ts");
		RprtNamesW.push_back("GW");
		RprtNamesW.push_back("LSrfo");
		RprtNamesW.push_back("LSrfi");
		RprtNamesW.push_back("LGWo");
		RprtNamesW.push_back("LGWi");
		RprtNamesW.push_back("SrfChn");
		if (oControl->toggle_SrfOVFmix == 1){
		  RprtNamesW.push_back("RSrfmix");  //NEW
		}
		RprtNamesW.push_back("GWChn");
 		RprtNamesW.push_back("ChnS"); //NEW       
		RprtNamesW.push_back("LChQUi");//NEW
		RprtNamesW.push_back("LChQDo");//NEW
		RprtNamesW.push_back("Q");
		RprtNamesW.push_back("LeakD");	
	    //extra GW yangx 2020-05
	    if(oControl->sw_extraGW){
		  RprtNamesW.push_back("ExtraGW");
		  RprtNamesW.push_back("ExGWtC");
		  RprtNamesW.push_back("LtExGW");
		  RprtNamesW.push_back("ExGWtL");
	    }

	//CREATE OUTPUTFILES
    try
    {      
      //only create nc files with output variables 
	  // Create the file.
      filename =  filepath + "Water_Fluxes_States.nc";
	  outvars = RprtNamesW;

		//get basic information
		Rsol = oBasin->getCellSize(); //m
		NLAT = oBasin->getNumRows();
		NLON = oBasin->getNumCols();
		START_LAT = oBasin->getSouthCoord();
		START_LON = oBasin->getWestCoord();
		NODATA = oBasin->getModelnodata();
		  
		//declare variable to store data
		REAL8 lats[NLAT],lons[NLON];
		REAL8 varout[NLAT][NLON];

		//prepare the lat and lon values at the central of each grid
		for (UINT4 lat = 0; lat < NLAT; lat++)
	      lats[lat] = START_LAT + Rsol * 0.5 + Rsol * lat;
		for (UINT4 lon = 0; lon < NLON; lon++)
		  lons[lon] = START_LON + Rsol * 0.5 + Rsol * lon;
		//put initial NLL values
		for(unsigned i = 0; i < NLAT; i++){
          for(unsigned j = 0; j < NLON; j++){
            varout[i][j] = NODATA;
	      }
        }
		//create the file		   
	  	NcFile file0(filename, NcFile::replace);
		  
		// Define the dimensions. NetCDF will hand back an ncDim object for
		// each.
		NcDim latDim = file0.addDim(LAT_NAME, NLAT);
		NcDim lonDim = file0.addDim(LON_NAME, NLON);
		NcDim recDim = file0.addDim(REC_NAME);  //adds an unlimited dimension
		   
		  // Define the coordinate variables.
		NcVar latVar = file0.addVar(LAT_NAME, ncDouble, latDim);
		NcVar lonVar = file0.addVar(LON_NAME, ncDouble, lonDim);
		   
		  // Define units attributes for coordinate vars. This attaches a
		  // text attribute to each of the coordinate variables, containing
		  // the units.
		latVar.putAtt(UNITS, DEGREES_NORTH);
		lonVar.putAtt(UNITS, DEGREES_EAST);
		  // Write the coordinate variable data to the file.
		latVar.putVar(lats);
		lonVar.putVar(lons);
		  
		  // Define the netCDF variables should be reported
		  // data.
		vector<NcDim> dimVector;
		dimVector.push_back(recDim);
		dimVector.push_back(latDim);
		dimVector.push_back(lonDim);

		vector<size_t> startp,countp;
		startp.push_back(0);
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(1);
		countp.push_back(NLAT);
		countp.push_back(NLON);
		//assign initial record with nodata value
		for (auto item : outvars) { 
	      NcVar rec = file0.addVar(item.c_str(), ncDouble, dimVector);
		  rec.putVar(startp,countp,varout);
		}
		  
		// Global attributes
		file0.putAtt("title", "Map of modelled information");
		file0.putAtt("institution", "Landscape ecohydrology, Leibniz - IGB, Berlin");
		file0.putAtt("source", "Model outputs of the EcH2O_iso model");
        file0.putAtt("history", "The same variables as the original PCRmaps");
		file0.putAtt("comment", "The first record of each variable is initialized as NULL");

		//created nc file
		cout << "*** SUCCESS Creating file " << filename << "!" << endl;
      // The file is automatically closed by the destructor. This frees
      // up any internal netCDF resources associated with the file, and
      // flushes any buffers.
      return 0;
    }
    catch(NcException& e)
    {
      e.what(); 
      return NC_ERR;
    }
}

//************************************************************************

int Report::UpdateOutputNC(const grid *input, string varname, string outTP){

   UINT4 NLAT, NLON;
   int TScount;
   string filename;
   //from Basin
   NLAT = oBasin->getNumRows();
   NLON = oBasin->getNumCols();  
   TScount = oControl->current_ts_count;
   REAL8 outdata[NLAT][NLON];   
   //get the outdata from input map
   for(unsigned i = 0; i < NLAT; i++){
     for(unsigned j = 0; j < NLON; j++){
       outdata[i][j] = input->matrix[i][j];
	 }
   }   
   //indexing the record location
   //c++ starts with index 0 AND the first nodata record should be replaced so
   TScount -= 1; 
   
   try
   {
    //open the right nc file according to the output type "outTP"
	if(outTP == "W") {
      filename =  oControl->path_ResultsFolder + "Water_Fluxes_States.nc";
	}else if(outTP == "TD"){
      filename =  oControl->path_ResultsFolder + "D_Fluxes_States.nc";
	}else if(outTP == "TO"){
      filename =  oControl->path_ResultsFolder + "O18_Fluxes_States.nc";	
	}else if(outTP == "TA"){
      filename =  oControl->path_ResultsFolder + "Age_Fluxes_States.nc";	
	}else if(outTP == "VG"){
      filename =  oControl->path_ResultsFolder + "VegDyn_Fluxes_States.nc";
	}
	//open the file
	NcFile file1(filename, NcFile::write);
	if(file1.isNull()){
		cout << "ERROR*** file " << filename << " doesn't exist!" << endl;
		return NC_ERR;
	}
/* 	// Get the  latitude and longitude coordinate variables 
	   NcVar latVar, lonVar;
	   latVar = dataFile.getVar("latitude");
	   if(latVar.isNull()) return NC_ERR;
	   lonVar = dataFile.getVar("longitude");
	   if(lonVar.isNull()) return NC_ERR; */  
	//Get the record variable
 	NcVar recVar;
    recVar = file1.getVar(varname.c_str());
    if(recVar.isNull()){
	  cout << "ERROR*** variable " << varname.c_str() << " doesn't exist!" << endl;
	  return NC_ERR;
    }

    //Get using the current_ts_count for the start index of new record  
    //cout << varname.c_str() << "existing number of records: " << TScount+1 << endl;
    vector<size_t> startp,countp;
    startp.push_back(TScount);
    startp.push_back(0);
    startp.push_back(0);
    countp.push_back(1);
    countp.push_back(NLAT);
    countp.push_back(NLON);
	recVar.putVar(startp,countp,outdata);
 
	// close of nc takes place in destructor
    return 0;
   }
   catch(NcException& e)
   {
      e.what();
      cout<<"FAILURE**You've probably done something wrong!!"<<endl;	  
      return NC_ERR;
   }	
  
}
