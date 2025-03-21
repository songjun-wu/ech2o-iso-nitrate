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
 * GenerateConfigTemplate.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: marcomaneta
 */

#include  <fstream>
#include <unistd.h>

#include "Sativa.h"

void GenerateConfigTemplate(const char *fn){

  ofstream ofOut;

  try{

    if (access(fn, F_OK) != -1) {

      cout << "File exists. Do you want to overwrite? (y, n):  " << endl;
      char c;
      cin.get(c);
      switch (c) {
      case 'y':
	break;
      case 'n':
	exit(EXIT_SUCCESS);
	break;
      default:
	cout << "Not a valid option. Bye" << endl;
	exit(EXIT_SUCCESS);

      }
    }

    ofOut.open(fn);
    if(!ofOut)
      throw std::ios::failure("Error opening file ");





    ofOut << "#ECH2O configuration file v1.7" << std::endl << std::endl;
    ofOut << "# Please, check Appendix A of Documentation" << std::endl;
    ofOut << "# for units of parameters and variables  " << std::endl;
    ofOut << "# (http://ech2o-iso.readthedocs.io/en/latest/Keywords.html)" << std::endl << std::endl;

    ofOut << "#" << endl << "#Folder section" << endl << "#" << endl << endl;

    ofOut << "Maps_Folder = ./Spatial" << endl;
    ofOut << "Clim_Maps_Folder = ./Climate" << endl;
    ofOut << "Output_Folder = ./Outputs" << endl << endl;

    ofOut << "#" << endl << "#Water tracking (isotopes and/or ages)" << endl;
    ofOut << "Tracking = 1" << endl ;
    ofOut << "TrackingConfig = ./configTrck.ini" << endl << endl; 

    ofOut << "#" << endl << "# Options section" << endl << "#" << endl << endl;
    
    ofOut << "MapTypes = csf" << endl;
    ofOut << "Species_State_Variable_Input_Method = tables # maps or tables" << endl << endl;

    ofOut << "#== Boolean switches" << endl;
    ofOut << "Vegetation_dynamics = 1" << endl;
    ofOut << "Reinfiltration = 1" << endl;
    ofOut << "Channel = 1" << endl ;
    ofOut << "Channel_Evaporation = 1 " << endl;
    ofOut << "# Exponential profiles: if set to 0, vertically-uniform with value equal to" << endl;
    ofOut << "# the top-of-profile open (see corresponding map inputs below)" << endl;
    ofOut << "Hydraulic_Conductivity_profile = 0" << endl ;
    ofOut << "Porosity_profile = 0" << endl << endl;
    ofOut << "#yangx 2022-05" << endl;
    ofOut << "Channelcell_partial_reinfiltration = 0" << endl << endl;
    ofOut << "#Deeper groundwater storage Yang et al., 2021 WRR" << endl;
    ofOut << "ExtraGW_Storage	= 0" << endl << endl;  //yangx 2020-05
    ofOut << "#Option to write out map information in NetCDF format" << endl;
    ofOut << "NetCDF_output_format = 0" << endl << endl;  //yangx 2021-02
    ofOut << "#option to write out the whole internal information" << endl;
    ofOut << "External_Coupling = 0" << endl << endl;  //yangx 2022-01

    ofOut << "# TOGGLE SWITCHES:" << endl;
    ofOut << "# Aerodynamic resistance choices: " << endl;
    ofOut << "# 0 = Penman Monteith option " << endl;
    ofOut << "# 1 = Thom and Oliver 1977 " << endl;
    ofOut << "Aerodyn_resist_opt = 0 " << endl << endl;

    ofOut << "# Soil resistance to vapor diffusion choices: " << endl;
    ofOut << "# 0 = No resistance" << endl;
    ofOut << "# 1 = Passerat de Silans et al. 1989" << endl;
    ofOut << "# 2 = Sellers et al. 1992" << endl;
    ofOut << "# 3 = Sakaguchi and Zeng 2009 (CLM 3.5)" << endl;
    ofOut << "Soil_resistance_opt = 3 " << endl << endl;

    ofOut << "# Surface overland flow tracer mixing with the top-layer soil water" << endl;
    ofOut << "# 0--> use the ratio of soil moisture content to soil porosity" << endl;
    ofOut << "# 1--> use the user specified ratio, a new parameter is then introduced" << endl;
    ofOut << "SrfOVFmix_mode = 1" << endl << endl ;

    ofOut << "#" << endl << "# Time variables section" << endl << "#" << endl;
    ofOut << "Simul_start = 0 # always 0" << endl;
    ofOut << "Simul_end = 31536000 # seconds (365 days)" << endl;
    ofOut << "Simul_tstep = 86400 # seconds (daily)" << endl;
    ofOut << "Clim_input_tstep = 86400 # seconds (daily)" << endl;
    ofOut << "Report_interval = 86400 # seconds (daily)" << endl ;
    ofOut << "ReportMap_interval = 86400 # seconds (daily)" << endl << endl;

    ofOut << "#" << endl;
    ofOut << "# Climate input information" << endl;
    ofOut << "# Maps in this section to be contained in folder pointed by Clim_Maps_Folder" << endl;
    ofOut << "#" << endl;
    ofOut << "ClimateZones = ClimZones.map" << endl;
    ofOut << "Snow_rain_temp_threshold = 2 " << " # Snow to rain temperatures threshold in degC" << endl;
    ofOut << "Isohyet_map = isohyet.map " << " # Precipitation multiplier map"<< endl;
    ofOut << "Precipitation = Precip.bin " << " # Precip rate in meters/second"<< endl;
    ofOut << "AirTemperature = Tavg.bin " << " # Average air temperature in degC" << endl;
    ofOut << "MaxAirTemp = Tmax.bin " << " # Maximum air temperature in degC" << endl;
    ofOut << "MinAirTemp = Tmin.bin " << " # Minimum air temperature in degC"<< endl;
    ofOut << "RelativeHumidity = RH.bin " << " # air relative humidity in kPa/kPa"<< endl;
    ofOut << "WindSpeed = windspeed.bin " << " # Wind speed in meters/second" << endl;
    ofOut << "IncomingLongWave = Ldown.bin " << " # Downwelling longwave radiation in W/sq.meter" << endl;
    ofOut << "IncomingShortWave = Sdown.bin " << " # Solar radiation in W/sq.meter" << endl << endl;

    ofOut << "#" << endl;
    ofOut << "# Spatial input information" << endl;
    ofOut << "# Maps below this line to be contained in folder pointed by Maps_Folder" << endl;
    ofOut << "#" << endl;

    ofOut << "#" << endl << "# Drainage network" << endl << "#" << endl;
    ofOut << "local_drain_direc = ldd.map" << endl;
    ofOut << "channel_width = chanwidth.map" << endl;
    ofOut << "channel_length = chanlength.map" << endl; // yangx 2020-11
    ofOut << "channel_gw_transfer_param = chanparam.map" << endl;
    ofOut << "channel_extragw_transfer_param = chanExtraparam.map" << endl; //yangx 2020-05
    ofOut << "mannings_n = chanmanningn.map" << endl << endl;
    ofOut << "# the ratio of top-layer soil water tracer signals e.g., 0 - all surface tracers, 1 - all top-layer soil water signals" << endl;
    ofOut << "overlandflow_mixing_ratio = srfovf_ratio.map" << endl << endl;

    ofOut << "#    " << endl;
    ofOut << "# Channel evaporation parameters - only if channel evaporation is on" << endl;
    ofOut << "Water_temperature = water_temp.map" << endl;
    ofOut << "Channel_roughness = chanrough.map" << endl << endl;

    ofOut << "#   " << endl;
    ofOut << "# Hydrologic Initial Conditions  " << endl;
    ofOut << "# Forest Initial states are included as maps or tables" << endl;
    ofOut << "#   " << endl;
    ofOut << "Streamflow = streamflow.map " << endl;
    ofOut << "snow_water_equivalent = swe.map " << endl;
    ofOut << "Soil_moisture_1 = Soil_moisture_1.map " << endl;
    ofOut << "Soil_moisture_2 = Soil_moisture_2.map " << endl;
    ofOut << "Soil_moisture_3 = Soil_moisture_3.map " << endl;
    ofOut << "Soil_temperature = soiltemp.map " << endl << endl;
    ofOut << "Groundwater_ExtraStorage = GW_ExtraStorage.map " << endl; //yangx 2020-05
    ofOut << "Fraction_Hydroactive_ExtraGW = ExtraGW_hydroFraction.map "  << endl << endl; //yangx 2020-05

    ofOut << "#   " << endl;
    ofOut << "#Soil parameters  " << endl;
    ofOut << "#   " << endl;
    ofOut << "DEM = DEM.map" << endl;
    ofOut << "Slope = slope.map " << endl;
    ofOut << "Total_area = total_proarea.map " << endl; //yangx 2020-11
    ofOut << "Top-of-profile_Horiz_Hydraulic_Conductivity = Keff.map " << endl;
    ofOut << "Horiz_Hydraulic_Conductivity_Profile_Coeff = kKsat.map " << endl;
    ofOut << "Vert_Horz_Anis_ratio = KvKh.map " << endl;
    ofOut << "Terrain_Random_Roughness = randrough.map " << endl;
    ofOut << "Top-of-profile_Porosity = poros.map " << endl;
    ofOut << "Porosity_Profile_Coeff = kporos.map " << endl;
    ofOut << "Air_entry_pressure = psi_ae.map " << endl;
    ofOut << "Brooks_Corey_lambda = BClambda.map " << endl;
    ofOut << "Residual_soil_moisture = theta_r.map " << endl;
    ofOut << "Soil_depth = depth_soil.map " << endl;
    ofOut << "Depth_soil_layer_1 = depth_soilL1.map " << endl;
    ofOut << "Depth_soil_layer_2 = depth_soilL2.map " << endl;
    ofOut << "Veget_water_use_param1 = Wc.map " << endl;
    ofOut << "Veget_water_use_param2 = Wp.map " << endl;
    //ofOut << "Fraction_roots_soil_layer_1 = rootfrac1.map " << endl;
    //ofOut << "Fraction_roots_soil_layer_2 = rootfrac2.map " << endl;
    //ofOut << "Root_profile_coeff = Kroot.map " << endl;
    ofOut << "Soil_bedrock_leakance = leakance.map " << endl << endl;
    ofOut << "Albedo = albedo.map" << endl;
    ofOut << "Surface_emissivity = emissivity.map" << endl;
    ofOut << "Dry_Soil_Heat_Capacity = soilheatcap.map" << endl;
    ofOut << "Dry_Soil_Therm_Cond = soilthermalK.map" << endl;
    ofOut << "Damping_depth = dampdepth.map" << endl;
    ofOut << "Temp_at_damp_depth = temp_damp.map" << endl;
    ofOut << "Snow_Melt_Coeff = snowmeltCoeff.map" << endl << endl;
    ofOut << "Sealed_area = sealed_arealproportions.map" << endl << endl;

    ofOut << "#   " << endl;
    ofOut << "#Forest Parameters and initial states " << endl;
    ofOut << "#   " << endl;
    ofOut << "ForestPatches = patches.map" << endl;
    ofOut << "Number_of_Species = 1 " << endl;
    ofOut << "Species_Parameters = SpeciesParams.tab " << endl << endl;
    ofOut << "#Tables below are only needed if Species_State_Variable_Input_Method = tables " << endl;
    ofOut << "Species_Proportion_Table = SpecsProp.tab " << endl;
    ofOut << "Species_StemDensity_Table = SpecsStemDens.tab " << endl;
    ofOut << "Species_LAI_Table = SpecsLAI.tab " << endl;
    ofOut << "Species_AGE_Table = SpecsAge.tab " << endl;
    ofOut << "Species_BasalArea_Table = SpeciesBasalArea.tab " << endl;
    ofOut << "Species_Height_table = SpeciesHeight.tab " << endl;
    ofOut << "Species_RootMass_table = SpecsRootDensity.tab " << endl << endl;

    ofOut << "#   " << endl;
    ofOut << "#Report map section " << endl;
    ofOut << "#   " << endl << endl ;
    ofOut << "Report_Long_Rad_Down = 0 " << endl;
    ofOut << "Report_Short_Rad_Down = 0 " << endl;
    ofOut << "Report_Precip = 0 " << endl;
    ofOut << "Report_Rel_Humidity = 0 " << endl;
    ofOut << "Report_Wind_Speed = 0 " << endl;
    ofOut << "Report_AvgAir_Temperature = 0 " << endl;
    ofOut << "Report_MinAir_Temperature = 0 " << endl;
    ofOut << "Report_MaxAir_Temperature = 0 " << endl << endl;

    ofOut << "Report_SWE = 0 " << endl;
    ofOut << "Report_Infilt_Cap = 0 " << endl;
    ofOut << "Report_Streamflow = 0 " << endl;
    ofOut << "Report_Saturation_Area = 0 " << endl;
    ofOut << "Report_Ponding = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_Average = 1 " << endl;
    ofOut << "Report_Soil_Water_Content_Up = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L1 = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L2 = 0 " << endl;
    ofOut << "Report_Soil_Water_Content_L3 = 0 " << endl;
    ofOut << "Report_WaterTableDepth = 0 " << endl;
    ofOut << "Report_Field_Capacity_L1 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L2 = 0 " << endl;
    ofOut << "Report_Field_Capacity_L3 = 0 " << endl;
    ofOut << "Report_Soil_Sat_Deficit = 0 " << endl;
    ofOut << "Report_Ground_Water = 0 " << endl;
    ofOut << "Report_Soil_Net_Rad = 0 " << endl;
    ofOut << "Report_Soil_LE = 0 " << endl;
    ofOut << "Report_Sens_Heat = 0 " << endl;
    ofOut << "Report_Grnd_Heat = 0 " << endl;
    ofOut << "Report_Snow_Heat = 0 " << endl;
    ofOut << "Report_Soil_Temperature = 0 " << endl;
    ofOut << "Report_Skin_Temperature = 0 " << endl << endl;

    ofOut << "Report_ChannelE = 0 " << endl;
    ofOut << "Report_Water_Temperature = 0 " << endl << endl;
	
    ofOut << "Report_Total_ET = 0 " << endl;
    ofOut << "Report_Transpiration_sum = 0 " << endl;
    ofOut << "Report_Einterception_sum = 0 " << endl;
    ofOut << "Report_Esoil_sum = 0 " << endl;
    ofOut << "Report_Net_Rad_sum = 0 " << endl ;
    ofOut << "Report_Canopy_Water_Stor_sum = 0 " << endl << endl;

    ofOut << "Report_Veget_frac = 0 " << endl;
    ofOut << "Report_Stem_Density = 0 " << endl;
    ofOut << "Report_RootFracL1_species = 0 " << endl;
    ofOut << "Report_RootFracL2_species = 0 " << endl;
    ofOut << "Report_Leaf_Area_Index = 1 " << endl;
    ofOut << "Report_Stand_Age = 0 " << endl;
    ofOut << "Report_Canopy_Conductance = 0 " << endl;
    ofOut << "Report_GPP = 0 " << endl;
    ofOut << "Report_NPP = 0 " << endl;
    ofOut << "Report_Basal_Area = 0 " << endl;
    ofOut << "Report_Tree_Height = 0 " << endl;
    ofOut << "Report_Root_Mass = 0 " << endl;
    ofOut << "Report_Canopy_Temp = 0 " << endl;
    ofOut << "Report_Canopy_NetR = 0 " << endl;
    ofOut << "Report_Canopy_LE_E = 0 " << endl;
    ofOut << "Report_Canopy_LE_T = 0 " << endl;
    ofOut << "Report_Canopy_Sens_Heat = 0 " << endl;
    ofOut << "Report_Canopy_Water_Stor = 0 " << endl;
    ofOut << "Report_species_ET = 0 " << endl;
    ofOut << "Report_Transpiration = 0 " << endl;
    ofOut << "Report_Einterception = 0 " << endl;
    ofOut << "Report_Esoil = 0 " << endl << endl;

    ofOut << "Report_GW_to_Channel = 0 " << endl;
    ofOut << "Report_Surface_to_Channel = 0 " << endl;
    ofOut << "Report_Infiltration = 0" << endl ;
    ofOut << "Report_Return_Flow_Surface = 0" << endl ;
    ofOut << "Report_Recharge_to_Layer3 = 0" << endl ;
    ofOut << "Report_Overland_Inflow = 0" << endl ;
    ofOut << "Report_Stream_Inflow = 0" << endl;
    ofOut << "Report_Groundwater_Inflow = 0 " << endl ;
    ofOut << "Report_Overland_Outflow = 0" << endl ;
    ofOut << "Report_Groundwater_Outflow = 0" << endl ;

    ofOut << "Report_GW_to_Channel_acc = 0 " << endl;
    ofOut << "Report_Surface_to_Channel_acc = 0 " << endl;
    ofOut << "Report_Infiltration_acc = 0" << endl ;
    ofOut << "Report_Return_Flow_Surface_acc = 0" << endl ;
    ofOut << "Report_Overland_Inflow_acc = 0" << endl ;
    ofOut << "Report_Stream_Inflow_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Inflow_acc = 0" << endl ;
    ofOut << "Report_Overland_Outflow_acc = 0" << endl ;
    ofOut << "Report_Groundwater_Outflow_acc = 0" << endl << endl;
    //yangx 2020-05
    ofOut << "Report_Bedrock_Leakage = 0 " << endl;
    ofOut << "Report_Extra_Ground_Water = 0 " << endl;
    ofOut << "Report_ExtraGW_to_Channel = 0 " << endl;
    ofOut << "Report_ExtraGW_to_Channel_acc = 0 " << endl;
    ofOut << "Report_ExtraGW_Inflow = 0 " << endl;
    ofOut << "Report_ExtraGW_Outflow = 0 " << endl;
    ofOut << "Report_ExtraGW_Inflow_acc = 0 " << endl;
    ofOut << "Report_ExtraGW_Outflow_acc = 0 " << endl;
    
    ofOut << "#   " << endl;
    ofOut << "#Report time series section " << endl;
    ofOut << "#   " << endl << endl;
    ofOut << "TS_mask = Tsmask.map " << endl <<"#" << endl;
    ofOut << "Ts_OutletDischarge = 1 " << endl;
    ofOut << "Ts_Long_Rad_Down = 0 " << endl;
    ofOut << "Ts_Short_Rad_Down = 0 " << endl;
    ofOut << "Ts_Precip = 0 " << endl;
    ofOut << "Ts_Rel_Humidity = 0 " << endl;
    ofOut << "Ts_Wind_Speed = 0 " << endl;
    ofOut << "Ts_AvgAir_Temperature = 0 " << endl;
    ofOut << "Ts_MinAir_Temperature = 0 " << endl;
    ofOut << "Ts_MaxAir_Temperature = 0 " << endl;
    ofOut << "Ts_SWE = 1 " << endl;
    ofOut << "Ts_Infilt_Cap = 0 " << endl;
    ofOut << "Ts_Streamflow = 0 " << endl;
    ofOut << "Ts_Ponding = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_Average = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_Up = 0 " << endl;
    ofOut << "Ts_Soil_Water_Content_L1 = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_L2 = 1 " << endl;
    ofOut << "Ts_Soil_Water_Content_L3 = 1 " << endl;
    ofOut << "Ts_WaterTableDepth = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L1 = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L2 = 0 " << endl;
    ofOut << "Ts_Field_Capacity_L3 = 0 " << endl;
    ofOut << "Ts_Soil_Sat_Deficit = 0 " << endl;
    ofOut << "Ts_Ground_Water = 0 " << endl;
    ofOut << "Ts_Soil_Net_Rad = 0 " << endl;
    ofOut << "Ts_Soil_LE = 0 " << endl;
    ofOut << "Ts_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Grnd_Heat = 0 " << endl;
    ofOut << "Ts_Snow_Heat = 0 " << endl;
    ofOut << "Ts_Soil_Temperature = 0 " << endl;
    ofOut << "Ts_Skin_Temperature = 0 " << endl << endl ;

    ofOut << "Ts_ChannelE = 0 " << endl;
    ofOut << "Ts_Water_Temperature = 0 " << endl << endl ;	

    ofOut << "Ts_Total_ET = 1 " << endl;
    ofOut << "Ts_Transpiration_sum = 1 " << endl;
    ofOut << "Ts_Einterception_sum = 1 " << endl;
    ofOut << "Ts_Esoil_sum = 1 " << endl;
    ofOut << "Ts_Net_Rad_sum = 0 " << endl ;
    ofOut << "Ts_Canopy_Water_Stor_sum = 0 " << endl << endl;

    ofOut << "Ts_Veget_frac = 0 " << endl;
    ofOut << "Ts_Stem_Density = 0 " << endl;
    ofOut << "Ts_RootFracL1_species = 0 " << endl;
    ofOut << "Ts_RootFracL2_species = 0 " << endl;
    ofOut << "Ts_Leaf_Area_Index = 1 " << endl;
    ofOut << "Ts_Stand_Age = 0 " << endl;
    ofOut << "Ts_Canopy_Conductance = 1 " << endl;
    ofOut << "Ts_GPP = 0 " << endl;
    ofOut << "Ts_NPP = 1 " << endl;
    ofOut << "Ts_Basal_Area = 0 " << endl;
    ofOut << "Ts_Tree_Height = 0 " << endl;
    ofOut << "Ts_Root_Mass = 0 " << endl;
    ofOut << "Ts_Canopy_Temp = 0 " << endl;
    ofOut << "Ts_Canopy_NetR = 0 " << endl;
    ofOut << "Ts_Canopy_LE_E = 0 " << endl;
    ofOut << "Ts_Canopy_LE_T = 0 " << endl;
    ofOut << "Ts_Canopy_Sens_Heat = 0 " << endl;
    ofOut << "Ts_Canopy_Water_Stor = 0 " << endl;
    ofOut << "Ts_species_ET = 0 " << endl;
    ofOut << "Ts_Transpiration = 0 " << endl;
    ofOut << "Ts_Einterception = 0 " << endl;
    ofOut << "Ts_Esoil = 0 " << endl << endl;

    ofOut << "Ts_GW_to_Channel = 0 " << endl;
    ofOut << "Ts_Surface_to_Channel = 0 " << endl;
    ofOut << "Ts_Infiltration = 0" << endl ;
    ofOut << "Ts_Return_Flow_Surface = 0" << endl ;
    ofOut << "Ts_Recharge_to_Layer3 = 0" << endl ;
    ofOut << "Ts_Overland_Inflow = 0" << endl ;
    ofOut << "Ts_Stream_Inflow = 0" << endl;
    ofOut << "Ts_Groundwater_Inflow = 0 " << endl ;
    ofOut << "Ts_Overland_Outflow = 0" << endl ;
    ofOut << "Ts_Groundwater_Outflow = 0" << endl ;

    ofOut << "Ts_GW_to_Channel_acc = 0 " << endl;
    ofOut << "Ts_Surface_to_Channel_acc = 0 " << endl;
    ofOut << "Ts_Infiltration_acc = 0" << endl ;
    ofOut << "Ts_Return_Flow_Surface_acc = 0" << endl ;
    ofOut << "Ts_Overland_Inflow_acc = 0" << endl ;
    ofOut << "Ts_Stream_Inflow_acc = 0" << endl ;
    ofOut << "Ts_Groundwater_Inflow_acc = 0" << endl ;
    ofOut << "Ts_Overland_Outflow_acc = 0" << endl ;
    ofOut << "Ts_Groundwater_Outflow_acc = 0" << endl;
    //yangx 2020-05
	
    ofOut << "Ts_Bedrock_Leakage = 0 " << endl;
    ofOut << "Ts_Extra_Ground_Water = 0 " << endl;
    ofOut << "Ts_ExtraGW_to_Channel = 0 " << endl;
    ofOut << "Ts_ExtraGW_to_Channel_acc = 0 " << endl;
    ofOut << "Ts_ExtraGW_Inflow = 0 " << endl;
    ofOut << "Ts_ExtraGW_Outflow = 0 " << endl;
    ofOut << "Ts_ExtraGW_Inflow_acc = 0 " << endl;
    ofOut << "Ts_ExtraGW_Outflow_acc = 0 " << endl;


    if (ofOut)
      ofOut.close();
  }
  catch(const std::exception &e){
    cout << "Failure writing configuration template file with  " << e.what() << endl;
    exit(EXIT_FAILURE);
  }
}

