/**
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/

#include "CreateHeatTransportBHEProcess.h"

#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "HeatTransportBHEProcess.h"
#include "HeatTransportBHEProcessData.h"
#include "BHE/BHEAbstract.h"
#include "BHE/BHE_1U.h"
#include "BHE/BHE_Net.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/Density/createFluidDensityModel.h"
#include "MaterialLib/Fluid/Viscosity/createViscosityModel.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"
#include "BaseLib/reorderVector.h"

namespace ProcessLib
{
    namespace HeatTransportBHE
    {
        std::unique_ptr<Process> createHeatTransportBHEProcess(
            MeshLib::Mesh& mesh,
            std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
            std::vector<ProcessVariable> const& variables,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            unsigned const integration_order,
            BaseLib::ConfigTree const& config, 
            std::map<std::string,
                     std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const& curves)
        {
            //! \ogs_file_param{prj__processes__process__type}
            config.checkConfigParameter("type", "HEAT_TRANSPORT_BHE");

            DBUG("Create HeatTransportBHEProcess.");

            // Process variable.

            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables}
            auto const pv_config = config.getConfigSubtree("process_variables");
            
            auto process_variables = findProcessVariables(
                variables, pv_config,
                {//! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__process_variables__process_variable}
                    "process_variable" });

            // solid phase thermal conductivity parameter.
            auto& thermal_conductivity_solid = findParameter<double>(
                config,
                //! \ogs_file_param_special prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_solid}
                "thermal_conductivity_solid", parameters, 1);

            DBUG("Use \'%s\' as solid phase thermal conductivity parameter.",
                thermal_conductivity_solid.name.c_str());

            // solid phase thermal conductivity parameter.
            auto& thermal_conductivity_fluid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_fluid}
                "thermal_conductivity_fluid", parameters, 1);

            DBUG("Use \'%s\' as fluid phase thermal conductivity parameter.",
                thermal_conductivity_fluid.name.c_str());

            // gas phase thermal conductivity parameter.
            auto& thermal_conductivity_gas = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_gas}
                "thermal_conductivity_gas", parameters, 1);

            DBUG("Use \'%s\' as gas phase thermal conductivity parameter.",
                thermal_conductivity_gas.name.c_str());

            // solid phase heat capacity parameter.
            auto& heat_capacity_solid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_solid}
                "heat_capacity_solid", parameters, 1);

            DBUG("Use \'%s\' as solid phase heat capacity parameter.", heat_capacity_solid.name.c_str());

            // fluid phase heat capacity parameter.
            auto& heat_capacity_fluid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_fluid}
                "heat_capacity_fluid", parameters, 1);

            DBUG("Use \'%s\' as fluid phase heat capacity parameter.", heat_capacity_fluid.name.c_str());

            // gas phase heat capacity parameter.
            auto& heat_capacity_gas = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__heat_capacity_gas}
                "heat_capacity_gas", parameters, 1);

            DBUG("Use \'%s\' as gas phase heat capacity parameter.", heat_capacity_gas.name.c_str());

            // solid phase density parameter.
            auto& density_solid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_solid}
                "density_solid", parameters, 1);

            DBUG("Use \'%s\' as solid phase density parameter.", density_solid.name.c_str());

            // fluid phase density parameter.
            auto& density_fluid = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_fluid}
                "density_fluid", parameters, 1);

            DBUG("Use \'%s\' as fluid phase density parameter.", density_fluid.name.c_str());

            // gas phase density parameter.
            auto& density_gas = findParameter<double>(
                config,
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__density_gas}
                "density_gas", parameters, 1);

            DBUG("Use \'%s\' as gas phase density parameter.", density_gas.name.c_str());

            // reading BHE parameters--------------------------------------------------------------
            std::vector<BHE::BHEAbstract*> vec_BHEs;
            BHE::BHE_Net BHE_network;

            // now read the BHE configurations
            auto const& bhe_configs =
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers}
                config.getConfigSubtree("borehole_heat_exchangers");

            for (
                auto const& bhe_conf :
                //! \ogs_file_param prj__processes__process___HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger}
                bhe_configs.getConfigSubtreeList("borehole_heat_exchanger"))
            {
                auto const bhe_id = bhe_conf.getConfigAttribute<int>("id");
                
                // read in the parameters
                using namespace BHE; 
                const std::string bhe_type_str = bhe_conf.getConfigParameter<std::string>("bhe_type");
                const std::string bhe_ply_name = bhe_conf.getConfigParameter<std::string>("bhe_polyline");
                const std::string bhe_bound_type_str = bhe_conf.getConfigParameter<std::string>("bhe_bound_type");
                const bool bhe_if_use_flow_rate_curve = false; 
                const double bhe_length = bhe_conf.getConfigParameter<double>("bhe_length");
                const double bhe_diameter = bhe_conf.getConfigParameter<double>("bhe_diameter");
                const double bhe_refrigerant_flow_rate = bhe_conf.getConfigParameter<double>("bhe_refrigerant_flow_rate");
                const double bhe_pipe_inner_radius = bhe_conf.getConfigParameter<double>("bhe_pipe_inner_radius");
                const double bhe_pipe_outer_radius = bhe_conf.getConfigParameter<double>("bhe_pipe_outer_radius");
                const double bhe_pipe_in_wall_thickness = bhe_conf.getConfigParameter<double>("bhe_pipe_in_wall_thickness");
                const double bhe_pipe_out_wall_thickness = bhe_conf.getConfigParameter<double>("bhe_pipe_out_wall_thickness");
                const std::size_t bhe_fluid_idx = bhe_conf.getConfigParameter<std::size_t>("bhe_fluid_idx");
                const double bhe_fluid_longitudinal_dispsion_length = bhe_conf.getConfigParameter<double>("bhe_fluid_longitudinal_dispsion_length"); 
                const double bhe_grout_density = bhe_conf.getConfigParameter<double>("bhe_grout_density");
                const double bhe_grout_porosity = bhe_conf.getConfigParameter<double>("bhe_grout_porosity");
                const double bhe_grout_heat_capacity = bhe_conf.getConfigParameter<double>("bhe_grout_heat_capacity");
                const double bhe_pipe_wall_thermal_conductivity = bhe_conf.getConfigParameter<double>("bhe_pipe_wall_thermal_conductivity");
                const double bhe_grout_thermal_conductivity = bhe_conf.getConfigParameter<double>("bhe_grout_thermal_conductivity");
                const double bhe_pipe_distance = bhe_conf.getConfigParameter<double>("bhe_pipe_distance");

                // optional parameters
                double bhe_power_in_watt_val;
                double bhe_intern_resistance; 
                double bhe_therm_resistance; 
                double bhe_R_fig; 
                double bhe_R_fog; 
                double bhe_R_gg1; 
                double bhe_R_gg2; 
                double bhe_R_gs; 
                double bhe_delta_T_val; 
                double bhe_switch_off_threshold; 

                // give default values to optional parameters
                // if the BHE is using external given thermal resistance values
                bool bhe_use_ext_therm_resis = false;
                if (auto const bhe_use_ext_therm_resis_conf =
                    config.getConfigParameterOptional<bool>("bhe_use_external_therm_resis"))
                {
                    DBUG("If using external given thermal resistance values : %s",
                        (*bhe_use_ext_therm_resis_conf) ? "true" : "false");
                    bhe_use_ext_therm_resis = *bhe_use_ext_therm_resis_conf;

                    // only when using it, read the two resistance values
                    if (bhe_use_ext_therm_resis)
                    {
                        bhe_intern_resistance = *bhe_conf.getConfigParameterOptional<double>("bhe_internal_resistance");

                        bhe_therm_resistance = bhe_conf.getConfigParameterOptional<double>("bhe_therm_resistance").get();
                    }
                }

                // if the BHE is using user defined thermal resistance values
                bool bhe_user_defined_therm_resis = false; 
                if (auto const bhe_user_defined_therm_resis_conf =
                    config.getConfigParameterOptional<bool>("bhe_user_defined_therm_resis"))
                {
                    DBUG("If appplying user defined thermal resistance values : %s",
                        (*bhe_user_defined_therm_resis_conf) ? "true" : "false");
                    bhe_user_defined_therm_resis = *bhe_user_defined_therm_resis_conf;

                    // only when using it, read the two resistance values
                    if (bhe_user_defined_therm_resis)
                    {
                        bhe_R_fig = bhe_conf.getConfigParameterOptional<double>("bhe_R_fig").get();
                        bhe_R_fog = bhe_conf.getConfigParameterOptional<double>("bhe_R_fog").get();
                        bhe_R_gg1 = bhe_conf.getConfigParameterOptional<double>("bhe_R_gg1").get();
                        bhe_R_gg2 = bhe_conf.getConfigParameterOptional<double>("bhe_R_gg2").get();
                        bhe_R_gs = bhe_conf.getConfigParameterOptional<double>("bhe_R_gs").get();
                    }
                }

                // convert BHE type
                BHE::BHE_TYPE bhe_type; 
                if (bhe_type_str =="BHE_TYPE_1U")
                    bhe_type = BHE_TYPE::BHE_TYPE_1U;
                else if (bhe_type_str.compare("BHE_TYPE_2U") == 0)
                    bhe_type = BHE_TYPE::BHE_TYPE_2U;
                else if (bhe_type_str.compare("BHE_TYPE_CXC") == 0)
                    bhe_type = BHE_TYPE::BHE_TYPE_CXC;
                else if (bhe_type_str.compare("BHE_TYPE_CXA") == 0)
                    bhe_type = BHE_TYPE::BHE_TYPE_CXA;

                // convert BHE boundary type
                BHE::BHE_BOUNDARY_TYPE bhe_bound_type; 
                if (bhe_bound_type_str.compare("FIXED_INFLOW_TEMP") == 0)
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_FIXED_INFLOW_TEMP;
                else if (bhe_bound_type_str.compare("FIXED_INFLOW_TEMP_CURVE") == 0)
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_FIXED_INFLOW_TEMP_CURVE;
                else if (bhe_bound_type_str.compare("POWER_IN_WATT") == 0)
                {
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_POWER_IN_WATT;
                    bhe_power_in_watt_val = bhe_conf.getConfigParameterOptional<double>("bhe_power_in_watt_value").get();
                    bhe_switch_off_threshold = bhe_conf.getConfigParameterOptional<double>("bhe_switch_off_threshold").get();
                }
                else if (bhe_bound_type_str.compare("POWER_IN_WATT_CURVE_FIXED_DT") == 0)
                {
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_DT;
                    bhe_switch_off_threshold = bhe_conf.getConfigParameterOptional<double>("bhe_switch_off_threshold").get();
                }
                else if (bhe_bound_type_str.compare("BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_DT") == 0)
                {
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_DT;
                    bhe_switch_off_threshold = bhe_conf.getConfigParameterOptional<double>("bhe_switch_off_threshold").get();
                }
                else if (bhe_bound_type_str.compare("BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE") == 0)
                {
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_BUILDING_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE;
                    bhe_switch_off_threshold = bhe_conf.getConfigParameterOptional<double>("bhe_switch_off_threshold").get();
                }
                else if (bhe_bound_type_str.compare("POWER_IN_WATT_CURVE_FIXED_FLOW_RATE") == 0)
                {
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_POWER_IN_WATT_CURVE_FIXED_FLOW_RATE;
                    bhe_switch_off_threshold = bhe_conf.getConfigParameterOptional<double>("bhe_switch_off_threshold").get();
                }
                else if (bhe_bound_type_str.compare("FIXED_TEMP_DIFF") == 0)
                {
                    bhe_bound_type = BHE_BOUNDARY_TYPE::BHE_BOUND_FIXED_TEMP_DIFF;
                    bhe_delta_T_val = bhe_conf.getConfigParameterOptional<double>("bhe_inout_delta_T_value").get();
                }

                // get the refrigerant properties from fluid property class
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__fluid}
                auto const& fluid_config = config.getConfigSubtree("fluid");
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_density}
                auto const& rho_conf = fluid_config.getConfigSubtree("refrigerant_density");
                auto bhe_refrigerant_density =
                    MaterialLib::Fluid::createFluidDensityModel(rho_conf);
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_viscosity}
                auto const& mu_conf = fluid_config.getConfigSubtree("refrigerant_viscosity");
                auto bhe_refrigerant_viscosity =
                    MaterialLib::Fluid::createViscosityModel(mu_conf);
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_specific_heat_capacity}
                auto const& cp_conf = fluid_config.getConfigSubtree("refrigerant_specific_heat_capacity");
                auto bhe_refrigerant_heat_capacity =
                    MaterialLib::Fluid::createSpecificFluidHeatCapacityModel(cp_conf);
                //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__material_property__refrigerant_thermal_conductivity}
                auto const& lambda_conf = fluid_config.getConfigSubtree("refrigerant_thermal_conductivity");
                auto bhe_regrigerant_heat_conductivity =
                    MaterialLib::Fluid::createFluidThermalConductivityModel(lambda_conf);

                MaterialLib::Fluid::FluidProperty::ArrayType vars; 
                vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::T)] = 298.15;
                vars[static_cast<int>(MaterialLib::Fluid::PropertyVariableType::p)] = 101325.0;

                // initialize the BHE class
                switch (bhe_type)
                {
                case BHE_TYPE::BHE_TYPE_1U:
                    BHE::BHE_1U * m_bhe_1u; 
                    m_bhe_1u = new BHE::BHE_1U(bhe_ply_name, bhe_bound_type, bhe_use_ext_therm_resis,
                        bhe_user_defined_therm_resis, curves, bhe_length,
                        bhe_diameter, bhe_refrigerant_flow_rate, bhe_pipe_inner_radius, 
                        bhe_pipe_outer_radius, bhe_pipe_in_wall_thickness, bhe_pipe_out_wall_thickness, 
                        bhe_refrigerant_viscosity->getValue(vars), bhe_refrigerant_density->getValue(vars), bhe_fluid_longitudinal_dispsion_length,
                        bhe_refrigerant_heat_capacity->getValue(vars), bhe_grout_density, bhe_grout_porosity, 
                        bhe_grout_heat_capacity, bhe_regrigerant_heat_conductivity->getValue(vars), bhe_pipe_wall_thermal_conductivity, 
                        bhe_grout_thermal_conductivity, bhe_pipe_distance, bhe_power_in_watt_val, 
                        bhe_delta_T_val, bhe_intern_resistance, bhe_therm_resistance,
                        bhe_R_fig, bhe_R_fog, bhe_R_gg1,
                        bhe_R_gg2, bhe_R_gs, bhe_if_use_flow_rate_curve, 
                        bhe_switch_off_threshold);

                    vec_BHEs.push_back(m_bhe_1u);
                    BHE_network.add_bhe_net_elem(m_bhe_1u);

                    // now adding a pipeline connecting the bottom of this BHE
                    BHE::BHE_Net_ELE_Pipe_Inner_1U * m_bhe_pipe_1u;
                    m_bhe_pipe_1u = new BHE::BHE_Net_ELE_Pipe_Inner_1U(m_bhe_1u->get_ele_name().append("_INNER_PIPE"), m_bhe_1u);
                    BHE_network.add_bhe_net_pipe(m_bhe_pipe_1u,
                        m_bhe_1u->get_ele_name(),
                        0,
                        m_bhe_1u->get_ele_name(),
                        0);

                    break;
                case BHE_TYPE::BHE_TYPE_2U:
                    // TODO
                    break; 
                case BHE_TYPE::BHE_TYPE_CXA:
                    // TODO
                    break; 
                case BHE_TYPE::BHE_TYPE_CXC:
                    // TODO
                    break;
                }

            }
            // end of reading BHE parameters-------------------------------------------------------

            HeatTransportBHEProcessData process_data{ thermal_conductivity_solid, 
                thermal_conductivity_fluid, 
                thermal_conductivity_gas, 
                heat_capacity_solid,
                heat_capacity_fluid,
                heat_capacity_gas,
                density_solid, 
                density_fluid,
                density_gas};

            SecondaryVariableCollection secondary_variables;

            NumLib::NamedFunctionCaller named_function_caller(
            { "HeatConduction_temperature" });

            ProcessLib::parseSecondaryVariables(config, secondary_variables,
                named_function_caller);

            return std::unique_ptr<Process>{new HeatTransportBHEProcess{
                mesh, std::move(jacobian_assembler), parameters, integration_order,
                std::move(process_variables), std::move(process_data),
                std::move(secondary_variables), std::move(named_function_caller) }};
        }

    }  // namespace HeatTransportBHE
}  // namespace ProcessLib
