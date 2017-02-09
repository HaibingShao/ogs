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
            BaseLib::ConfigTree const& config)
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
                //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__thermal_conductivity_solid}
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
