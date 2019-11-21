/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "lambertTargetingTools.h"

using namespace tudat;
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::basic_mathematics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::unit_conversions;
using namespace tudat::coordinate_conversions;
using namespace tudat::estimatable_parameters;
using namespace tudat::statistics;
using namespace tudat::json_interface;
using namespace tudat::mission_segments;

std::string getCurrentRootPath( )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string currentPath = filePath_.substr( 0, filePath_.length( ) -
                                                ( std::string( "lambertTargetingTools.cpp" ).length( ) ) );

    return currentPath;
}

std::shared_ptr< Ephemeris > getLambertProblemResult(
        const NamedBodyMap& bodyMap,
        const std::string targetBody,
        const double departureTime,
        const double arrivalTime )
{
    // Gravitational parameter of the Sun
    double centralBodyGravitationalParameter = bodyMap.at( "Sun" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set initial and final positions of Lambert targeter.
    Eigen::Vector6d initialState = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( departureTime ) -
            bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris( departureTime );
    Eigen::Vector6d finalState = bodyMap.at( targetBody )->getStateInBaseFrameFromEphemeris( arrivalTime ) -
            bodyMap.at( "Sun" )->getStateInBaseFrameFromEphemeris( arrivalTime );


    // Setup Lambert targeter algorithm
    LambertTargeterIzzo lambertTargeter(
                initialState.segment( 0, 3 ), finalState.segment( 0, 3 ),
                arrivalTime - departureTime, centralBodyGravitationalParameter );

    // Compute initial Cartesian state of Lambert arc
    Eigen::Vector6d lambertArcInitialState = initialState;
    lambertArcInitialState.segment( 3, 3 ) = lambertTargeter.getInertialVelocityAtDeparture( );

    // Compute Keplerian state of Lambert arc
    Eigen::Vector6d lambertArcKeplerianElements = convertCartesianToKeplerianElements(
                lambertArcInitialState, centralBodyGravitationalParameter );

    // Setup Kepler ephemeris model
    std::shared_ptr< Ephemeris > keplerEphemeris = std::make_shared< KeplerEphemeris >(
                lambertArcKeplerianElements, departureTime, centralBodyGravitationalParameter );

    // Create interpolator for Cartesian states
    return keplerEphemeris;

}

std::map< double, Eigen::VectorXd > getLambertArcHistory(
        const std::shared_ptr< Ephemeris > lambertArcStateModel,
        const std::map< double, Eigen::VectorXd >& numericalStateHistory )
{
    std::map< double, Eigen::VectorXd > lambertArcStates;
    for( auto stateIterator : numericalStateHistory )
    {
        lambertArcStates[ stateIterator.first ] = lambertArcStateModel->getCartesianState(
                    stateIterator.first );
    }
    return lambertArcStates;
}

void printEphemerides(
        const double departureEpoch, const double arrivalEpoch, const NamedBodyMap& bodyMap,
        const std::string targetBody )
{
    std::string outputFolder = getCurrentRootPath( ) + "SimulationOutput";

    double interval = 3600.0;
    std::map< double, Eigen::VectorXd > earthStates;
    std::map< double, Eigen::VectorXd > targetStates;

    double currentTime = departureEpoch - tudat::physical_constants::JULIAN_YEAR;
    while( currentTime < arrivalEpoch + tudat::physical_constants::JULIAN_YEAR )
    {
        earthStates[ currentTime ] = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( currentTime );
        targetStates[ currentTime ] = bodyMap.at( targetBody )->getStateInBaseFrameFromEphemeris( currentTime );
        currentTime += interval;
    }

    input_output::writeDataMapToTextFile(
                earthStates, "earth_states.dat", outputFolder );
    input_output::writeDataMapToTextFile(
                targetStates, "target_states.dat", outputFolder );

}

//! Create propagator settings for unperturbed arc
std::shared_ptr< TranslationalStatePropagatorSettings< > > getUnperturbedPropagatorSettings(
        const NamedBodyMap& bodyMap,
        const nlohmann::json& jsonObject,
        const Eigen::Vector6d& initialState,
        const double terminationTime )
{
    // Define central body (Sun) and body undergoing accelerations (Spacecraft)
    std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    std::vector< std::string > centralBodies = { "Sun" };

    // Define accelerations: Sun only
    simulation_setup::SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    basic_astrodynamics::AccelerationMap accelerationsMap_ =
            simulation_setup::createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    // Retrieve list of dependent variables from json object
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariableList =
            getValue< std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > >(
                jsonObject, "dependentVariableSettings" );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariableList, false );

    // Create propagator settings
    return std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap_, bodiesToIntegrate,
                initialState, std::make_shared< PropagationTimeTerminationSettings >(
                    terminationTime, true ), cowell, dependentVariablesToSave );
}

//! Create propagator settings for perturbed arc
std::shared_ptr< TranslationalStatePropagatorSettings< > > getPerturbedPropagatorSettings(
        const NamedBodyMap& bodyMap,
        const nlohmann::json& jsonObject,
        const Eigen::Vector6d& initialState,
        const double terminationTime )
{
    // Define central body (SSB) and body undergoing accelerations (Spacecraft)
    std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    std::vector< std::string > centralBodies = { "Sun" };

    // Retrieve acceleration settings from json object
    simulation_setup::SelectedAccelerationMap accelerationSettingsMap =
            getValue< SelectedAccelerationMap >( jsonObject, "accelerations" );
    basic_astrodynamics::AccelerationMap accelerationsMap_ =
            simulation_setup::createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    // Retrieve list of dependent variables from json object
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariableList =
            getValue< std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > >( jsonObject, "dependentVariableSettings" );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariableList, false );

    // Create propagator settings
    return std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap_, bodiesToIntegrate,
                initialState, std::make_shared< PropagationTimeTerminationSettings >(
                    terminationTime, true ), cowell, dependentVariablesToSave );
}

void writePropagationResultsToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::shared_ptr< Ephemeris > lambertArcStateModel,
        const std::string& fileOutputIdentifier )
{
    std::string outputFolder = getCurrentRootPath( ) + "SimulationOutput";

    input_output::writeDataMapToTextFile(
                dynamicsSimulator->getEquationsOfMotionNumericalSolution( ), fileOutputIdentifier + "_numerical_states.dat", outputFolder );
    input_output::writeDataMapToTextFile(
                dynamicsSimulator->getDependentVariableHistory( ), fileOutputIdentifier + "_dependent.dat", outputFolder );

    std::map< double, Eigen::VectorXd > lambertArcStates =
            getLambertArcHistory( lambertArcStateModel, dynamicsSimulator->getEquationsOfMotionNumericalSolution( ) );
    input_output::writeDataMapToTextFile(
                lambertArcStates, fileOutputIdentifier + "_lambert_states.dat", outputFolder );
}

std::shared_ptr< EstimatableParameterSet< > >  getSensitityParameterSet(
        const std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings,
        const NamedBodyMap& bodyMap )
{
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "Spacecraft", propagatorSettings->getInitialStates( ), "SSB" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Spacecraft", radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );

    return createParametersToEstimate( parameterNames, bodyMap );
}

std::shared_ptr< SingleArcDynamicsSimulator< > > propagateTrajectory(
        const double initialTime, const double finalTime, const nlohmann::json& jsonObject,
        const NamedBodyMap& bodyMap, const std::shared_ptr< Ephemeris > lambertArcStateModel,
        const std::string& fileOutputIdentifier, const bool usePerturbations )
{
    // Compute initial state
    Eigen::Vector6d lambertArcInitialState = lambertArcStateModel->getCartesianState( initialTime );

    // Get propagator settings for perturbed forward and backward arcs
    std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings;
    if( usePerturbations )
    {
        propagatorSettings = getPerturbedPropagatorSettings( bodyMap, jsonObject, lambertArcInitialState, finalTime );
    }
    else
    {
        propagatorSettings = getUnperturbedPropagatorSettings( bodyMap, jsonObject, lambertArcInitialState, finalTime );
    }

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings;
    updateFromJSON( integratorSettings, jsonObject, Keys::integrator );
    integratorSettings->initialTime_ = initialTime;

    if( initialTime > finalTime )
    {
        integratorSettings->initialTimeStep_ *= -1.0;
    }

    // Propagate forward perturbed arc and save results
    std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator = std::make_shared< SingleArcDynamicsSimulator< > >(
                bodyMap, integratorSettings, propagatorSettings );
    writePropagationResultsToFile( dynamicsSimulator, lambertArcStateModel, fileOutputIdentifier );

    return dynamicsSimulator;

}

std::shared_ptr< SingleArcVariationalEquationsSolver< > > propagateVariationalEquations(
        const double initialTime, const double finalTime, const nlohmann::json& jsonObject,
        const NamedBodyMap& bodyMap, const std::shared_ptr< Ephemeris > lambertArcStateModel )
{
    // Compute initial state
    Eigen::Vector6d lambertArcInitialState = lambertArcStateModel->getCartesianState( initialTime );

    // Get propagator settings for perturbed forward and backward arcs
    std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
            getPerturbedPropagatorSettings( bodyMap, jsonObject, lambertArcInitialState, finalTime );


    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings;
    updateFromJSON( integratorSettings, jsonObject, Keys::integrator );
    integratorSettings->initialTime_ = initialTime;

    if( initialTime > finalTime )
    {
        integratorSettings->initialTimeStep_ *= -1.0;
    }

    std::shared_ptr< EstimatableParameterSet< > > parametersForWhichToComputeSensitivity =
            getSensitityParameterSet( propagatorSettings, bodyMap );
    printEstimatableParameterEntries( parametersForWhichToComputeSensitivity );

    std::shared_ptr< SingleArcVariationalEquationsSolver< > > variationalEquationsSolver =
            std::make_shared< SingleArcVariationalEquationsSolver< > >(
                bodyMap, integratorSettings,
                propagatorSettings, parametersForWhichToComputeSensitivity,
                true, nullptr, false, true, false );

    return variationalEquationsSolver;

}
