#ifndef LAMBERTTARGETINGTOOLS_H
#define LAMBERTTARGETINGTOOLS_H

#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"
#include "Tudat/JsonInterface/jsonInterface.h"

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

std::string getCurrentRootPath( );

/*!
 * Run a Lambert targeter instance from Earth to target body subjected to the departure and arrival times.
 * \param bodyMap Body map containing all simulation bodies
 * \param targetBody String describing the Lambert targeter target body
 * \param departureTime The departure epoch [s]
 * \param arrivalTime The arrival epoch in seconds [s]
 * \return Shared pointer to Ephemeris containing the Lambert targeter states.
 */
std::shared_ptr< Ephemeris > getLambertProblemResult(
        const NamedBodyMap& bodyMap,
        const std::string targetBody,
        const double departureTime,
        const double arrivalTime );

std::map< double, Eigen::VectorXd > getLambertArcHistory(
        const std::shared_ptr< Ephemeris > lambertArcStateModel,
        const std::map< double, Eigen::VectorXd >& numericalStateHistory );

void printEphemerides(
        const double departureEpoch, const double arrivalEpoch, const NamedBodyMap& bodyMap,
        const std::string targetBody );

std::shared_ptr< TranslationalStatePropagatorSettings< > > getUnperturbedPropagatorSettings(
        const NamedBodyMap& bodyMap,
        const nlohmann::json& jsonObject,
        const Eigen::Vector6d& initialState,
        const double terminationTime );

std::shared_ptr< TranslationalStatePropagatorSettings< > > getPerturbedPropagatorSettings(
        const NamedBodyMap& bodyMap,
        const nlohmann::json& jsonObject,
        const Eigen::Vector6d& initialState,
        const double terminationTime );

/*! Function to write numerical propagation results to three different output files:
 *  numerical Cartesian states, dependent variables, and the Lambert targeter Cartesian states.
 *
 * THIS FUNCTION OUTPUTS THE FOLLOWING FILES:
 *
 *  -  fileOutputIdentifier_numerical_states.dat Numerical integrated state history of the spacecraft
 *  -  fileOutputIdentifier_dependent.dat Dependent variable history obtained from the numerical propagation
 *  -  fileOutputIdentifierlambert_states.dat States evaluated from the Lambert arc model, at the same epochs as the
 *                                            numerical propagated states
 *
 *  Where 'fileOutputIdentifier' is a string replaced by the value of the associated input variable to this function.
 *  All files are written to the src2/SimulationOutput directory.
 *
 * \param dynamicsSimulator Shared pointer to the dynamics simulator.
 * \param lambertArcStateModel Shared pointer to the Lambert target arc state model.
 * \param filePrefix String that is prepended to the name of the output files.
 */
void writePropagationResultsToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::shared_ptr< Ephemeris > lambertArcStateModel,
        const std::string& fileOutputIdentifier );

std::shared_ptr< EstimatableParameterSet< > >  getSensitityParameterSet(
        const std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings,
        const NamedBodyMap& bodyMap );

/*!
 *  Numerically propagate a satellite subjected to the dynamical model as specified in the JSON object.
 *
 *  NOTE: THIS FUNCTION PROVIDES OUTPUT FILES BY CALLING THE writePropagationResultsToFile FUNCTION
 *
 * \param initialTime Initial time of the propagation [s]
 * \param finalTime Final time of the propagation [s]
 * \param jsonObject JSON object containing the simulation settings
 * \param bodyMap Body map containing all simulation bodies
 * \param lambertArcStateModel Lambert targeter model containing the computed states
 * \param fileOutputIdentifier Identifier string that is prepended to the three output files
 * \param usePerturbations Boolean value specifying whether or not perturbations should be included
 * \return Shared pointer to the single arc dynamics simulator
 */
std::shared_ptr< SingleArcDynamicsSimulator< > > propagateTrajectory(
        const double initialTime, const double finalTime, const nlohmann::json& jsonObject,
        const NamedBodyMap& bodyMap, const std::shared_ptr< Ephemeris > lambertArcStateModel,
        const std::string& fileOutputIdentifier, const bool usePerturbations );

/*!
 * Numerically propagate the variational equations subject to the dynamical model as specified in the JSON object.
 * \param initialTime Initial time of the propagation [s]
 * \param finalTime Final time of the propagation [s]
 * \param jsonObject JSON object containing the simulation settings
 * \param bodyMap Body map containing all simulation bodies
 * \param lambertArcStateModel Lambert targeter model containing the computed states
 * \return Shared pointer to the single arc variational equations solver
 */
std::shared_ptr< SingleArcVariationalEquationsSolver< > > propagateVariationalEquations(
        const double initialTime, const double finalTime, const nlohmann::json& jsonObject,
        const NamedBodyMap& bodyMap, const std::shared_ptr< Ephemeris > lambertArcStateModel );

#endif // LAMBERTTARGETINGTOOLS_H
