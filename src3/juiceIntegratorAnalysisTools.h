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

#ifndef JUICEINTEGRATORANALYSISTOOLS_H
#define JUICEINTEGRATORANALYSISTOOLS_H

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

/*! Function to get the approxomate pericenter time of a propagated orbit. E.g. in a flyby, time of closest approach.
 *
 * \param numericalSolution Map containing the numerical propagation results.
 * \return Epoch of pericenter.
 */
double getClosestApproachTime( std::map< double, Eigen::VectorXd > numericalSolution );

/*! Function that returns the integrator settings with a fixed time step RKF78 integrator.
 * \param initialTime Initial time of the propagation.
 * \param timeStep Size of the integrator time step.
 * \return Shared pointer to an IntegratorSettings object containing the RKF78 fixed time step settings.
 */
std::shared_ptr< IntegratorSettings< double > > getFixedStepSizeIntegratorSettings(
        const double initialTime, const double timeStep );

/*! Function to write numerical propagation results to two different output files:
 *  numerical Cartesian states and dependent variables
 *
 * THIS FUNCTION OUTPUTS THE FOLLOWING FILES INTO THE src3/SimulationOutput FOLDER:
 *
 *  -  fileOutputIdentifier_numerical_states.dat Numerical integrated state history of the spacecraft
 *  -  fileOutputIdentifier_dependent.dat Dependent variable history obtained from the numerical propagation
 *
 *  Where 'fileOutputIdentifier' is a string replaced by the value of the associated input variable to this function.
 *
 * \param dynamicsSimulator Shared pointer to the dynamics simulator.
 * \param fileOutputIdentifier String that is prepended to the name of the output files.
 */
void writePropagationResultsToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::string& fileOutputIdentifier );

/*! Function to retrieve the acceleration map pertaining to the unperturbed orbit model (just the central body's point gravity).
 *
 * \param centralBody String holding the current central body.
 * \param bodyMap Map with bodies in the simulation. Must contain all bodies that are propagated or exert accelerations.
 * \return AccelerationMap object containing the accelerations on Juice for the unperturbed model.
 */
AccelerationMap getUnperturbedAccelerations(
        const std::string& centralBody, const NamedBodyMap& bodyMap );

/*! Function to retrieve the acceleration map pertaining to the perturbed orbit model (3rd bodies + aerodynamics,
 * spherical harmonics).
 *
 * \param centralBody String holding the current central body.
 * \param bodyMap Map with bodies in the simulation. Must contain all bodies that are propagated or exert accelerations.
 * \return AccelerationMap object containing the accelerations on Juice for the perturbed model.
 */
AccelerationMap getPerturbedAccelerations(
        const std::string& centralBody, const NamedBodyMap& bodyMap );

std::map< double, Eigen::VectorXd > getDifferenceWrtKeplerOrbit(
        const std::map< double, Eigen::VectorXd >& numericalSolution,
        const double centralBodyGravitationalParameter );

/*! Function to write the numerically propagated orbit and its difference w.r.t. to a benchmark solution to a file.
 *
 *  THIS FUNCTION OUTPUTS THE FOLLOWING FILES INTO THE src3/SimulationOutput FOLDER:
 *
 *  -  fileOutputIdentifier_numerical_states.dat Numerical integrated state history of the spacecraft
 *  -  fileOutputIdentifier_benchmark_difference.dat State difference history between numerical orbit and benchmark
 *  -  fileOutputIdentifier_dependent_states.dat Dependent variable history (if any) for propagation
 *
 *  Where 'fileOutputIdentifier' is a string replaced by the value of the associated input variable to this function.
 *
 * \param dynamicsSimulator Shared pointer to a SingleArcDynamicsSimulator object.
 * \param fileOutputIdentifier String that is prepended to the name of the output files.
 * \param benchmarkInterpolator Shared pointer to a OneDimensionalInterpolator object that is used to interpolate the benchmark
 * solution to the epochs of the propagated orbit.
 */
void writePropagationResultsAndBenchmarkDifferenceToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::string& fileOutputIdentifier,
        const std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkInterpolator );

/*! Function to write the numerically propagated orbit and its difference w.r.t. to an analytical Kepler orbit to a file.
 *
 *  THIS FUNCTION OUTPUTS THE FOLLOWING FILES INTO THE src3/SimulationOutput FOLDER:
 *
 *  -  fileOutputIdentifier_numerical_states.dat Numerical integrated state history of the spacecraft
 *  -  fileOutputIdentifier_keplerian_difference.dat State difference history between numerical and analytical (= Kepler) orbits
 *  -  fileOutputIdentifier_dependent_states.dat Dependent variable history (if any) for propagation
 *
 *  Where 'fileOutputIdentifier' is a string replaced by the value of the associated input variable to this function.
 *
 * \param dynamicsSimulator Shared pointer to a SingleArcDynamicsSimulator object.
 * \param fileOutputIdentifier String that is prepended to the name of the output files.
 * \param centralBodyGravitationalParameter Value of the gravitational parameter of the central body in the simulation.
 */
void writePropagationResultsAndAnalyticalSolutionToFile(
        std::shared_ptr< SingleArcDynamicsSimulator< > >& dynamicsSimulator,
        const std::string& fileOutputIdentifier,
        const double centralBodyGravitationalParameter );

#endif // JUICEINTEGRATORANALYSISTOOLS_H
