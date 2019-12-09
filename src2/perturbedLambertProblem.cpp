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


int main( )
{

    ///////////////////////////////////////////////////////////////////
    ////////  LOAD DATA FROM JSON FILES
    ///////////////////////////////////////////////////////////////////

    std::string inputJsonFile = "lambertProblemSettings.json";
    boost::filesystem::path inputFilePath_ = getPathForJSONFile( inputJsonFile, getCurrentRootPath( ) + "JsonInput/" );
    boost::filesystem::current_path( inputFilePath_.parent_path( ) );
    nlohmann::json jsonObject = getDeserializedJSON( inputJsonFile );

    // Retrieve basic simulation settings
    double departureEpoch = getValue< double >( jsonObject, "departureEpoch" );
    double arrivalEpoch = getValue< double >( jsonObject, "arrivalEpoch" );
    std::string targetBody = getValue< std::string >( jsonObject, "targetBody" );

    // Load Spice data
    std::shared_ptr< SpiceSettings > spiceSettings;
    updateFromJSON( spiceSettings, jsonObject, Keys::spice );
    loadSpiceKernels( spiceSettings );


    ///////////////////////////////////////////////////////////////////
    ////////     CREATE BODIES
    ///////////////////////////////////////////////////////////////////

    // Load body settings, and create bodies
    NamedBodyMap bodyMap;
    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettingsMap;
    std::string globalFrameOrigin = getValue< std::string >( jsonObject, Keys::globalFrameOrigin );
    std::string globalFrameOrientation = getValue< std::string >( jsonObject, Keys::globalFrameOrientation );
    updateBodiesFromJSON< >( jsonObject, bodyMap, bodySettingsMap, globalFrameOrigin, globalFrameOrientation,
                             spiceSettings );


    ///////////////////////////////////////////////////////////////////
    ////////     COMPUTE LAMBERT ARC
    ///////////////////////////////////////////////////////////////////

    std::shared_ptr< Ephemeris > lambertArcStateModel =
            getLambertProblemResult( bodyMap, targetBody, departureEpoch, arrivalEpoch );

    ///////////////////////////////////////////////////////////////////
    ////////    RUN CODE FOR QUESTION 1
    ///////////////////////////////////////////////////////////////////

    if( getValue< bool >( jsonObject, "runQuestion1" ) )
    {
        // Numerically propagate trajectory without perturbations
        propagateTrajectory(
                    departureEpoch, arrivalEpoch, jsonObject, bodyMap, lambertArcStateModel, "Q1", false );
    }

    ///////////////////////////////////////////////////////////////////
    ////////    RUN CODE FOR QUESTION 2
    ///////////////////////////////////////////////////////////////////

    if( getValue< bool >( jsonObject, "runQuestion2" ) )
    {
        // Iterate over 3 required cases for 'buffer time' (no buffer time; 1 hour buffer time; 2 day buffer time), and run
        // numerical integration WITH perturbations
        for( int i = 0; i < 3; i++ )
        {
            //! STUDENT CODE TASK: set buffer time for different runs (cases i, ii and iii in question 2, for i = 0, 1 and 2)
            double currentBufferTime; // = ...
            double departureEpochWithBuffer = departureEpoch + currentBufferTime;
            double arrivalEpochWithBuffer = arrivalEpoch - currentBufferTime;

            // Propagate  trajectory from beginning to end (with modified initial and final states)
            propagateTrajectory(
                        departureEpochWithBuffer, arrivalEpochWithBuffer, jsonObject, bodyMap, lambertArcStateModel,
                        "Q2a" + std::to_string( i ), true  );


            //! STUDENT CODE TASK: propagate trajectory forward/backward from midpoint (for case ii only, see question 2b)
            if( i == 1 )
            {
                //! Propagate forward from midpoint
                //propagateTrajectory( ... );

                //! Propagate backward from midpoint
                //propagateTrajectory( ... );

            }
        }
    }

    ///////////////////////////////////////////////////////////////////
    ////////    RUN CODE FOR QUESTION 3
    ///////////////////////////////////////////////////////////////////

    if( getValue< bool >( jsonObject, "runQuestion3" ) )
    {
        // Set full start and end times.
        double currentBufferTime = 2.0 * 86400.0;
        double departureEpochWithBuffer = departureEpoch + currentBufferTime;
        double arrivalEpochWithBuffer = arrivalEpoch - currentBufferTime;

        // Set arc length
        int numberOfArcs = 10;
        double arcLength = ( arrivalEpochWithBuffer - departureEpochWithBuffer ) / static_cast< double >( numberOfArcs );

        // Initialize Delta V variables
        double totalDeltaV = 0.0;

        // Variable used to store the required Delta V at the END of each arc. HINT: use this in your algorithm!
        Eigen::Vector3d previousArcFinalDeltaV = Eigen::Vector3d::Constant( TUDAT_NAN );

        // Iterate over all arcs
        for( int arcIndex = 0; arcIndex < numberOfArcs; arcIndex ++ )
        {
            // Set initial and final times for current arc
            double currentArcInitialTime = departureEpochWithBuffer + static_cast< double >( arcIndex ) * arcLength;
            double currentArcFinalTime = departureEpochWithBuffer + static_cast< double >( arcIndex + 1 ) * arcLength;

            ///////////////////////////////////////////////////////////////////
            ////////    RUN CODE FOR QUESTION 3a
            ///////////////////////////////////////////////////////////////////

            //! STUDENT CODE TASK: propagate trajectory for current arc (with corresponding call to propagateTrajectory function)
            //propagateTrajectory( ... );

            ///////////////////////////////////////////////////////////////////
            ////////    RUN CODE FOR QUESTION 3b
            ///////////////////////////////////////////////////////////////////

            //! Perhaps you could also leave this to the students (again, makes them look into the tools.cpp file), while giving them a hint
            std::shared_ptr< SingleArcVariationalEquationsSolver< > > variationalEquationsSolver = propagateVariationalEquations(
                        currentArcInitialTime, currentArcFinalTime, jsonObject, bodyMap, lambertArcStateModel );
            std::map< double, Eigen::MatrixXd > currentArcStateTransitionMatrixHistory =
                    variationalEquationsSolver->getNumericalVariationalEquationsSolution( ).at( 0 );

            // Get final state transition matrrix
            Eigen::MatrixXd finalStateTransitionMatrix =
                    currentArcStateTransitionMatrixHistory.rbegin( )->second;

            // Retrieve numerical and Lambert arc history
            std::map< double, Eigen::VectorXd > numericalStateHistory =
                    variationalEquationsSolver->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::VectorXd > lambertArcStateHistory =
                    getLambertArcHistory( lambertArcStateModel, numericalStateHistory );

            // Compute final state deviation (between numerically propagated state and Lambert solution)
            Eigen::VectorXd finalStateDeviation =
                    numericalStateHistory.rbegin( )->second - lambertArcStateHistory.rbegin( )->second;

            //! STUDENT CODE TASK: compute required velocity change at beginning of arc to meet required final state
            Eigen::Vector3d initialVelocityChange; // = ...

            // Set the correction to the initial state
            Eigen::Vector6d initialStateCorrection = Eigen::Vector6d::Zero( );
            initialStateCorrection.segment( 3, 3 ) = initialVelocityChange;

            //! STUDENT CODE TASK: compute Delta V arc arc start and end; increment total mission Delta V
            totalDeltaV; // = totalDeltaV + ...

            // Get initial state for corrected arc
            Eigen::Vector6d correctedInitialState =
                    lambertArcStateModel->getCartesianState( currentArcInitialTime ) + initialStateCorrection;

            // Set integrator settings
            std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings;
            updateFromJSON( integratorSettings, jsonObject, Keys::integrator );
            integratorSettings->initialTime_ = currentArcInitialTime;

            // Repropagate with corrected initial state
            {
                // Set propagation settings for corrected arc
                std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                        getPerturbedPropagatorSettings( bodyMap, jsonObject, correctedInitialState, currentArcFinalTime );

                // Propagate forward perturbed arc and save results
                std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator = std::make_shared< SingleArcDynamicsSimulator< > >(
                            bodyMap, integratorSettings, propagatorSettings );
                writePropagationResultsToFile( dynamicsSimulator, lambertArcStateModel, "Q3_arc_" + std::to_string( arcIndex ) + "_corrected" );
            }

        }
    }

    ///////////////////////////////////////////////////////////////////
    ////////    RUN CODE FOR QUESTION 4
    ///////////////////////////////////////////////////////////////////

    if( getValue< bool >( jsonObject, "runQuestion4" ) )
    {
        // Set full start and end times.
        double currentBufferTime = 2.0 * 86400.0;
        double departureEpochWithBuffer = departureEpoch + currentBufferTime;
        double arrivalEpochWithBuffer = arrivalEpoch - currentBufferTime;

        // Set arc length
        int numberOfArcs = 10;
        double arcLength = ( arrivalEpochWithBuffer - departureEpochWithBuffer ) / static_cast< double >( numberOfArcs );


        std::vector< int > consideredArcIndices = { 0, 4 };
        for( unsigned int i = 0; i < consideredArcIndices.size( ); i++ )
        {
            ////////////////////////////////////////////////////////////////////////
            ////////    GET PROPAGATION AND INTEGRATION SETTINGS FOR CURRENT ARC
            ////////////////////////////////////////////////////////////////////////

            int arcIndex = consideredArcIndices.at( i );

            // Compute start and end time for current arc
            double currentArcInitialTime = departureEpochWithBuffer + static_cast< double >( arcIndex ) * arcLength;
            double currentArcFinalTime = departureEpochWithBuffer + static_cast< double >( arcIndex + 1 ) * arcLength;

            // Get propagator settings for perturbed forward and backward arcs
            Eigen::Vector6d arcInitialState = lambertArcStateModel->getCartesianState( currentArcInitialTime );
            std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings =
                    getPerturbedPropagatorSettings( bodyMap, jsonObject, arcInitialState, currentArcFinalTime );

            // Set integrator settings
            std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings;
            updateFromJSON( integratorSettings, jsonObject, Keys::integrator );
            integratorSettings->initialTime_ = currentArcInitialTime;

            ////////////////////////////////////////////////////////////////////////
            ////////    PROPAGATE NOMINAL TRAJECTORY AND VARIATIONAL EQUATIONS
            ////////////////////////////////////////////////////////////////////////

            SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
                        bodyMap, integratorSettings, propagatorSettings,
                        getSensitityParameterSet( propagatorSettings, bodyMap ),
                        true, nullptr, false, true, false );

            // Retrieve state transition matrix history (Phi(t,t0) for list of values of t)
            std::map< double, Eigen::MatrixXd > stateTransitionResult =
                    variationalEquationsSimulator.getNumericalVariationalEquationsSolution( ).at( 0 );

            // Retrieve propagated state history (x(t) for list of values of t)
            std::map< double, Eigen::VectorXd > nominalIntegrationResult =
                    variationalEquationsSimulator.getDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );

            // Retrieve nominal initial state value (e.g. initial state with Delta x_0 = 0)
            Eigen::Vector6d originalInitialState = nominalIntegrationResult.begin( )->second;

            ////////////////////////////////////////////////////////////////////////
            ////////    START ANALYSIS ALGORIHTM FOR QUESTION 4
            ////////////////////////////////////////////////////////////////////////        

            // This vector will hold the maximum permitted initial state perturbations for which the linearization is valid (for the
            // current arc. The vector is initialized to 0, and each of its 6 entries is computed in the 6 iterations of the
            // coming for loop (that runs over the iteration variable 'entry')
            Eigen::Vector6d permittedPerturbations = Eigen::Vector6d::Zero( );

            // Iterate over all initial state entries
            for( unsigned int entry = 0; entry < 6; entry++ )
            {

                //! STUDENT CODE TASK: Define (iterative) algorithm to compute current entry of 'permittedPerturbations'
                //! General structure: define an initial state perturbation (perturbedInitialState variable),
                //! compute epsilon_x (see assignment), and iterate your algorithm until convergence.

                //while( ... ) OR for ( ... ) OR your custom iteration scheme
                {
                    //! STUDENT CODE TASK define initial state perturbation for current iteration.
                    Eigen::Vector6d initialStatePerturbation; // = ...

                    // Reset propagator settings with perturbed initial state
                    Eigen::Vector6d perturbedInitialState = arcInitialState + initialStatePerturbation;
                    propagatorSettings->resetInitialStates( perturbedInitialState );

                    // Create simulation object and propagate dynamics with perturbed initial state
                    SingleArcDynamicsSimulator< > dynamicsSimulator(
                                bodyMap, integratorSettings, propagatorSettings, true, false, false );

                    // Retrieve state histroy computed directly from perturbed initial state
                    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                }


            }

            // Print permitted initial state perturbations computed by algorihm
            std::cout<<"Permitted perturbations: "<<arcIndex<<" "<<permittedPerturbations.transpose( )<<std::endl<<std::endl;
        }
    }

    return 0;
}
