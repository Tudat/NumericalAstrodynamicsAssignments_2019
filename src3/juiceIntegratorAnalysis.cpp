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

#include "juiceIntegratorAnalysisTools.h"

int main()
{
    ///////////////////////////////////////////////////////////////////
    ////////  LOAD DATA FROM JSON FILES
    ///////////////////////////////////////////////////////////////////

    std::string inputJsonFile = "juiceIntegratorAnalysisSettings.json";
    boost::filesystem::path inputFilePath = getPathForJSONFile( inputJsonFile, getCurrentRootPath( ) + "JsonInput/" );
    boost::filesystem::current_path( inputFilePath.parent_path( ) );
    nlohmann::json jsonObject = getDeserializedJSON( inputJsonFile );

    std::shared_ptr< SpiceSettings > spiceSettings;
    updateFromJSON(spiceSettings, jsonObject, Keys::spice);
    loadSpiceKernels(spiceSettings);

    ///////////////////////////////////////////////////////////////////
    ////////     CREATE BODIES
    ///////////////////////////////////////////////////////////////////

    NamedBodyMap bodyMap;
    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettingsMap;
    std::string globalFrameOrigin = getValue< std::string >( jsonObject, Keys::globalFrameOrigin );
    std::string globalFrameOrientation = getValue< std::string >( jsonObject, Keys::globalFrameOrientation );
    updateBodiesFromJSON< >( jsonObject, bodyMap, bodySettingsMap, globalFrameOrigin, globalFrameOrientation,
                             spiceSettings );

    ///////////////////////////////////////////////////////////////////
    ////////     DEFINE SIMULATION SETTIONS COMMON FOR ALL RUNS
    ///////////////////////////////////////////////////////////////////

    const std::vector< std::string > bodiesToIntegrate = { "JUICE" };
    std::vector< std::string > centralBodiesPerPhase = { "Callisto", "Ganymede" };
    std::vector< double > initialTimesPerPhase = {
        getValue< double >( jsonObject, "flybyInitialTime" ),
        getValue< double >( jsonObject, "orbitInitialTime" ) };
    std::vector< double > propagationTimesPerPhase = { 8.0 * 3600.0, 24.0 * 3600.0 };

     //! STUDENT CODE TASK: define variable step-size integrator according to your student number
    RungeKuttaCoefficients::CoefficientSets coefficientSet;// =

    ///////////////////////////////////////////////////////////////////
    ////////    RUN CODE FOR QUESTION 1    ////////////////////////////
    ///////////////////////////////////////////////////////////////////

    if( getValue< bool >( jsonObject, "runQuestion1" ) )
    {
        // Step fixed sizes to use
        std::vector< double > stepSizes = { 25.0, 50.0, 100.0, 200.0 };

        // Run code for flyby and orbit phase
        for( unsigned int currentPhase = 0; currentPhase < centralBodiesPerPhase.size( ); currentPhase++)
        {
            std::string currentCentralBody = centralBodiesPerPhase.at( currentPhase );
            std::vector< std::string > centralBodies = { currentCentralBody };

            // Create initial state and time
            double currentPhaseStartTime = initialTimesPerPhase.at( currentPhase );
            double currentPhaseEndTime = currentPhaseStartTime + propagationTimesPerPhase.at( currentPhase );
            Eigen::Vector6d initialState = getBodyCartesianStateAtEpoch(
                        "JUICE", currentCentralBody, globalFrameOrientation, "NONE", currentPhaseStartTime );

            // Create accelerations
            AccelerationMap accelerationMap = getUnperturbedAccelerations(
                        currentCentralBody, bodyMap );

            // Define propagator settings
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
            propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies, accelerationMap, bodiesToIntegrate,
                        initialState, std::make_shared< PropagationTimeTerminationSettings >( currentPhaseEndTime, true ) );

            // Propagate for each required time step
            for( unsigned int stepSizeCounter = 0; stepSizeCounter < stepSizes.size( ); stepSizeCounter++ )
            {
                // Integrator settings
                std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                        getFixedStepSizeIntegratorSettings( currentPhaseStartTime, stepSizes.at( stepSizeCounter ) );

                 //! STUDENT CODE TASK: create dynamics simulator
                std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;

                // Write numerical results, and difference w.r.t. Kepler orbit, to file
                std::string fileOutputIdentifier =
                        "Q1_StepIndex" + std::to_string( stepSizeCounter ) + "_PhaseIndex" + std::to_string( currentPhase );
                writePropagationResultsAndAnalyticalSolutionToFile(
                            dynamicsSimulator, fileOutputIdentifier,
                            bodyMap.at( currentCentralBody )->getGravityFieldModel( )->getGravitationalParameter( ) );
            }
        }
    }

    ///////////////////////////////////////////////////////////////////
    ////////    RUN CODE FOR QUESTION 2 AND 3   ////////////////////
    ///////////////////////////////////////////////////////////////////

    // Check which questions to run
    std::vector< int > questionsToRun;
    if( getValue< bool >( jsonObject, "runQuestion2" ) )
    {
        questionsToRun.push_back( 2 );
    }

    if( getValue< bool >( jsonObject, "runQuestion3" ) )
    {
        // Question 3 can only be run if question 2 runs first
        if( questionsToRun.at( 0 ) != 2 )
        {
            throw std::runtime_error( "Error, cannot run Q3 without Q2" );
        }
        questionsToRun.push_back( 3 );
    }

    // Define list of integrator tolerances
    std::vector< double > integratorTolerances = { 1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6 };

    // Declare maps that will contain benchmark solutions (generated in question 2)
    std::map< double, Eigen::VectorXd > question2FlybyBenchmark;

    // Iterate over question 2 (and 3 if needed)
    for( unsigned int i = 0; i < questionsToRun.size( ); i++ )
    {
        // Retrieve current question
        unsigned int currentQuestion = questionsToRun.at( i );

        // Run only flyby for question 3
        int numberOfPhases =  currentQuestion == 2 ? centralBodiesPerPhase.size( ) : 1;

        // Iterate over flyby phase and orbit phase (for question 2 only)
        for( int currentPhase = 0; currentPhase < numberOfPhases; currentPhase++)
        {
            std::string currentCentralBody = centralBodiesPerPhase.at( currentPhase );
            std::vector< std::string > centralBodies = { currentCentralBody };

            // Set initial time (question 2)
            double currentPhaseStartTime = initialTimesPerPhase.at( currentPhase );

            // Set initial time (question 3)
            if( currentQuestion == 3 )
            {
                currentPhaseStartTime = getClosestApproachTime( question2FlybyBenchmark );
            }

            // Define final time, centeal body and initial state
            double currentPhaseEndTime = currentPhaseStartTime + propagationTimesPerPhase.at( currentPhase );
            Eigen::Vector6d initialState = getBodyCartesianStateAtEpoch(
                        "JUICE", currentCentralBody, globalFrameOrientation, "NONE", currentPhaseStartTime );

            // Create accelerations with perturbation
            AccelerationMap accelerationMap = getPerturbedAccelerations(
                        currentCentralBody, bodyMap );

            //! STUDENT CODE TASK: create list of dependent variables (if needed)
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;

            // Create object with list of dependent variables
            std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
                    std::make_shared< DependentVariableSaveSettings >( dependentVariablesList, false );

            // Define propagator settings with perturbations
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > perturbedPropagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies, accelerationMap, bodiesToIntegrate,
                        initialState, std::make_shared< PropagationTimeTerminationSettings >( currentPhaseEndTime, true ),
                        cowell, dependentVariablesToSave );

            // Create accelerations WITHOUT perturbation
            AccelerationMap accelerationMapUnperturbed = getUnperturbedAccelerations(
                        currentCentralBody, bodyMap );

            // Define propagator settings WITHOUT perturbations
            std::shared_ptr< TranslationalStatePropagatorSettings< double > > unperturbedPropagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies, accelerationMapUnperturbed, bodiesToIntegrate,
                        initialState, std::make_shared< PropagationTimeTerminationSettings >( currentPhaseEndTime, true ),
                        cowell, dependentVariablesToSave );

            // Define integrator settings for benchmark
            double benchmarkTimeStep = ( currentPhase == 0 ) ? 1.0 : 5.0;
            std::shared_ptr< IntegratorSettings< double > > benchmarkIntegratorSettings =
                    getFixedStepSizeIntegratorSettings( currentPhaseStartTime, benchmarkTimeStep );

            // Propagate benchmark dynamics
            //! STUDENT CODE TASK: create dynamics simulator
            std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;

            // Store benchmark results, if currently running question 2
            if( currentQuestion == 2  )
            {
                if( currentPhase == 0 )
                {
                    question2FlybyBenchmark = dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
                }
            }

            // Create interpolator for benchmark results
            std::shared_ptr< OneDimensionalInterpolator< double, Eigen::VectorXd > > benchmarkInterpolator =
                    createOneDimensionalInterpolator(
                        dynamicsSimulator->getEquationsOfMotionNumericalSolution( ),
                        std::make_shared< LagrangeInterpolatorSettings >( 8 ) );

            // Perform integration of perturbed dynamics with different tolerances
            for( unsigned int j = 0; j < integratorTolerances.size(); j++)
            {
                //! STUDENT CODE TASK: create integrator settings for current run;
                std::shared_ptr< IntegratorSettings< double > > integratorSettings;

                //! STUDENT CODE TASK: create dynamics simulator (perturbed dynamics)
                std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;

                // Write numerical results, and difference w.r.t. benchmark, to file
                std::string fileOutputIdentifier = "Q" + std::to_string(currentQuestion) +
                        "_ToleranceIndex" + std::to_string( j ) + "_PhaseIndex" + std::to_string( currentPhase );
                writePropagationResultsAndBenchmarkDifferenceToFile(
                            dynamicsSimulator, fileOutputIdentifier, benchmarkInterpolator );

                //! STUDENT CODE TASK: create dynamics simulator (unperturbed dynamics)
                std::shared_ptr< SingleArcDynamicsSimulator< > > unperturbedDynamicsSimulator;

                // Write numerical results, and difference w.r.t. Kepler orbit, to file
                writePropagationResultsAndAnalyticalSolutionToFile(
                            unperturbedDynamicsSimulator, fileOutputIdentifier + "_unperturbed",
                            bodyMap.at( currentCentralBody )->getGravityFieldModel( )->getGravitationalParameter( ) );
            }
        }
    }


    //    ///////////////////////////////////////////////////////////////////
    //    ////////    RUN CODE FOR QUESTION 4    ////////////////////////////
    //    ///////////////////////////////////////////////////////////////////

    if( getValue< bool >( jsonObject, "runQuestion4" ) )
    {
        std::vector< std::string > centralBodies = { "Ganymede" };
        double currentPhaseStartTime = initialTimesPerPhase.at( 1 );

        // Create accelerations
        AccelerationMap accelerationMap = getUnperturbedAccelerations(
                    "Ganymede", bodyMap );

        // Iterate over 20 individual steps, separated by 600 s
        double stepPerRun = 600.0;
        for( int i = 0; i < 20; i ++ )
        {
            double currentStartTime = currentPhaseStartTime + static_cast< double >( i ) *
                    stepPerRun;

            Eigen::Vector6d initialState = getBodyCartesianStateAtEpoch(
                        "JUICE", "Ganymede", globalFrameOrientation, "NONE", currentStartTime );
            {
                double timeStep = 300.0;

                // Define propagator settings
                std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            centralBodies, accelerationMap, bodiesToIntegrate,
                            initialState, std::make_shared< PropagationTimeTerminationSettings >(
                                currentStartTime + timeStep, true ) );


                std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                        getFixedStepSizeIntegratorSettings( currentStartTime, timeStep );

                //! STUDENT CODE TASK: create dynamics simulator that propagates for single time step
                std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;

                //! STUDENT CODE TASK: compute difference w.r.t. analytical (Keplerian) result
            }

            {
                double timeStep = 10.0;

                // Define propagator settings
                std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            centralBodies, accelerationMap, bodiesToIntegrate,
                            initialState, std::make_shared< PropagationTimeTerminationSettings >(
                                currentStartTime + timeStep, true ) );

                // Create Euler integrator settings
                std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                        std::make_shared< IntegratorSettings< > >( euler, currentStartTime, timeStep );

                //! STUDENT CODE TASK: create dynamics simulator that propagates for single time step
                std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator;

                //! STUDENT CODE TASK: compute difference w.r.t. analytical (Keplerian) result
            }
        }
    }


    // Exit main function and application
    return 0;
}
