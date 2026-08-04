# Define the common variables.
$patientIdsFile = "Input/PatientIds.txt"
$drugTargetProteinsFile = "Input/DrugTargets.txt"
$essentialProtFile = "Input/GBM_dependencies.txt"
$numberOfRunsPerNetwork = 100

# Define the specific variables.
$specificVariables = @(
    [Tuple]::Create(
        "Parameters/GreedyAlgNetControlParameters-MaxPath3.json",
        "Output/NetworksFromUpWithoutDown-InBetween1",
        "Output/AnalysesFromUpWithoutDown-InBetween1-MaxPath3"
    )
)

# Go over each set of specific variables.
foreach ($specificVariable in $specificVariables) {

    $parametersFile = $specificVariable.Item1
    $patientNetworksFolder = $specificVariable.Item2
    $outputFolder = $specificVariable.Item3
    
    # Go over each patient ID.
    foreach ($patientId in [System.IO.File]::ReadLines($patientIdsFile)) {
        # Go over each run index.
        for ($runIndex = 1; $runIndex -lt $numberOfRunsPerNetwork + 1; $runIndex++) {
            # Call the program.
            Programs/GreedyAlgNetControl/GreedyAlgNetControl --Mode "Cli" --Edges "$($patientNetworksFolder)/$($patientId).txt" --Targets "$($essentialProtFile)" --Preferred "$($drugTargetProteinsFile)" --Parameters "$($parametersFile)" --Output "$($outputFolder)/$($patientId)-$($runIndex).json"
        }
    }

}
