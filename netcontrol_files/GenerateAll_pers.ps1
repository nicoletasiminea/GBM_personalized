# Define the common variables.
$patientIdsFile = "Input/PatientIds.txt"
$mainNetworkFile = "Input/HumanOISKS.txt"
$drugTargetProteinsFile = "Input/DrugTargets.txt"
$essentialProtFile = "Input/GBM_dependencies.txt"

# Define the specific variables.
$specificVariables = @(
    [Tuple]::Create(
        "Parameters/InBetweenNetGenerationParameters-InBetween1.json",
        "Output/NetworksFromUp-InBetween1",
        "Output/NetworksFromUpWithoutDown-InBetween1",
        "Input/SetA",
        "Input/SetB",        
        "Input/SetMDown"
    )
)

# Go over each set of specific variables.
foreach ($specificVariable in $specificVariables) {

    # Get the paths.
    $parametersFile = $specificVariable.Item1
    $initialOutputFolder = $specificVariable.Item2
    $parsedOutputFolder = $specificVariable.Item3
    $firstProteinsFolder = $specificVariable.Item4
    $secondProteinsFolder = $specificVariable.Item5
    $mdownProteinsFolder = $specificVariable.Item6
    
    # Go over each patient ID.
    foreach ($patientId in [System.IO.File]::ReadLines($patientIdsFile)) {
        # Call the program.
        Programs/InBetweenNetGeneration/InBetweenNetGeneration --Mode "Cli" --MainNetwork "$($mainNetworkFile)" --DownstreamNodes "$($secondProteinsFolder)/$($patientId).txt" --UpstreamNodes "$($firstProteinsFolder)/$($patientID).txt" --Parameters "$($parametersFile)" --Output "$($initialOutputFolder)/$($patientId).txt"
    }
    
    # Go over each patient ID.
    foreach ($patientId in [System.IO.File]::ReadLines($patientIdsFile)) {
        # Define the interactions to keep.
        $interactionsToKeep = New-Object Collections.Generic.List[System.String]
        # Load the proteins encoded by genes with mutations/ down-reglated.
        $mdownProteins = [System.IO.File]::ReadLines("$($mdownProteinsFolder)/$($patientId).txt") -as [Collections.Generic.List[System.String]]
        # Go over each line in the file.
        foreach ($interaction in [System.IO.File]::ReadLines("$($initialOutputFolder)/$($patientId).txt")) {
            # Get the two proteins.
            $proteins = $interaction.Split(";")
            # Check if there are not two proteins within the interaction.
            if ($proteins.Count -ne 2) {
                # Continue.
                continue
            }
            # Check if any protein is a protein encoded by genes that exhibit mutations or are down-regulated.
            if ($mdownProteins.Contains($proteins[0])){
                continue
            }
            # Add the interaction to the list.
            $interactionsToKeep.Add($interaction)
        }
        # Create the parsed output file.
        New-Item -Path "$($parsedOutputFolder)/$($patientId).txt" -Value ([System.String]::Join("`n", $interactionsToKeep))
    }

}
