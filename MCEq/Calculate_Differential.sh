#Move the created folder, heights.txt, cosines.txt, and energies.txt into the same folder as nu_source.h



#Sorts so that the data is grouped with height decreasing for constant coszen and energy, then energy decreasing for constant coszen 
sort -n -r -k1,1 -k3,3 -k2,2 -o Sort_Init.txt ATM_Particle_Flux_January.txt

#Adds a blank line whenever the energy changes, separating the data into coszen-energy blocks
awk -v i=3 'NR>1 && $i!=p { print "" }{ p=$i } 1' "Sort_Init.txt" > tmp && mv tmp "Sort_Init.txt"

#Calculates the height differential of the production for each coszen-energy block
awk -F "   " -v EProd=0 -v PrevEFlux=0 -v MuProd=0 -v PrevMuFlux=0 '{
  if (NF != 0) {
    EProd = ($4-PrevEFlux)/(log(10)*6.778/99*$2);
    MuProd = ($5-PrevMuFlux)/(log(10)*6.778/99*$2);
    print $1"   "$2"   "$3"   "EProd"   "MuProd >> "Flux_Differential.txt";
    PrevEFlux = $4;
    PrevMuFlux = $5;
  } else {
    EProd = 0;
    MuProd = 0;
    PrevEFlux = 0;
    PrevMuFlux = 0;
    print "" >> "Flux_Differential.txt";
  }  
}' Sort_Init.txt


# Create the subfolder if it doesn't exist
mkdir -p "czenfiles"

# Run the awk command to split the files along coszen and move it to the subfolder
# the subfiles will have indexed names
awk 'BEGIN { split_file_prefix = "coszen_" 
		split_file_extension = ".txt" 
		split_file_counter = 0 }
{
    if (NF >= 4) {
        if ($1 != previous_value) {
            close(output_file)
            split_file_counter++
            output_file = "czenfiles/" split_file_prefix split_file_counter split_file_extension
        }
        print > output_file
        previous_value = $1
    }
}
END {
    close(output_file)
}' Flux_Differential.txt

#creates files containing the cosines, energies, and heights the data was generated on
awk '{ 
    col1[$1] = 1
    col2[$2] = 1
    col3[$3] = 1
}
END {
    for (value in col1) {
        print value > "cosines.txt"
    }
    for (value in col2) {
        print value > "heights.txt"
    }
    for (value in col3) {
        print value > "energies.txt"
    }
}' Flux_Differential.txt

#sorts the files
sort -n -k1,1 -o energies.txt energies.txt
sort -n -k1,1 -o heights.txt heights.txt
sort -n -r -k1,1 -o cosines.txt cosines.txt

rm Sort_Init.txt
rm Flux_Differential.txt
