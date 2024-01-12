#include <string>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unordered_map>
#include <utility>
#include <set>
#include <vector>
#include <tuple>
#include <iostream>
#include <chrono>
#include <algorithm>

double lengthMean(char **peptides, int pep_size){
    if (pep_size == 0) return 0.0;  // consider no input situation
    double tot = 0.0;
    for (int i = 0; i <pep_size; i++){
        std::string peptide = peptides[i];
        tot += peptide.size();
    }
    return tot / pep_size;
}

double XCorr(std::vector<double> digit_spectrum, std::string peptide, double pre_mass, double resolution, 
    std::vector<std::pair<size_t, double>> mods, std::unordered_map<char, double> aa_map){
    double score = 0.0;
    const double PROTON = 1.007276;
    double inv_resolution = 1.0 / resolution; 
    std::vector<double> mod_mass_array(peptide.size());
    for (int i = 0; i < peptide.size(); i++){
        mod_mass_array[i] = aa_map[peptide[i]];
    }
    for (auto mod: mods){
        mod_mass_array[mod.first] += mod.second;
    }
    double current_b_ion = PROTON;
    for(double aa: mod_mass_array){
        if (resolution > aa && aa > -resolution) continue;  // avoid same peak counted twice for the ptm equals to the aa mass
        current_b_ion += aa;
        size_t index_b = current_b_ion * inv_resolution;  // b ion index in digit spectrum
        size_t index_y = (pre_mass + 2 * PROTON - current_b_ion) * inv_resolution; // y ion index in digit spectrum
        
        if (index_b < 0 || index_y < 0){  // consider negtive ptm mass
            return -1.0;
        }
        if (index_b < digit_spectrum.size()){
            score += digit_spectrum[index_b];
        }
        if (index_y < digit_spectrum.size()){
            score += digit_spectrum[index_y];
        }
    }
    return score;
}

std::tuple<double, size_t, double> xlmsPTMSearch(std::vector<double> digit_spectrum, double resolution, std::string peptide, std::vector<std::pair<size_t, double>> known_mods,
    std::set<int> unchange_pos, double xl_mass, int xl_site, double backbone_mass, double pre_mass, std::unordered_map<char, double> aa_map, std::unordered_map<char, std::vector<double>> unimod_map){
    
    std::tuple<double, size_t, double> res(-1.0, -1, 0.0); // score, position, mod_mass
    
    // make unchange positions from terminal aa to the known mods position

    for (auto mod: known_mods){
        int mod_idx = mod.first;
        int direction = (xl_site < mod.first) ? 1 : -1;  // could be -1/+1
        while(mod_idx >=0 && mod_idx < peptide.size()){
            unchange_pos.insert(mod_idx);
            mod_idx += direction;
        }
    }
    unchange_pos.insert(xl_site);

    for (int i = 0; i < peptide.size(); i++){
        if (unchange_pos.find(i) == unchange_pos.end()){  // aa need to be in unimod_map

            //consider regular ones
            if (unimod_map.find(peptide[i]) != unimod_map.end()){
                for (double possible_mass : unimod_map[peptide[i]]){
                    std::vector<std::pair<size_t, double>> temp_known_mods = known_mods;
                    temp_known_mods.push_back(std::make_pair(i, possible_mass));
                    double xl_site_mass = pre_mass - backbone_mass;
                    for (auto temp_known_mod : temp_known_mods){
                        xl_site_mass -= temp_known_mod.second;
                    }
                    if (xl_site_mass < xl_mass){  // make sure a valid cross-linked counterpart peptide
                        continue;
                    }
                    temp_known_mods.push_back(std::make_pair(xl_site, xl_site_mass));
                    double score = XCorr(digit_spectrum, peptide, pre_mass, resolution, temp_known_mods, aa_map);
                    if (score > std::get<0>(res)){
                        res = std::make_tuple(score, i, possible_mass);
                    }
                }
            }

            //consider n term
            if (i == 0 && unimod_map.find('[') != unimod_map.end()){
                for (double possible_mass : unimod_map['[']){
                    std::vector<std::pair<size_t, double>> temp_known_mods = known_mods;
                    temp_known_mods.push_back(std::make_pair(i, possible_mass));
                    double xl_site_mass = pre_mass - backbone_mass;
                    for (auto temp_known_mod : temp_known_mods){
                        xl_site_mass -= temp_known_mod.second;
                    }
                    if (xl_site_mass < xl_mass){  // make sure a valid cross-linked counterpart peptide
                        continue;
                    }
                    temp_known_mods.push_back(std::make_pair(xl_site, xl_site_mass));
                    double score = XCorr(digit_spectrum, peptide, pre_mass, resolution, temp_known_mods, aa_map);
                    if (score > std::get<0>(res)){
                        res = std::make_tuple(score, i, possible_mass);
                    }

                }
            
            }

            //consider c term
            if (i == peptide.size() - 1 && unimod_map.find(']') != unimod_map.end()){
                for (double possible_mass : unimod_map[']']){
                    std::vector<std::pair<size_t, double>> temp_known_mods = known_mods;
                    temp_known_mods.push_back(std::make_pair(i, possible_mass));
                    double xl_site_mass = pre_mass - backbone_mass;
                    for (auto temp_known_mod : temp_known_mods){
                        xl_site_mass -= temp_known_mod.second;
                    }
                    if (xl_site_mass < xl_mass){  // make sure a valid cross-linked counterpart peptide
                        continue;
                    }
                    temp_known_mods.push_back(std::make_pair(xl_site, xl_site_mass));
                    double score = XCorr(digit_spectrum, peptide, pre_mass, resolution, temp_known_mods, aa_map);
                    if (score > std::get<0>(res)){
                        res = std::make_tuple(score, i, possible_mass);
                    }

                }
            
            }

        }

    }
    return res;
}

std::tuple<double, size_t, double> linearPTMSearch(std::vector<double> digit_spectrum, double resolution, std::string peptide, std::vector<std::pair<size_t, double>> known_mods,
    std::set<int> unchange_pos, double backbone_mass, std::unordered_map<char, double> aa_map, std::unordered_map<char, std::vector<double>> unimod_map){
    
    std::tuple<double, size_t, double> res(-1.0, -1, 0.0); // score, position, mod_mass


    for (auto mod: known_mods){
        unchange_pos.insert(mod.first);
    }

    for (int i = 0; i < peptide.size(); i++){
        if (unchange_pos.find(i) == unchange_pos.end()){  // aa need to be in unimod_map

            //consider regular ones
            if (unimod_map.find(peptide[i]) != unimod_map.end()){
                for (double possible_mass : unimod_map[peptide[i]]){
                    std::vector<std::pair<size_t, double>> temp_known_mods = known_mods;
                    temp_known_mods.push_back(std::make_pair(i, possible_mass));
                    double mod_pep_mass = backbone_mass;
                    for (auto temp_known_mod : temp_known_mods){
                        mod_pep_mass += temp_known_mod.second;
                    }

                    double score = XCorr(digit_spectrum, peptide, mod_pep_mass, resolution, temp_known_mods, aa_map);
                    if (score > std::get<0>(res)){
                        res = std::make_tuple(score, i, possible_mass);
                    }
                }
            }

            //consider n term
            if (i == 0 && unimod_map.find('[') != unimod_map.end()){
                for (double possible_mass : unimod_map['[']){
                    std::vector<std::pair<size_t, double>> temp_known_mods = known_mods;
                    temp_known_mods.push_back(std::make_pair(i, possible_mass));
                    double mod_pep_mass = backbone_mass;
                    for (auto temp_known_mod : temp_known_mods){
                        mod_pep_mass += temp_known_mod.second;
                    }

                    double score = XCorr(digit_spectrum, peptide, mod_pep_mass, resolution, temp_known_mods, aa_map);
                    if (score > std::get<0>(res)){
                        res = std::make_tuple(score, i, possible_mass);
                    }

                }
            
            }

            //consider c term
            if (i == peptide.size() - 1 && unimod_map.find(']') != unimod_map.end()){
                for (double possible_mass : unimod_map[']']){
                    std::vector<std::pair<size_t, double>> temp_known_mods = known_mods;
                    temp_known_mods.push_back(std::make_pair(i, possible_mass));
                    double mod_pep_mass = backbone_mass;
                    for (auto temp_known_mod : temp_known_mods){
                        mod_pep_mass += temp_known_mod.second;
                    }

                    double score = XCorr(digit_spectrum, peptide, mod_pep_mass, resolution, temp_known_mods, aa_map);
                    if (score > std::get<0>(res)){
                        res = std::make_tuple(score, i, possible_mass);
                    }

                }
            
            }

        }

    }
    return res;
}



std::vector<int> charStarToPos(char *a){
    std::vector<int> res;
    // Create a copy of the input string because strtok modifies the original string
    char* a_copy = new char[strlen(a) + 1];
    strcpy(a_copy, a);

    // Split the string using strtok
    char* token = strtok(a_copy, "$");

    while (token != nullptr) {
        // Convert the token to a int using atoi function
        int number = atoi(token);

        // Print the extracted double
        res.push_back(number);

        // Get the next token
        token = strtok(nullptr, "$");
    }

    // Clean up the allocated memory
    delete[] a_copy;
    return res;
}


std::vector<double> charStarToMass(char *a){
    std::vector<double> res;
    // Create a copy of the input string because strtok modifies the original string
    char* a_copy = new char[strlen(a) + 1];
    strcpy(a_copy, a);

    // Split the string using strtok
    char* token = strtok(a_copy, "$");

    while (token != nullptr) {
        // Convert the token to a double using atof function
        double number = atof(token);

        // Print the extracted double
        res.push_back(number);

        // Get the next token
        token = strtok(nullptr, "$");
    }

    // Clean up the allocated memory
    delete[] a_copy;
    return res;
}

double calculatePeptMass(std::string str, std::unordered_map<char, double> aa_map){
    const double WATER_MASS = 18.01056;
    double mass = WATER_MASS;
    for (char &ch: str){
        mass += aa_map[ch];
    }
    return mass;
}

double partialSequenceMz(std::string partial_seq, int b_y_indicator, std::unordered_map<char, double> aa_map){
    // generate charge 1 state m/z in the spectrum
    const double PROTON = 1.007276;
    const double WATER_MASS = 18.01056;
    double mass = PROTON +  WATER_MASS;
    if (b_y_indicator == 0){
        // subtrcut water mass
        mass -= WATER_MASS;
    }
    for (char &ch: partial_seq){
        mass += aa_map[ch];
    }
    return mass;
}




extern "C" {
    void coarseScore(double *digit_spectrum, int size_spec, double *mass_1st, int *check_b_y_ion,  int *align_position, 
        char **alias_peptides, char **xl_sites, int *pattern_length, int size_pep, double pre_mass, char **aa_name, 
        double *aa_mass, int size_aa, double xl_mass, double resolution, double min_ptm, double max_ptm, double *score_ls){
        const double PROTON = 1.007276;
        std::unordered_map<char, double> aa_map;
        for (int i = 0; i < size_aa; i++){
            aa_map[aa_name[i][0]] = aa_mass[i];  //use [0] to convert char* to char
        }

        //construct vector form spectrum
        std::vector<double> vec_digit_spectrum(digit_spectrum, digit_spectrum + size_spec);
        
        // double mass_without_mod = calculatePeptMass(peptide, aa_map);
        
        for (int i = 0; i < size_pep; i++){
            std::string alias_peptide = alias_peptides[i];
            std::vector<int> vec_xl_sites = charStarToPos(xl_sites[i]);
            double alias_mass = calculatePeptMass(alias_peptide, aa_map);
            double theo_align_mz = 0.0;
            std::string partial_seq = "";
            double shift_term_mass = 0.0;
            double shift_xl_mass = 0.0;
            bool peptide_validity = false;
            int term_pos = 0;  // shift term position in alias peptide
            double current_score = 0.0;
            
            // regular linear peptide PTM, put the extra shift on both N-term & C-term
            if (vec_xl_sites.size() == 0){
                // peptide_validity = true;  // need to remove this line afterwards!!! 
                double shift_n_term_mass = 0.0;
                double shift_c_term_mass = 0.0;
                if (check_b_y_ion[i] == 0){  // b ion series situation

                    partial_seq = alias_peptide.substr(0, align_position[i]);
                    theo_align_mz = partialSequenceMz(partial_seq, check_b_y_ion[i], aa_map);
                    shift_n_term_mass = mass_1st[i] - theo_align_mz;
                    shift_c_term_mass = pre_mass - alias_mass - shift_n_term_mass;
                    if (max_ptm > shift_n_term_mass && shift_n_term_mass > min_ptm && max_ptm > shift_c_term_mass && shift_c_term_mass > min_ptm){
                        peptide_validity = true;
                    // check the validity of peptide sequence by ptm mass
                    }
                }
                else{  // y ion series situation
                    partial_seq = alias_peptide.substr(align_position[i], alias_peptide.size() - align_position[i]);
                    theo_align_mz = partialSequenceMz(partial_seq, check_b_y_ion[i], aa_map);
                    shift_c_term_mass = mass_1st[i] - theo_align_mz;
                    shift_n_term_mass = pre_mass - alias_mass - shift_c_term_mass;
                    if (max_ptm > shift_n_term_mass && shift_n_term_mass > min_ptm && max_ptm > shift_c_term_mass && shift_c_term_mass > min_ptm){
                        peptide_validity = true;
                    // check the validity of peptide sequence by ptm mass
                    }    
                }

                if (peptide_validity){  // satisfy the validity check
                    std::vector<std::pair<size_t, double>> mods;
                    mods.push_back(std::make_pair(0, shift_n_term_mass));
                    mods.push_back(std::make_pair(alias_peptide.size() - 1, shift_c_term_mass));
                    double score = XCorr(vec_digit_spectrum, alias_peptide, pre_mass, resolution, mods, aa_map);
                    if (score > current_score){
                        current_score = score;
                    }
                }
            }

            // cross-linked peptide PTM
            for (auto site : vec_xl_sites){
                peptide_validity = false;
                if (check_b_y_ion[i] == 0){  // b ion series situation
                    if (align_position[i] + pattern_length[i] <= site){
                        partial_seq = alias_peptide.substr(0, align_position[i]);
                        theo_align_mz = partialSequenceMz(partial_seq, check_b_y_ion[i], aa_map);
                        shift_term_mass = mass_1st[i] - theo_align_mz;
                        shift_xl_mass = pre_mass - alias_mass - shift_term_mass;
                        if (shift_xl_mass > xl_mass + 500.0 && max_ptm > shift_term_mass && shift_term_mass > min_ptm){
                            peptide_validity = true;
                            term_pos = 0;
                        // check the validity of peptide sequence by ptm mass and cross-linked mass
                        }

                    }else if (align_position[i] > site){
                        partial_seq = alias_peptide.substr(align_position[i], alias_peptide.size() - align_position[i]);
                        theo_align_mz = pre_mass + 2 * PROTON - partialSequenceMz(partial_seq, 1 - check_b_y_ion[i], aa_map);  // b ion
                        shift_term_mass = theo_align_mz - mass_1st[i];
                        shift_xl_mass = pre_mass - alias_mass - shift_term_mass;
                        if (shift_xl_mass > xl_mass + 500.0 && max_ptm > shift_term_mass && shift_term_mass > min_ptm){
                            peptide_validity = true;
                            term_pos = alias_peptide.size() - 1;
                        // check the validity of peptide sequence by ptm mass and cross-linked mass
                        }                    

                    }

                }
                else{  // y ion series situation

                    if (site < align_position[i] - pattern_length[i]){
                        partial_seq = alias_peptide.substr(align_position[i], alias_peptide.size() - align_position[i]);
                        theo_align_mz = partialSequenceMz(partial_seq, check_b_y_ion[i], aa_map);
                        shift_term_mass = mass_1st[i] - theo_align_mz;
                        shift_xl_mass = pre_mass - alias_mass - shift_term_mass;
                        if (shift_xl_mass > xl_mass + 500.0 && max_ptm > shift_term_mass && shift_term_mass > min_ptm){
                            peptide_validity = true;
                            term_pos = alias_peptide.size() - 1;
                        // check the validity of peptide sequence by ptm mass and cross-linked mass
                        }    
                    }else if (align_position[i] <= site){
                        partial_seq = alias_peptide.substr(0, align_position[i]);
                        theo_align_mz = pre_mass + 2 * PROTON - partialSequenceMz(partial_seq, 1 - check_b_y_ion[i], aa_map);
                        shift_term_mass = theo_align_mz - mass_1st[i];
                        shift_xl_mass = pre_mass - alias_mass - shift_term_mass;
                        if (shift_xl_mass > xl_mass + 500.0 && max_ptm > shift_term_mass && shift_term_mass > min_ptm){
                            peptide_validity = true;
                            term_pos = 0;
                        // check the validity of peptide sequence by ptm mass and cross-linked mass
                        }  

                    }
                    

                }
                if (peptide_validity){  // satisfy the validity check
                    std::vector<std::pair<size_t, double>> mods;
                    mods.push_back(std::make_pair(site, shift_xl_mass));
                    mods.push_back(std::make_pair(term_pos, shift_term_mass));
                    double score = XCorr(vec_digit_spectrum, alias_peptide, pre_mass, resolution, mods, aa_map);
                    if (score > current_score){
                        current_score = score;
                    }
                }

            }
            score_ls[i] = current_score;

            // std::vector<std::pair<size_t, double>> bold_mod;  // for no ptm score only
            // bold_mod.push_back(std::make_pair(site, pre_mass - mass_without_mod));  // only consider link site mod
            // double specific_site_current_score = XCorr(vec_digit_spectrum, peptide, pre_mass, resolution, bold_mod, aa_map);

        }

    }

    char* matchPTM(double *digit_spectrum, int size_spec, char **peptides, char **xl_sites, char **unchange_pos, int size_pep, double pre_mass, 
        char **aa_name, double *aa_mass, int size_aa, char **ptm_site, char **ptm_mass, int size_ptm, double xl_mass, double resolution, int max_unkonwn){
        //construct animo acid map
        // std::cout<<"initialize the aa_map and unimap"<<std::endl;
        // std::chrono::steady_clock sc;   // create an object of `steady_clock` class
        // auto start = sc.now();  

        std::unordered_map<char, double> aa_map;
        for (int i = 0; i < size_aa; i++){
            aa_map[aa_name[i][0]] = aa_mass[i];  //use [0] to convert char* to char
        }

        //construct vector form spectrum
        std::vector<double> vec_digit_spectrum(digit_spectrum, digit_spectrum + size_spec);
        
        //construct unimod ptm unimod map, aa:{mass2, mass2, mass3}
        std::unordered_map<char, std::vector<double>> unimod_map;  // note that n-term should be '[', c-term should be ']', input mass form 1.3$2.2$4.5
        for(int i = 0; i < size_ptm; i++){
            std::vector<double> site_ptm = charStarToMass(ptm_mass[i]);
            char site = ptm_site[i][0];
            unimod_map[site] = site_ptm;
        }
        // auto end = sc.now();       // end timer (starting & ending is done by measuring the time at the moment the process started & ended respectively)
        // auto time_span = static_cast<std::chrono::duration<double>>(end - start);   // measure time span between start & end
        // std::cout<<"Operation took: "<<time_span.count()<<" second"<<std::endl;

        std::string best_one = "";
        std::string best_under_mass_restriction = "";
        bool linearPeptide = true;
        // std::unordered_map<double, std::string> candidates;  // define candidates (peptide mass: 'sequence$score$xl_site$site:mass$site:mass')
        std::vector<std::tuple<double, std::string, std::string>> candidates;  // define candidates [(peptide mass, backbone seq, 'sequence$score$xl_site$site:mass$site:mass'), (), ()]
        double best_score_for_all = 0.0;
        double best_candidate_mass = 0.0;  // used for finding the other cross-linked peptide
        // start to compute possible ptm for each peptide sequence
        // std::cout<<"compute possible ptm"<<std::endl;
        double mean_length = lengthMean(peptides, size_pep);
        for (int i = 0; i < size_pep; i++){
            // std::cout<<"compute possible ptm1"<<std::endl;
            std::string peptide = peptides[i];
            // std::cout<<"peptide is "<<peptide<<std::endl;
            double mass_without_mod = calculatePeptMass(peptide, aa_map);
            // std::cout<<"mass is "<<mass_without_mod<<std::endl;
            

            std::vector<int> vec_xl_sites = charStarToPos(xl_sites[i]);
            std::vector<int> vec_unchange_pos = charStarToPos(unchange_pos[i]);
            
            double proportion = 0.3; // for faster implementation, require large proportion of tag found
            if (vec_unchange_pos.size() / double(peptide.size()) < proportion) continue;

            std::set<int> set_unchange_pos;
            // std::cout<<"compute possible ptm2"<<std::endl;
            for (auto pos: vec_unchange_pos){
                set_unchange_pos.insert(pos);
            }
            
            int current_site = -1;
            std::vector<std::pair<size_t, double>> current_ptms;  //initilize empty ptms
            double current_score = 0.0;  //initialize the current score with no ptm case
            // std::cout<<"initialized current score "<<current_score<<std::endl;
            // std::cout<<"compute possible ptm3"<<std::endl;


            // regular linear peptide PTM search
            if (vec_xl_sites.size() == 0){
                std::vector<std::pair<size_t, double>> bold_mod;  // for no ptm score only
                double specific_site_current_score = XCorr(vec_digit_spectrum, peptide, mass_without_mod, resolution, bold_mod, aa_map);
                // std::cout<<"bold score is "<<specific_site_current_score<<std::endl;
                
                std::vector<std::pair<size_t, double>> temp_mod;

                // add initial bold peptide information to candidates
                double total_mass = mass_without_mod;
                std::string str_candidate = peptide;

                str_candidate += "$" + std::to_string(specific_site_current_score);
                str_candidate += "$" + std::to_string(current_site);
                candidates.push_back(std::make_tuple(total_mass, peptide, str_candidate));



                int current_number_mod = 0;
                do{
                    // std::cout<<"max unkown is "<<current_number_mod<<std::endl;
                    current_number_mod ++;
                    // std::cout<<"before "<<set_unchange_pos.size()<<std::endl;
                    std::tuple<double, size_t, double> linear_result = linearPTMSearch(vec_digit_spectrum, resolution, peptide, 
                    temp_mod, set_unchange_pos, mass_without_mod, aa_map, unimod_map);
                    // std::cout<<"after "<<set_unchange_pos.size()<<std::endl;
                    // std::cout<<"linear search "<<std::get<0>(linear_result)<<" "<<std::get<1>(linear_result)<<" "<<std::get<2>(linear_result)<<std::endl;
                    // std::cout<<"compute possible ptm4"<<std::endl;
                    double penality_index = mean_length / peptide.size();
                    // std::cout<<"current peptide is "<< peptide<<" score is "<< std::get<0>(linear_result)<< "penality is "<<1.0 - std::pow(current_number_mod / double(peptide.size()), penality_index)<<" mod "<<std::get<1>(linear_result)<<" "<<std::get<2>(linear_result)<<std::endl;
                    if (std::get<0>(linear_result) * (1.0 - std::pow(current_number_mod / double(peptide.size()), penality_index)) > specific_site_current_score * 0.98){  // provide 2% cushion
                        specific_site_current_score = std::get<0>(linear_result) * (1.0 - std::pow(current_number_mod / double(peptide.size()), penality_index));
                        temp_mod.push_back(std::make_pair(std::get<1>(linear_result), std::get<2>(linear_result)));
                        
                        // add each new peptide with ptm information to candidates
                        total_mass = mass_without_mod;
                        str_candidate = peptide;
                        str_candidate += "$" + std::to_string(specific_site_current_score);
                        str_candidate += "$" + std::to_string(current_site);
                        for (auto pri: temp_mod){
                            total_mass += pri.second;
                            str_candidate += "$" + std::to_string(pri.first) + ":" + std::to_string(pri.second);
                        }
                        candidates.push_back(std::make_tuple(total_mass, peptide, str_candidate));

                    }
                    else{
                        break;
                    }
                    // std::cout<<"compute possible ptm5"<<std::endl;
                }
                while (peptide.size() - set_unchange_pos.size() - temp_mod.size() - 1 > 0 && max_unkonwn > current_number_mod);

                if (specific_site_current_score > current_score){
                    current_ptms = temp_mod;
                    current_score = specific_site_current_score;
                }
            }


            // cross-linked peptide PTM search
            for (auto site : vec_xl_sites){
                linearPeptide = false;
                // std::cout<<"site pos is "<<site<<std::endl;

                if (set_unchange_pos.find(site) == set_unchange_pos.end()){  // link site should not be changed

                    std::vector<std::pair<size_t, double>> bold_mod;  // for no ptm score only
                    bold_mod.push_back(std::make_pair(site, pre_mass - mass_without_mod));  // only consider link site mod
                    double specific_site_current_score = XCorr(vec_digit_spectrum, peptide, pre_mass, resolution, bold_mod, aa_map);
                    // std::cout<<"bold score is "<<specific_site_current_score<<std::endl;
                    
                    std::vector<std::pair<size_t, double>> temp_mod;

                    // add initial bold peptide information to candidates
                    double total_mass = mass_without_mod;
                    std::string str_candidate = peptide;

                    str_candidate += "$" + std::to_string(specific_site_current_score);
                    str_candidate += "$" + std::to_string(site);
                    candidates.push_back(std::make_tuple(total_mass, peptide, str_candidate));



                    int current_number_mod = 0;
                    do{
                        // std::cout<<"max unkown is "<<current_number_mod<<std::endl;
                        current_number_mod ++;
                        // std::cout<<"before "<<set_unchange_pos.size()<<std::endl;
                        std::tuple<double, size_t, double> linear_result = xlmsPTMSearch(vec_digit_spectrum, resolution, peptide, 
                        temp_mod, set_unchange_pos, xl_mass, site, mass_without_mod, pre_mass, aa_map, unimod_map);
                        // std::cout<<"after "<<set_unchange_pos.size()<<std::endl;
                        // std::cout<<"linear search "<<std::get<0>(linear_result)<<" "<<std::get<1>(linear_result)<<" "<<std::get<2>(linear_result)<<std::endl;
                        // std::cout<<"compute possible ptm4"<<std::endl;
                        double penality_index = mean_length / peptide.size();
                        // std::cout<<"current peptide is "<< peptide<<" score is "<< std::get<0>(linear_result)<< "penality is "<<1.0 - std::pow(current_number_mod / double(peptide.size()), penality_index)<<" mod "<<std::get<1>(linear_result)<<" "<<std::get<2>(linear_result)<<std::endl;
                        if (std::get<0>(linear_result) * (1.0 - std::pow(current_number_mod / double(peptide.size()), penality_index)) > specific_site_current_score * 0.98){  // provide 2% cushion
                            specific_site_current_score = std::get<0>(linear_result) * (1.0 - std::pow(current_number_mod / double(peptide.size()), penality_index));
                            temp_mod.push_back(std::make_pair(std::get<1>(linear_result), std::get<2>(linear_result)));
                            
                            // add each new peptide with ptm information to candidates
                            total_mass = mass_without_mod;
                            str_candidate = peptide;
                            str_candidate += "$" + std::to_string(specific_site_current_score);
                            str_candidate += "$" + std::to_string(site);
                            for (auto pri: temp_mod){
                                total_mass += pri.second;
                                str_candidate += "$" + std::to_string(pri.first) + ":" + std::to_string(pri.second);
                            }
                            candidates.push_back(std::make_tuple(total_mass, peptide, str_candidate));

                        }
                        else{
                            break;
                        }
                        // std::cout<<"compute possible ptm5"<<std::endl;
                    }
                    while (peptide.size() - set_unchange_pos.size() - temp_mod.size() - 1 > 0 && max_unkonwn > current_number_mod);

                    if (specific_site_current_score > current_score){
                        current_site = site;
                        current_ptms = temp_mod;
                        current_score = specific_site_current_score;
                    }

                }

            }
            // calculate the best one situation
            double total_mass_1 = mass_without_mod;
            std::string str_candidate_1 = peptide;
            

            str_candidate_1 += "$" + std::to_string(current_score);
            str_candidate_1 += "$" + std::to_string(current_site);
            for (auto pri: current_ptms){
                total_mass_1 += pri.second;
                str_candidate_1 += "$" + std::to_string(pri.first) + ":" + std::to_string(pri.second);
            }

            // std::cout<<"current peptide: "<< str_candidate_1 << " "<< current_score <<std::endl;
            if (current_score > best_score_for_all){
                best_score_for_all = current_score;
                best_one = str_candidate_1;
                best_candidate_mass = total_mass_1;
            }
        }
        // try to match the cross-linked peptides by mass restriction
        if (!linearPeptide){
            std::sort(candidates.begin(), candidates.end());
            
            int idx_l = 0;  // left (small mass) index
            int idx_r = candidates.size() - 1;  // right (large mass) index
            while(idx_l <= idx_r){
                if (std::get<0>(candidates[idx_l]) + std::get<0>(candidates[idx_r]) + xl_mass > pre_mass + 2 * resolution){
                    idx_r--;
                }
                else if(std::get<0>(candidates[idx_l]) + std::get<0>(candidates[idx_r]) + xl_mass < pre_mass - 2 * resolution){
                    idx_l++;
                }
                else if(std::get<1>(candidates[idx_l]) != std::get<1>(candidates[idx_r])){  // report two peptide results if they come from different backbone
                
                    best_under_mass_restriction = std::get<2>(candidates[idx_l]) + ";" + std::get<2>(candidates[idx_r]);
                    char *chars = new char[best_under_mass_restriction.size() + 1];
                    strcpy(chars, best_under_mass_restriction.c_str());
                    return chars;
                    
                }
                else{
                    idx_l++;  // or idx_r--
                }
            }
        }else{  // linear peptide mass restriction, in liear case, peptide mass has to be equal to precursor mass
            best_one = "";
            std::sort(candidates.begin(), candidates.end());
            int idx_l = 0;  // left (small mass) index

            while (idx_l < candidates.size()){
                if (std::get<0>(candidates[idx_l]) > pre_mass +  resolution){
                    break;
                }else if (std::get<0>(candidates[idx_l]) < pre_mass -  resolution){
                    idx_l ++;
                }else{
                    best_under_mass_restriction = std::get<2>(candidates[idx_l]);
                    char *chars = new char[best_under_mass_restriction.size() + 1];
                    strcpy(chars, best_under_mass_restriction.c_str());
                    return chars;
                }
            }
        }

        // report one peptide result that in the end cannot satisfy mass restriction
        char *chars = new char[best_one.size() + 1];
        strcpy(chars, best_one.c_str());
        return chars;

    }

    void freeMemory(char *chars){
        delete [] chars;
    }


    void digitizeSpectrum( double *mz, double *intensity, int size_spec, double *digit_spectrum, double resolution){
        // Important! Size cannot be zero! SEQUEST way to proceed spectra.
        const int resize = mz[size_spec - 1] / resolution + 1;

        //double * res = calloc(resize, sizeof(double));
        
        // normalize in local range
        const int num_step = 10;
        int step_unit  = resize / 10 + 1;
        int start = 0;
        int end = 0;
        int step = 0;

        do {
            double lower_bound = step * step_unit * resolution;
            double upper_bound = lower_bound + step_unit * resolution;
            start = end;
            while (end < size_spec && mz[end] <= upper_bound) {
                ++end;
            }
            // from start to end (excluded) is the range for normalization
            if (start == end) {
                continue;
            }
            double local_max_intensity = 0.0;
            for (int i = start; i < end; ++i) {
                if (intensity[i] > local_max_intensity) {
                    local_max_intensity = intensity[i];
                }
            }

            // fill processed peaks
            double normalized_factor = 50.0 / sqrt(local_max_intensity);
            for (int i = start; i < end; ++i) {
                double normalized_intensity = sqrt(intensity[i]) * normalized_factor;
                int index =  int(mz[i] / resolution);
                // take the max value
                // if (normalized_intensity > processed_peaks[index]) {
                //    processed_peaks[index] = normalized_intensity;
                // }
                // take the sum value
                digit_spectrum[index] += normalized_intensity;
            }
        } while ((++step) < num_step);
        std::vector<double> leftshifted(resize);

        // add flanking peaks by default

        for (size_t i = 0; i + 1 < resize; ++i) {
            leftshifted[i] = digit_spectrum[i + 1] * 0.5;
        }
        leftshifted[resize - 1] = 0;

        std::vector<double> rightshifted(resize);

        for (int i = resize - 1; i > 0; --i) {
            rightshifted[i] = digit_spectrum[i - 1] * 0.5;
        }
        rightshifted[0] = 0;

        for (int i = 0; i < resize; ++i) {
            digit_spectrum[i] += (leftshifted[i] + rightshifted[i]);
        }
    }


    void constructTag( const double *mz, const double *intensity, int size_spec, double da, int *start_idx, int *end_idx,
        char **char_name, int size_tag, char **aa_name, double *aa_mass, int size_aa){
        if (size_spec == 0){
            return;
        }
        // std::cout<<"here is ok 1"<<std::endl;
        std::vector<std::tuple<int, int, char*>> temp;  // store the temporary tags (start, end, name)
        double max_mass = aa_mass[0];
        double min_mass = aa_mass[0];
        for (int i = 0; i < size_aa; i++){
            if (max_mass < aa_mass[i]){
                max_mass = aa_mass[i];
            }
            if (min_mass > aa_mass[i]){
                max_mass = aa_mass[i];
            }
        }
        for (int i = 0; i < size_spec - 1; i++){
            for (int j = i + 1; j < size_spec; j++){
                if (mz[j] - mz[i] > max_mass + da){
                    break;
                }
                if (mz[j] - mz[i] + da < min_mass){
                    continue;
                }
                for( int k = 0; k < size_aa; k++){
                    if((mz[j] - da - mz[i]) <= aa_mass[k] && (mz[j] + da - mz[i]) >= aa_mass[k]){
                        std::tuple<int, int, char*> tag(i, j, aa_name[k]); 
                        // start_idx[tag_idx] = i;
                        // end_idx[tag_idx] = j;
                        // char_name[tag_idx] = aa_name[k];
                        temp.push_back(tag);

                    }
                }

            }
        }
        // std::cout<<"here is ok 2"<<std::endl;
        std::set<std::tuple<int, int>> contin_tags;  // consecutive tag1s (tag2) should not be counted twice
        // std::cout<<"temp size is "<<temp.size()<<std::endl;
        if (temp.size() <= 2) return;
        for (int i = 0; i < temp.size() - 1; i++){
            int start = std::get<0>(temp[i]);
            int mid = std::get<1>(temp[i]);
            for(int j = i + 1; j < temp.size(); j++){
                if (mid < std::get<0>(temp[j])){
                    break;
                }
                else if (mid == std::get<0>(temp[j])){
                    std::tuple<int, int> contin_tag(start, std::get<1>(temp[j])); 
                    contin_tags.insert(contin_tag);
                }

            }
        }
        // std::cout<<"here is ok 3"<<std::endl;
        int tag_idx = 0;

        for (auto i: temp){
            // std::cout<<"did i failed here?"<<std::endl;
            if(tag_idx < size_tag){
                std::tuple<int, int> flag(std::get<0>(i), std::get<1>(i));
                if ('(' != std::get<2>(i)[0]){
                    start_idx[tag_idx] = std::get<0>(i);
                    end_idx[tag_idx] = std::get<1>(i);
                    char_name[tag_idx] = std::get<2>(i);
                    tag_idx++;
                }
                else if(contin_tags.find(flag) == contin_tags.end()){
                    start_idx[tag_idx] = std::get<0>(i);
                    end_idx[tag_idx] = std::get<1>(i);
                    char_name[tag_idx] = std::get<2>(i);
                    tag_idx++;
                }
                

            }

        }
        // std::cout<<"here is ok 4"<<std::endl;

    }

    void bitap(char* text, char* pattern, int size_pattern, size_t* result, int size_res, int k){
        int idx = 0;
        int m = size_pattern;
        unsigned long *R;
        unsigned long patternMask[CHAR_MAX + 1];
        int i, d;
        if (pattern[0] == '\0') return;
        if (m > 31) return;  // depends on the architecture of different PCs, 31 is conservative for all.
        if (k >= m) return;

        R = new unsigned long[(k + 1) * sizeof *R];
        for (i = 0; i <= k; ++i)
            R[i] = ~1;

        for (i = 0; i <= CHAR_MAX; ++i)
            patternMask[i] = ~0;

        for (i = 0; i < m; ++i)
            patternMask[pattern[i]] &= ~(1UL << i);

        for (i = 0; text[i] != '\0'; ++i)
        {
            unsigned long oldRd1 = R[0];

            R[0] |= patternMask[text[i]];
            R[0] <<= 1;

            for (d = 1; d <= k; ++d)
            {
                unsigned long tmp = R[d];

                R[d] = (oldRd1 & (R[d] | patternMask[text[i]])) << 1;
                oldRd1 = tmp;
            }

            if (0 == (R[k] & (1UL << m)))
            {
                result[idx] = (i - m) + 1;
                idx++;
                if (idx >= size_res){
                    break;
                }
            }
        }

        delete [] R;
    }

}


