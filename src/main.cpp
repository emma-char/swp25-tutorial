// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>
#include <iostream>
#include "fastq_conversion.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include <unordered_map>

using namespace seqan3::literals;

int main(int argc, char ** argv)
{
    // Configuration
    configuration config{};

    // Parser
    sharg::parser parser{"Fastq-to-Fasta-Converter", argc, argv};

    // General information.
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    // Positional option: The FASTQ file to convert.
    parser.add_positional_option(config.fasta_input,
                                 sharg::config{.description = "The FASTQ file to convert.",
                                               .validator = sharg::input_file_validator{{"fa", "fasta"}}});
    
    
    // Open: Output FASTA file. Default: print to terminal - handled in fastq_conversion.cpp.
    parser.add_option(config.fasta_output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "The output FASTA file.",
                                    .default_message = "Print to terminal (stdout)",
                                    .validator = sharg::output_file_validator{}});

    // Flag: Verose output.
    parser.add_flag(
        config.verbose,
        sharg::config{.short_id = 'v', .long_id = "verbose", .description = "Give more detailed information."});
    

    try
    {
        parser.parse(); // Trigger command line parsing.
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        return -1;
    }


    seqan3::sequence_file_input fin_from_filename{config.fasta_input};

    for (auto & record : fin_from_filename)
    {
        seqan3::debug_stream << "ID:  " << record.id() << '\n';
        seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';

        auto minimisers = record.sequence() | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
        seqan3::debug_stream << minimisers << '\n';

        std::unordered_map< uint64_t, int > u;

        //for (auto & i : minimisers)
        //u.insert(std::make_pair<uint64_t, int>(i, ));


    }
    /*
    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
    // auto minimisers = record.sequence() | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
    // results in: [10322096095657499240, 10322096095657499142, 10322096095657499224]
    // representing the k-mers [GTAC, TCGA, GACG]
    // seqan3::debug_stream << minimisers << '\n';

    // Get hash values
    uint64_t seed = 0x8F'3F'73'B5'CF'1C'9A'DE; // The default seed from minimiser_hash
    // Use XOR on all minimiser values
    auto hash_values = minimisers | std::views::transform(
        [seed](uint64_t i)
            {
                return i ^ seed;
            });
    seqan3::debug_stream << hash_values << '\n'; // results in: [182, 216, 134]

    /*std::unordered_map< uint64_t, int > u;


    for (auto & i: minimisers){
        u.insert(std::make_pair<uint64_t,int>(minimisers, x)); 
    }
     */
            



    // convert_fastq(config); // Call fastq to fasta converter.

    //if (config.verbose) // If flag is set.
    //    std::cerr << "Conversion was a success. Congrats!\n";

    return 0;
}
