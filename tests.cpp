
// auto time_function = [](auto&& func) -> double {
//     // collect min, max, average times 
//     auto start = std::chrono::high_resolution_clock::now();
//     std::invoke(std::forward<decltype(func)>(func));
//     auto end = std::chrono::high_resolution_clock::now();
//     return std::chrono::duration<double, std::milli>(end - start).count(); // Duration in milliseconds
// };

/* time vs. target_len (data type) */
// void PWM_TYPE()
// {
//     std::vector<int> motif_lengths { 8, 10, 12, 14, 16 };
//     int num_sequences { 50 };
//     int sequence_length { 1'000 };

//     for (const auto& ml : motif_lengths) {
//         auto motif_lengths { ml };
//         Data data { motif_lengths, num_sequences, sequence_length };

//         Serial<float> run { data };
//         run.find_motifs(ml, 1);

//         Serial<double> run { data };
//         run.find_motifs(ml, 1);

//         Serial<long double> run { data };
//         run.find_motifs(ml, 1);
//     }
// }

/* time vs. # sequences (len_sequence) */

/* time vs. pseudocount (# sequences) */

/* time vs. len_sequence (target_len) */

/* #iterations vs. correctness ( #seqs) */
