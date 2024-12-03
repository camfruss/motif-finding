#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_SEQ_LEN 2000
#define MAX_SEQS 100
#define ACGT_CHARS "ACGT"
const int base_ordering[20] = {0,  -1, 1,  -1, -1, -1, 2,  -1, -1, -1,
                               -1, -1, -1, -1, -1, -1, -1, -1, -1, 3};

double sumLogProbs(double a, double b) {
    if (a > b) {
        return a + log1p(exp(b - a));
    } else {
        return b + log1p(exp(a - b));
    }
}

// Reads sequences from a FASTA file
int read_fasta(const char *filename, char sequences[MAX_SEQS][MAX_SEQ_LEN]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        fprintf(stderr, "Failed to open: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    char line[MAX_SEQ_LEN];
    int seq_count = 0;
    int index = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            if (index > 0) {
                sequences[seq_count][index] = '\0';
                seq_count++;
                index = 0;
            }
        } else {
            char *pos = strchr(line, '\n');
            if (pos) *pos = '\0';  // Remove newline
            strcpy(&sequences[seq_count][index], line);
            index += strlen(line);
        }
    }
    sequences[seq_count][index] = '\0';
    fclose(file);
    return seq_count + 1;
}

// Computes the consensus motif from a set of motifs
void majority(char motif_sequences[MAX_SEQS][MAX_SEQ_LEN], int num_seqs,
              int motif_len, char *consensus) {
    int freqs[MAX_SEQ_LEN][4] = {0};
    int base_ordering[256] = {0};
    for (int i = 0; i < 4; i++) {
        base_ordering[(int)ACGT_CHARS[i]] = i;
    }

    for (int i = 0; i < num_seqs; i++) {
        for (int j = 0; j < motif_len; j++) {
            char base = motif_sequences[i][j];
            freqs[j][base_ordering[(int)base]]++;
        }
    }

    for (int i = 0; i < motif_len; i++) {
        int max_idx = 0;
        for (int j = 1; j < 4; j++) {
            if (freqs[i][j] > freqs[i][max_idx]) {
                max_idx = j;
            }
        }
        consensus[i] = ACGT_CHARS[max_idx];
    }
    consensus[motif_len] = '\0';
}

// Adds or subtracts a motif from the probability distribution
void add_motif(double pi[MAX_SEQ_LEN][4], const char *sequence, int start,
               int k, int add_subtr) {
    for (int i = 0; i < k; i++) {
        char base = sequence[start + i];
        pi[i][base_ordering[(int)base - 65]] += add_subtr;
    }
}

// Computes the probability of a motif
double motif_prob(double model[MAX_SEQ_LEN][4], const char *sequence, int start,
                  int k, double denom) {
    double pr = 0.0;
    for (int offset = 0; offset < k; offset++) {
        char base = sequence[start + offset];
        pr += log(model[offset][base_ordering[(int)base - 65]]);
    }
    pr -= (k * denom);
    return pr;
}

// Calculates the likelihood of a set of motifs given their model
double likelihood(double model[MAX_SEQ_LEN][4], int starts[MAX_SEQS],
                  int num_seqs, int k, char sequences[MAX_SEQS][MAX_SEQ_LEN],
                  double denom) {
    double pr = 0.0;

    // Sum the motif probabilities for each sequence
    for (int i = 0; i < num_seqs; i++) {
        pr += motif_prob(model, sequences[i], starts[i], k, denom);
    }

    return pr;
}

void print_pi(double pi[MAX_SEQ_LEN][4], int elements) {
    printf("-------\n");
    for (int i = 0; i < elements; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%f ", pi[i][j]);
        }
        printf("\n");
    }
    printf("-------\n");
}

void print_sequences(char sequences[MAX_SEQS][MAX_SEQ_LEN], int num_seqs) {
    printf("Sequences loaded:\n");
    for (int i = 0; i < num_seqs; i++) {
        printf("Sequence %d: %s\n", i + 1, sequences[i]);
    }
}

// The Gibbs sampling algorithm
void gibbs_sampling(char sequences[MAX_SEQS][MAX_SEQ_LEN], int num_seqs, int k,
                    double epsilon, int starts[MAX_SEQS],
                    double pi[MAX_SEQ_LEN][4]) {
    print_sequences(sequences, num_seqs);
    double small_denom = log(num_seqs - 1 + 4 * epsilon);
    double big_denom = log(num_seqs + 4 * epsilon);

    // Initialize random starts
    for (int i = 0; i < num_seqs; i++) {
        // starts[i] = rand() % (strlen(sequences[i]) - k + 1);
        starts[i] = 0;
    }

    // Initialize probability distribution
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < 4; j++) {
            pi[i][j] = epsilon;
        }
    }

    // Add motifs to the distribution
    for (int i = 0; i < num_seqs; i++) {
        add_motif(pi, sequences[i], starts[i], k, 1);
    }

    // Store the best model and likelihood
    int best_starts[MAX_SEQS];
    double best_model[MAX_SEQ_LEN][4];
    double best_likelihood = -INFINITY;

    // Perform Gibbs sampling with random sampling
    int first = 1;
    int time_since_change = 0;
    while (time_since_change < 100) {
        for (int i = 0; i < num_seqs; i++) {
            // Remove the current motif from the distribution
            add_motif(pi, sequences[i], starts[i], k, -1);

            // Calculate the probabilities for each possible start position
            double motif_distr[MAX_SEQ_LEN] = {0};
            for (int j = 0; j < (strlen(sequences[i]) - k + 1); j++) {
                motif_distr[j] =
                    motif_prob(pi, sequences[i], j, k, small_denom);
            }

            double logsum = motif_distr[0];
            for (int j = 1; j < (strlen(sequences[i]) - k + 1); j++) {
                logsum = sumLogProbs(logsum, motif_distr[j]);
            }

            // Normalize each probability distribution by subtracting the logsum
            // and exponentiating
            for (int j = 0; j < (strlen(sequences[i]) - k + 1); j++) {
                motif_distr[j] = exp(motif_distr[j] - logsum);
            }

            // Sample a new start position based on the probability distribution
            double rand_val = (double)rand() / RAND_MAX;
            double cumulative_prob = 0.0;
            int new_start = 0;
            for (int j = 0; j < (strlen(sequences[i]) - k + 1); j++) {
                cumulative_prob += motif_distr[j];
                if (rand_val <= cumulative_prob) {
                    new_start = j;
                    break;
                }
            }

            // Update the start position
            starts[i] = new_start;

            // Add the new motif back into the distribution
            add_motif(pi, sequences[i], starts[i], k, 1);
        }

        // Compute the likelihood of the new model

        double new_likelihood =
            likelihood(pi, starts, num_seqs, k, sequences, big_denom);
        if (new_likelihood > best_likelihood) {
            best_likelihood = new_likelihood;
            memcpy(best_starts, starts, num_seqs * sizeof(int));
            memcpy(best_model, pi, k * 4 * sizeof(double));
            time_since_change =
                0;  // Reset the counter since we found an improvement
        }

        time_since_change++;
    }

    // Copy the best start positions and model back to the result
    memcpy(starts, best_starts, num_seqs * sizeof(int));
    memcpy(pi, best_model, k * 4 * sizeof(double));
}

int main(int argc, char **argv) {
    srand(time(NULL));

    // Parse arguments
    char *filename = "template/motif.fasta";
    int k = 10;
    double epsilon = 0.5;

    // Read sequences
    char sequences[MAX_SEQS][MAX_SEQ_LEN];
    int num_seqs = read_fasta(filename, sequences);

    // Perform Gibbs sampling
    int starts[MAX_SEQS];
    double pi[MAX_SEQ_LEN][4];
    gibbs_sampling(sequences, num_seqs, k, epsilon, starts, pi);

    // Print results
    for (int i = 0; i < num_seqs; i++) {
        const char *sequence = sequences[i];
        int start = starts[i];
        char motif[k +
                   1];  // Buffer to hold the motif (including null terminator)

        // Extract the motif from the sequence
        strncpy(motif, &sequence[start], k);
        motif[k] = '\0';  // Null-terminate the string

        // Print the motif and its start position
        printf("Sequence %d: Start position = %d, Motif = %s\n", i, start,
               motif);
    }

    return 0;
}
