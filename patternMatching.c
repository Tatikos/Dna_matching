/**
 * @file patternMatching.c
 * @brief DNA Sequence Pattern Matching using Brute Force and Karp-Rabin algorithms
 * @author Student Name
 * @date September 2025
 * 
 * This program implements two pattern matching algorithms for DNA sequences:
 * 1. Brute Force algorithm (-bf)
 * 2. Karp-Rabin algorithm (-kr)
 * 
 * Usage: ./patternMatching -alg DNASequenceFile.txt patternFile.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/** Maximum sequence size that can be handled */
#define N 512000

/** Modulo value for Karp-Rabin hash function */
#define MOD INT_MAX

/**
 * @brief Reads a DNA sequence from a file into a character array
 * @param filename Name of the file to read from
 * @param sequence Array to store the sequence
 * @param maxSize Maximum size of the sequence
 * @return Length of the sequence read, or -1 on error
 */
int readSequence(const char* filename, char* sequence, int maxSize) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Cannot open file %s\n", filename);
        return -1;
    }
    
    int length = 0;
    int ch;
    
    while ((ch = fgetc(file)) != EOF && ch != '\n' && length < maxSize - 1) {
        if (ch == 'A' || ch == 'T' || ch == 'C' || ch == 'G' || 
            ch == 'a' || ch == 't' || ch == 'c' || ch == 'g') {
            // Convert to uppercase for consistency
            if (ch >= 'a' && ch <= 'z') {
                ch = ch - 'a' + 'A';
            }
            sequence[length++] = ch;
        }
    }
    
    sequence[length] = '\0';
    fclose(file);
    
    if (length == maxSize - 1) {
        printf("Warning: Sequence may have been truncated\n");
    }
    
    return length;
}

/**
 * @brief Implements brute force pattern matching algorithm
 * @param text The DNA sequence text to search in
 * @param pattern The pattern to search for
 * @param textLen Length of the text
 * @param patternLen Length of the pattern
 * @return Number of matches found
 */
int bruteForceSearch(const char* text, const char* pattern, int textLen, int patternLen) {
    int matches = 0;
    int i, j;
    
    // Search for pattern in text
    for (i = 0; i <= textLen - patternLen; i++) {
        j = 0;
        
        // Check if pattern matches at current position
        while (j < patternLen && text[i + j] == pattern[j]) {
            j++;
        }
        
        // If we matched the entire pattern
        if (j == patternLen) {
            matches++;
        }
    }
    
    return matches;
}

/**
 * @brief Calculates hash value for a string using rolling hash
 * @param str String to hash
 * @param len Length of the string
 * @return Hash value
 */
long long calculateHash(const char* str, int len) {
    long long hash = 0;
    long long base = 1;
    int i;
    
    // Calculate hash = (str[0]*2^(len-1) + str[1]*2^(len-2) + ... + str[len-1]*2^0) % MOD
    for (i = len - 1; i >= 0; i--) {
        hash = (hash + ((long long)str[i] * base) % MOD) % MOD;
        if (i > 0) {
            base = (base * 2) % MOD;
        }
    }
    
    return hash;
}

/**
 * @brief Recalculates hash value by removing old character and adding new one
 * @param oldChar Character being removed
 * @param oldHash Previous hash value
 * @param newChar Character being added
 * @param patternLen Length of the pattern
 * @return New hash value
 */
long long rehash(char oldChar, long long oldHash, char newChar, int patternLen) {
    long long powerOf2 = 1;
    int i;
    
    // Calculate 2^(patternLen-1)
    for (i = 0; i < patternLen - 1; i++) {
        powerOf2 = (powerOf2 * 2) % MOD;
    }
    
    // rehash(a, h, b) = ((h - a*2^(M-1)) * 2 + b) % MOD
    long long newHash = oldHash - ((long long)oldChar * powerOf2) % MOD;
    if (newHash < 0) {
        newHash += MOD;
    }
    newHash = (newHash * 2 + newChar) % MOD;
    
    return newHash;
}

/**
 * @brief Verifies if pattern actually matches at given position
 * @param text The text to check
 * @param pattern The pattern to match
 * @param pos Position in text to start checking
 * @param patternLen Length of pattern
 * @return 1 if match, 0 otherwise
 */
int verifyMatch(const char* text, const char* pattern, int pos, int patternLen) {
    int i;
    for (i = 0; i < patternLen; i++) {
        if (text[pos + i] != pattern[i]) {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Implements Karp-Rabin pattern matching algorithm
 * @param text The DNA sequence text to search in
 * @param pattern The pattern to search for
 * @param textLen Length of the text
 * @param patternLen Length of the pattern
 * @return Number of matches found
 */
int karpRabinSearch(const char* text, const char* pattern, int textLen, int patternLen) {
    if (patternLen > textLen) {
        return 0;
    }
    
    int matches = 0;
    long long patternHash = calculateHash(pattern, patternLen);
    long long textHash = calculateHash(text, patternLen);
    int i;
    
    // Check first window
    if (patternHash == textHash && verifyMatch(text, pattern, 0, patternLen)) {
        matches++;
    }
    
    // Roll through the rest of the text
    for (i = 1; i <= textLen - patternLen; i++) {
        textHash = rehash(text[i - 1], textHash, text[i + patternLen - 1], patternLen);
        
        if (patternHash == textHash && verifyMatch(text, pattern, i, patternLen)) {
            matches++;
        }
    }
    
    return matches;
}

/**
 * @brief Prints usage information
 * @param programName Name of the program
 */
void printUsage(const char* programName) {
    printf("Usage: %s -alg DNASequenceFile.txt patternFile.txt\n", programName);
    printf("Where alg can be:\n");
    printf("  -bf  : Brute Force algorithm\n");
    printf("  -kr  : Karp-Rabin algorithm\n");
}

/**
 * @brief Main function
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return 0 on success, 1 on error
 */
int main(int argc, char *argv[]) {
    // Check command line arguments
    if (argc != 4) {
        printf("Error: Invalid number of arguments\n");
        printUsage(argv[0]);
        return 1;
    }
    
    char* algorithm = argv[1];
    char* dnaFile = argv[2];
    char* patternFile = argv[3];
    
    // Validate algorithm argument
    if (strcmp(algorithm, "-bf") != 0 && strcmp(algorithm, "-kr") != 0) {
        printf("Error: Invalid algorithm. Use -bf for Brute Force or -kr for Karp-Rabin\n");
        printUsage(argv[0]);
        return 1;
    }
    
    // Allocate memory for sequences
    char* dnaSeq = (char*)malloc(N * sizeof(char));
    char* patSeq = (char*)malloc(N * sizeof(char));
    
    if (dnaSeq == NULL || patSeq == NULL) {
        printf("Error: Memory allocation failed\n");
        if (dnaSeq) free(dnaSeq);
        if (patSeq) free(patSeq);
        return 1;
    }
    
    // Read DNA sequence
    int dnaLen = readSequence(dnaFile, dnaSeq, N);
    if (dnaLen == -1) {
        printf("Error: Failed to read DNA sequence file\n");
        free(dnaSeq);
        free(patSeq);
        return 1;
    }
    
    if (dnaLen >= N - 1) {
        printf("Error: DNA sequence too large\n");
        free(dnaSeq);
        free(patSeq);
        return 1;
    }
    
    // Read pattern sequence
    int patLen = readSequence(patternFile, patSeq, N);
    if (patLen == -1) {
        printf("Error: Failed to read pattern file\n");
        free(dnaSeq);
        free(patSeq);
        return 1;
    }
    
    if (patLen >= N - 1) {
        printf("Error: Pattern sequence too large\n");
        free(dnaSeq);
        free(patSeq);
        return 1;
    }
    
    if (patLen == 0) {
        printf("Error: Empty pattern\n");
        free(dnaSeq);
        free(patSeq);
        return 1;
    }
    
    // Perform pattern matching based on selected algorithm
    int matches = 0;
    
    if (strcmp(algorithm, "-bf") == 0) {
        matches = bruteForceSearch(dnaSeq, patSeq, dnaLen, patLen);
    } else if (strcmp(algorithm, "-kr") == 0) {
        matches = karpRabinSearch(dnaSeq, patSeq, dnaLen, patLen);
    }
    
    // Output result
    printf("The pattern was found: %d times\n", matches);
    
    // Cleanup
    free(dnaSeq);
    free(patSeq);
    
    return 0;
}