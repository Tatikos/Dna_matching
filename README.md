/**
 * @mainpage DNA Pattern Matching - Assignment 1
 * 
 * @section intro_sec Introduction
 * 
 * This program implements two different algorithms for pattern matching in DNA sequences:
 * - Brute Force algorithm
 * - Karp-Rabin algorithm
 * 
 * @section compile_sec Compilation
 * 
 * To compile the program, use:
 * ```
 * gcc -o patternMatching patternMatching.c
 * ```
 * 
 * @section usage_sec Usage
 * 
 * The program is executed from the command line with the following syntax:
 * ```
 * ./patternMatching -alg DNASequenceFile.txt patternFile.txt
 * ```
 * 
 * Where:
 * - `alg` can be either `-bf` for Brute Force or `-kr` for Karp-Rabin algorithm
 * - `DNASequenceFile.txt` contains the DNA sequence to search in
 * - `patternFile.txt` contains the pattern to search for
 * 
 * @section examples_sec Examples
 * 
 * ```
 * ./patternMatching -bf dnaSequence.txt patSequence.txt
 * ./patternMatching -kr swinefluDNA.txt normal-dna.txt
 * ```
 * 
 * @section algorithms_sec Algorithms
 * 
 * @subsection bf_sec Brute Force Algorithm
 * 
 * The brute force algorithm performs an exhaustive search by comparing the pattern
 * against every possible position in the text. For each position, it compares
 * character by character until either a mismatch is found or the entire pattern
 * is matched.
 * 
 * Time Complexity: O(n*m) where n is the text length and m is the pattern length.
 * 
 * @subsection kr_sec Karp-Rabin Algorithm
 * 
 * The Karp-Rabin algorithm uses a rolling hash function to efficiently compare
 * the pattern with substrings of the text. It consists of two phases:
 * 
 * 1. Preprocessing: Calculate hash value of the pattern
 * 2. Execution: Calculate hash values of text substrings and compare with pattern hash
 * 
 * The hash function used is:
 * hash(s) = (s[0]*2^(m-1) + s[1]*2^(m-2) + ... + s[m-1]*2^0) mod INT_MAX
 * 
 * For rolling hash calculation:
 * rehash(old_char, old_hash, new_char) = ((old_hash - old_char*2^(m-1)) * 2 + new_char) mod INT_MAX
 * 
 * Time Complexity: Average O(n+m), Worst case O(n*m)
 * 
 * @section files_sec File Format
 * 
 * Input files should contain DNA sequences on a single line, consisting only of
 * the characters A, T, C, G (case insensitive). The program automatically converts
 * lowercase letters to uppercase.
 * 
 * @section limitations_sec Limitations
 * 
 * - Maximum sequence length: 512,000 characters (defined by N constant)
 * - Only supports standard DNA nucleotides (A, T, C, G)
 * - Input files should have sequences on a single line
 * 
 * @section author_sec Author Information
 * 
 * - Course: EPL232 - Programming Techniques and Tools
 * - Department: Computer Science, University of Cyprus
 * - Assignment: 1 - DNA Sequence Processing and Applications
 * - Semester: Fall 2025
 * 
 * @section notes_sec Implementation Notes
 * 
 * - The program uses dynamic memory allocation for sequences
 * - Proper error handling for file operations and memory allocation
 * - Hash collision handling in Karp-Rabin through character-by-character verification
 * - Modular design with separate functions for each algorithm
 * - Comprehensive documentation using Doxygen format
 */