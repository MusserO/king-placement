#include <iostream>
#include <vector>
#include <set> // Still needed for orbit calculation in precomputation
#include <unordered_map>
#include <cstdint>      // For uint64_t
#include <numeric>      // For std::iota
#include <algorithm>    // For std::min, std::max, std::sort
#include <chrono>       // For timing
#include <stdexcept>    // For std::runtime_error
#include <functional>   // For std::function
#include <cmath>        // For std::ceil, std::log2 (optional fallback)
#include <fstream>      // << Added for file output
#include <iomanip>      // << Added for std::hex output formatting
#include <string>


// Compiler-specific intrinsics for finding the index of the least significant bit
#ifdef _MSC_VER
#include <intrin.h> // For _BitScanForward64
#endif

// --- Configuration ---
// Define board size N here. Change this value and recompile to test different sizes (N <= 8).
const int N = 8;
// Maximum expected Grundy value + 1 for optimized MEX calculation buffer size. Adjust if needed.
const int MAX_GRUNDY_VALUE_PLUS_1 = 16;
const std::string OUTPUT_FILE = "memo_states.bin"; // Output file for memoization table
const bool SAVE_GRUNDY_ZERO_STATES = true; // Save Grundy value 0 states to the output file
const bool SAVE_NON_ZERO_GRUNDY_STATES = false; // Save states with non-zero Grundy values to the output file
// --- End Configuration ---

// Derived constants
const int N_SQUARES = N * N;
static_assert(N > 0 && N <= 8, "Board size N must be between 1 and 8 for uint64_t bitboard.");


// Type aliases
using Bitboard = uint64_t;
// Memoization maps the *canonical* board state to its Grundy value
using MemoTable = std::unordered_map<Bitboard, int>;

// Globals
MemoTable memo; // Memoization table
Bitboard attack_masks[N_SQUARES]; // Precomputed attack masks for placing a king
Bitboard FULL_BOARD_MASK; // Bitmask representing all squares available

// --- Bit Manipulation Helpers ---

/**
 * @brief Gets the index (0..N*N-1) of the least significant bit (LSB) set in the bitboard.
 * Uses compiler intrinsics for speed if available.
 * @param bb The bitboard. Assumes bb is not 0.
 * @return The 0-based index of the LSB.
 */
inline int bit_scan_forward(Bitboard bb) {
#ifdef __GNUC__ // GCC or Clang
    return __builtin_ctzll(bb); // Count Trailing Zeros (unsigned long long)
#elif defined(_MSC_VER) // Microsoft Visual C++
    unsigned long index;
    unsigned char isNonzero = _BitScanForward64(&index, bb);
    return isNonzero ? (int)index : -1;
#else
    if (bb == 0) return -1;
    Bitboard mask = 1ULL;
    for (int i = 0; i < 64; ++i) {
        if (bb & mask) return i;
        mask <<= 1;
        if (mask == 0) break;
    }
    return -1;
#endif
}

/**
 * @brief Creates a bitmask with only the bit for square (r, c) set.
 */
inline Bitboard square_to_bitmask(int r, int c) {
    int index = r * N + c;
    if (r < 0 || r >= N || c < 0 || c >= N) {
         throw std::out_of_range("Square coordinates out of range in square_to_bitmask");
    }
    return 1ULL << index;
}

// --- Precomputation Functions ---

/**
 * @brief Precomputes the attack masks for each square.
 */
void precompute_attack_masks() {
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            int square_index = r * N + c;
            Bitboard mask = 0ULL;
            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    int nr = r + dr;
                    int nc = c + dc;
                    if (nr >= 0 && nr < N && nc >= 0 && nc < N) {
                        mask |= square_to_bitmask(nr, nc);
                    }
                }
            }
            attack_masks[square_index] = mask;
        }
    }
}


/**
 * @brief Performs all necessary precomputations.
 */
constexpr uint64_t make_full_board_mask(int n_squares) {
    return (n_squares >= 64) ? 0xFFFFFFFFFFFFFFFFULL :
           (n_squares > 0)   ? ((1ULL << n_squares) - 1) :
                               0ULL;
}

void precompute_all() {
    FULL_BOARD_MASK = make_full_board_mask(N_SQUARES);
    precompute_attack_masks();
}

// --- Bitboard Symmetry Transformations ---
// Helper to get bit index
inline int get_bit_index(int r, int c) {
    return r * N + c;
}

/**
 * @brief Applies a specific symmetry transformation (0-7) to a full bitboard.
 * Operates based on the '1' = available convention.
 */
Bitboard transform_bitboard(Bitboard bb, int type) {
    // This loop-based version is correct but potentially slow.
    // Optimizing this with bit-twiddling hacks is a major speedup area.
    Bitboard result = 0ULL;
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < N; ++c) {
            int original_index = get_bit_index(r, c);
            if ((bb >> original_index) & 1ULL) {
                int tr, tc;
                switch (type) {
                    case 0: tr = r; tc = c; break; // Identity
                    case 1: tr = c; tc = N - 1 - r; break; // Rotate 90
                    case 2: tr = N - 1 - r; tc = N - 1 - c; break; // Rotate 180
                    case 3: tr = N - 1 - c; tc = r; break; // Rotate 270
                    case 4: tr = r; tc = N - 1 - c; break; // Reflect Horiz (Vertical Axis)
                    case 5: tr = N - 1 - r; tc = c; break; // Reflect Vert (Horizontal Axis)
                    case 6: tr = c; tc = r; break; // Reflect Diag (y=x)
                    case 7: tr = N - 1 - c; tc = N - 1 - r; break; // Reflect AntiDiag
                    default: throw std::runtime_error("Invalid transformation type");
                }
                result |= (1ULL << get_bit_index(tr, tc));
            }
        }
    }
    return result;
}

/**
 * @brief Finds the canonical representation of a state (numerically smallest bitboard).
 */
Bitboard get_canonical_state(Bitboard bb) {
    if (bb == 0ULL) return 0ULL;

    Bitboard canonical_bb = bb;
    for (int i = 1; i < 8; ++i) {
        Bitboard symmetric_bb = transform_bitboard(bb, i);
        if (symmetric_bb < canonical_bb) {
            canonical_bb = symmetric_bb;
        }
    }
    return canonical_bb;
}


// --- Optimized MEX Calculation ---
/**
 * @brief Calculates the Minimum Excluded value (MEX) using a boolean vector.
 * @param seen_values A vector where index i is true if Grundy value i was seen.
 * @return The smallest non-negative integer index that is false.
 */
int calculate_mex_opt(const std::vector<bool>& seen_values) {
    for (int i = 0; i < seen_values.size(); ++i) {
        if (!seen_values[i]) {
            return i;
        }
    }
    // This indicates Grundy values exceeded the buffer size.
    std::cerr << "Warning: MEX calculation exceeded buffer size ("
              << seen_values.size() << "). Returning buffer size." << std::endl;
    return seen_values.size();
}

// --- Grundy Calculation (Recursive with Full Symmetry + Optimized MEX) ---

/**
 * @brief Recursively calculates the Grundy value for a given board state.
 * Uses memoization keyed by the canonical state and optimized MEX calculation.
 * @param board_mask Bitboard where '1' represents an available square.
 * @return The Grundy value (Nim-sum) for the state.
 */
int grundy(Bitboard board_mask) {
    // Base case: no squares available
    if (board_mask == 0) {
        return 0;
    }

    // Calculate canonical state for memoization key
    Bitboard canonical_key = get_canonical_state(board_mask);

    // Check memoization table using the canonical key
    auto it = memo.find(canonical_key);
    if (it != memo.end()) {
        return it->second; // Return cached result
    }

    // Use a boolean vector for fast MEX calculation
    std::vector<bool> seen_grundy_values(MAX_GRUNDY_VALUE_PLUS_1, false);
    bool move_found = false; // Track if any move is possible

    // Iterate through all available squares ('1' bits) in the original mask
    Bitboard remaining_mask = board_mask;
    while (remaining_mask != 0) {
        move_found = true; // A move is possible
        // Isolate the least significant bit (LSB)
        Bitboard pos_bb = remaining_mask & -remaining_mask;
        // Find the index (0..N*N-1) of this available square
        int index = bit_scan_forward(pos_bb);
         if (index < 0 || index >= N_SQUARES) {
             throw std::runtime_error("Invalid index obtained from bit_scan_forward");
         }

        // Calculate the board mask after placing a king at 'index'
        // '& ~attack_masks' removes the king and its neighbors
        Bitboard next_mask = board_mask & ~attack_masks[index];

        // Recursively find the Grundy value of the resulting state
        int next_g = grundy(next_mask);

        // Mark this Grundy value as seen (if within buffer range)
        if (next_g < MAX_GRUNDY_VALUE_PLUS_1) {
            seen_grundy_values[next_g] = true;
        } else {
            // Handle potentially large Grundy value
             std::cerr << "Warning: Calculated Grundy value " << next_g
                       << " exceeds MAX_GRUNDY_VALUE_PLUS_1 (" << MAX_GRUNDY_VALUE_PLUS_1
                       << ") for state " << std::hex << canonical_key << std::dec << std::endl; // Print key in hex
             // Option: Dynamically resize (might impact performance consistency)
             // if (next_g >= seen_grundy_values.size()) {
             //    seen_grundy_values.resize(next_g + 1, false);
             // }
             // seen_grundy_values[next_g] = true;
        }


        // Clear the LSB to move to the next available square
        remaining_mask ^= pos_bb;
    }

    // Calculate MEX using the optimized function
    int result = calculate_mex_opt(seen_grundy_values);

    // Store the result in the memoization table using the canonical key
    memo[canonical_key] = result;
    return result;
}


// --- Main Function ---

int main() {
    // Perform all necessary precomputations once
    precompute_all();

    std::cout << "Calculating Grundy value for " << N << "x" << N << " board..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    memo.clear(); // Ensure memoization table is empty before starting
    int grundy_value = grundy(FULL_BOARD_MASK); // Start calculation with the full board

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration_ms = end_time - start_time;

    // --- Output Results ---
    std::cout << "\n--- Results for " << N << "x" << N << " ---" << std::endl;
    std::cout << "The Grundy value is: " << grundy_value << std::endl;
    std::cout << "Calculation time: " << duration_ms.count() << " ms" << std::endl;
    std::cout << "Number of memoized states: " << memo.size() << std::endl;
    std::cout << "-----------------------------" << std::endl;

    std::cout << "Saving memoization table to " << OUTPUT_FILE << "..." << std::endl;
    size_t saved_states = 0;

    // --- Optimization: If only saving Grundy=0 states, don't save the value (it's always 0) ---
    bool only_save_zero = SAVE_GRUNDY_ZERO_STATES && !SAVE_NON_ZERO_GRUNDY_STATES;

    if (OUTPUT_FILE.substr(OUTPUT_FILE.find_last_of(".") + 1) == "bin") {
        std::ofstream memo_outfile(OUTPUT_FILE, std::ios::binary);
        if (!memo_outfile.is_open()) {
            std::cerr << "Error: Could not open " << OUTPUT_FILE << " for writing!" << std::endl;
        } else {
            for (const auto& pair : memo) {
                Bitboard canonical_state = pair.first;
                uint8_t grundy_byte = static_cast<uint8_t>(pair.second);

                if ((grundy_byte == 0 && !SAVE_GRUNDY_ZERO_STATES) ||
                    (grundy_byte != 0 && !SAVE_NON_ZERO_GRUNDY_STATES)) {
                    continue;
                }

                memo_outfile.write(reinterpret_cast<const char*>(&canonical_state), sizeof(Bitboard));
                if (!only_save_zero) {
                    memo_outfile.write(reinterpret_cast<const char*>(&grundy_byte), 1);
                }
                ++saved_states;
            }
            memo_outfile.close();
            std::cout << "Memoization table successfully written to " << OUTPUT_FILE << "." << std::endl;
        }
    } else {
        std::ofstream memo_outfile(OUTPUT_FILE);
        if (!memo_outfile.is_open()) {
            std::cerr << "Error: Could not open " << OUTPUT_FILE << " for writing!" << std::endl;
        } else {
            memo_outfile << "# Memoization Table: Canonical State (Hex Bitmask) -> Grundy Value\n";
            memo_outfile << "# Board Size N = " << N << "\n";
            memo_outfile << "# Bitmask convention: '1' = available square\n";
            memo_outfile << "# SAVE_GRUNDY_ZERO_STATES = " << (SAVE_GRUNDY_ZERO_STATES ? "true" : "false") << "\n";
            memo_outfile << "# SAVE_NON_ZERO_GRUNDY_STATES = " << (SAVE_NON_ZERO_GRUNDY_STATES ? "true" : "false") << "\n";
            memo_outfile << "----------------------------------------------------------\n";
            for (const auto& pair : memo) {
                Bitboard canonical_state = pair.first;
                int grundy_value = pair.second;
                if ((grundy_value == 0 && !SAVE_GRUNDY_ZERO_STATES) ||
                    (grundy_value != 0 && !SAVE_NON_ZERO_GRUNDY_STATES)) {
                    continue;
                }
                memo_outfile << "0x" << std::hex << std::setw(N_SQUARES <= 32 ? 8 : 16) << std::setfill('0') << canonical_state;
                if (!only_save_zero) {
                    memo_outfile << std::dec << " " << grundy_value;
                }
                memo_outfile << "\n";
                ++saved_states;
            }
            memo_outfile.close();
            std::cout << "Memoization table successfully written to " << OUTPUT_FILE << "." << std::endl;
        }
    }
    std::cout << "Saved Grundy=0 states: " << (SAVE_GRUNDY_ZERO_STATES ? "yes" : "no") << std::endl;
    std::cout << "Saved Grundy>0 states: " << (SAVE_NON_ZERO_GRUNDY_STATES ? "yes" : "no") << std::endl;
    if (only_save_zero) {
        std::cout << "Optimization: Grundy value is not stored in the file (always 0 for all saved states)." << std::endl;
    }
    std::cout << "Number of states saved: " << saved_states << std::endl;
    return 0;
}
