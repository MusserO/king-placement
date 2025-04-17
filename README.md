# King Placement Game

This project implements a solver and web interface for the **King Placement Game**. 

## King Placement Game rules:
- The game is played on a standard `8 x 8` chessboard where two players take turns placing kings. 
- A king can **not** be placed on a square attacked by another king already on the board.
- As in chess, a king attacks all adjacent squares, including diagonals (the 8 squares around it).
- The player who **cannot make a legal move** on their turn **loses**.

> **Note:** You can play this game on a physical chessboard using pawns to represent kings. 8 kings per player is always enough.

## Features

- **Grundy Value Calculation**: Efficiently computes Grundy values (also known as nim-values) for all board states using memoization and symmetry optimizations.
- **Precomputed Data**: Stores calculated game states in `memo_states.bin` for faster future access.
- **Web Interface**: Play the game against AI or with a friend in local two-player mode through a web-based UI.
- **AI Player**: Choose between an easy AI that plays randomly and a hard AI that uses the precomputed game states to play optimally.

## Project Structure
 ├── solver.cpp # C++ implementation of the solver  
 ├── index.html # Web-based UI for the game  
 ├── memo_states.bin # Precomputed memoization table  

## How to Use

### 1. Compile and Run the Solver
Compile the C++ solver to precompute Grundy values and generate `memo_states.bin`.

```bash
g++ -std=c++17 -O3 solver.cpp -o solver
./solver
```

This will output the Grundy values and save the memoization table in a file. Takes around two minutes to complete.

### 2. Run the Web Interface
Open index.html in a browser to play the game against the AI. The AI uses the precomputed memo_states.bin for optimal moves.

**Try it online:** [mussero.github.io/king-placement/](https://mussero.github.io/king-placement/)

## Attribution
The `wK.svg` and `bK.svg` files used in the web interface are made by Colin M.L. Burnett and are licensed under the [GPLv2+](https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt) license.

## License
This project is licensed under the MIT License. See LICENSE for details.