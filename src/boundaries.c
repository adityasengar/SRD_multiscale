#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "initialise.h" // For L_y, L_x, L_z, nearwalls

/**
 * @brief Determines which cells partly overlap with walls.
 * These cells will be filled with ghost particles during the collision steps
 * to enforce no-slip boundary conditions.
 */
void determinenearwalls()
{
  int cellx, celly, cellz;
  for (cellx = 0; cellx < L_x; cellx++) {
    for (celly = 0; celly < L_y + 1; celly++) { // L_y + 1 to include the wall at y=L_y
      for (cellz = 0; cellz < L_z; cellz++) {
        nearwalls[cellx][celly][cellz] = 0; // Default: not near a wall
        // Mark cells at y=0 or y=L_y as near walls
        if (celly == 0 || celly == L_y) {
          nearwalls[cellx][celly][cellz] = 1;
        }
      }
    }
  }
}
