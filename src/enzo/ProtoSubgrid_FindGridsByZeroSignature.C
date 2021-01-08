/***********************************************************************
/
/  PROTOSUBGRID CLASS (DIVIDE THIS GRID BY FINDING ZEROS IN THE SIGNATURE)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
 
int ProtoSubgrid::FindGridsByZeroSignature(int dim, int &NumberOfNewGrids,
				     int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2])
{
  /* Error check */
 
  if (dim >= GridRank) {
    ENZO_VFAIL("Passed dim(%"ISYM") > GridRank(%"ISYM")\n", dim, GridRank)
  }
 
  if (Signature[dim] == NULL) {
    ENZO_VFAIL("Signature %"ISYM" not yet computed.\n", dim)
  }

  if (GridDimension[dim] < MinimumSubgridEdge) {
    ENZO_FAIL("The subgrid is already too small!");
  }

  /* Initialize */
 
  int i = 0;
  NumberOfNewGrids = 0;
 
  /* Loop over signature. */
 
  while (i < GridDimension[dim]) {
 
    /* Look for the start of a new subgrid. */
 
    if (Signature[dim][i] != 0) {

      /* If we do not have enough zones on the right. */

      if (GridDimension[dim] - i <= 2 * MinimumSubgridEdge) {
        if (GridDimension[dim] - i <= MinimumSubgridEdge) {
          GridEnds[NumberOfNewGrids][0] = StartIndex[dim] + GridDimension[dim] - MinimumSubgridEdge;
        }
        else {
          GridEnds[NumberOfNewGrids][0] = StartIndex[dim] + i;
        }
        GridEnds[NumberOfNewGrids++][1] = StartIndex[dim] + GridDimension[dim] - 1;
        /*printf("HJ DEBUG: FindGridsByZeroSignature 0, GridEnds[NumberOfNewGrids][0]=%"ISYM", "
               "GridEnds[NumberOfNewGrids][1]=%"ISYM", StartIndex[dim]=%"ISYM", EndIndex[dim]=%"ISYM"\n",
               GridEnds[NumberOfNewGrids - 1][0], GridEnds[NumberOfNewGrids - 1][1], StartIndex[dim], EndIndex[dim]);*/
        break;
      }

      GridEnds[NumberOfNewGrids][0] = StartIndex[dim] + i;

      /* Now find the end of the subgrid. */

      i += MinimumSubgridEdge - 1;
      if (Signature[dim][i] == 0) {
        GridEnds[NumberOfNewGrids++][1] = StartIndex[dim] + i;
      }
      else {
        while (i < GridDimension[dim] && Signature[dim][i] != 0) i++;
        if (GridDimension[dim] - i < MinimumSubgridEdge) i = GridDimension[dim];
        GridEnds[NumberOfNewGrids++][1] = StartIndex[dim] + i - 1;
      }

      /*printf("HJ DEBUG: FindGridsByZeroSignature 1, GridEnds[NumberOfNewGrids][0]=%"ISYM", "
             "GridEnds[NumberOfNewGrids][1]=%"ISYM", StartIndex[dim]=%"ISYM", EndIndex[dim]=%"ISYM"\n",
             GridEnds[NumberOfNewGrids - 1][0], GridEnds[NumberOfNewGrids - 1][1], StartIndex[dim], EndIndex[dim]);*/

      if (NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS) {
        ENZO_VFAIL("PE %"ISYM" NumberOfNewGrids > MAX_NUMBER_OF_SUBGRIDS in "
                   "ProtoSubgrid_FindGridsByZeroSignature\n", MyProcessorNumber)
      }

    } // end: if (Signature[dim][i] != 0)
 
    /* Next zone in signature. */
 
    i++;
  } // end: while (i < GridDimension[dim])
 
  return SUCCESS;
}
