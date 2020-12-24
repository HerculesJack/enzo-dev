/***********************************************************************
/
/  PROTOSUBGRID CLASS (SHRINK FLAGGING CELL COVERAGE TO MINIMUM POSSIBLE)
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
 
extern "C" void FORTRAN_NAME(copy3dint)(int *source, int *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
 
int ProtoSubgrid::ShrinkToMinimumSize()
{
  int dim, i, MoveFlag = FALSE, NewGridDim[MAX_DIMENSION],
      Start[MAX_DIMENSION], End[MAX_DIMENSION];
  FLOAT CellWidth;
 
  /* First, Compute all the signatures. */
 
  for (dim = 0; dim < GridRank; dim++)
    this->ComputeSignature(dim);
 
  /* Clear index values. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim] = 0;
    End[dim] = 0;
  }
 
  /* For each dim, find first and last location of a non-zero signature. */
 
  for (dim = 0; dim < GridRank; dim++) {

    /* Check the subgrid size. */

    if (StarFeedbackSmallGridFatalError && GridDimension[dim] < MinimumSubgridEdge) {
      ENZO_FAIL("The subgrid is already too small!");
    }

    /* Look for the start. */
 
    i = 0;
    while (i < GridDimension[dim] && Signature[dim][i] == 0)
      i++;
    if (i != 0) MoveFlag = TRUE;
    Start[dim] = i + StartIndex[dim];
 
    if (i == GridDimension[dim]) {
      printf("HJ Warning: No flagged cells in ProtoSubgrid!\n");
      return NO_FLAGGED;
    }
 
    /* Look for the end. */
 
    i = GridDimension[dim] - 1;
    while (i >= 0 && Signature[dim][i] == 0)
      i--;
    End[dim] = i + StartIndex[dim];
    if (i != GridDimension[dim]-1) MoveFlag = TRUE;

    /* Check the subgrid size. */

    if ((End[dim] - Start[dim]) < MinimumSubgridEdge) {
      int m_left, m_right, n_total, n_left, n_right;
      m_left = Start[dim] - StartIndex[dim]; // # of zones remaining at the left & right
      m_right = StartIndex[dim] + GridDimension[dim] - 1 - End[dim];
      n_total = MinimumSubgridEdge + Start[dim] - End[dim] - 1;
      n_left = n_total / 2; // # of zones to extend at the left & right
      n_right = n_total - n_left;
      if (n_left <= m_left && n_right <= m_right) {
        Start[dim] -= n_left;
        End[dim] += n_right;
      }
      else if (n_left > m_left) {
        Start[dim] -= m_left;
        End[dim] += n_total - m_left;
        /*printf("HJ DEBUG: ShrinkToMinimumSize 0, m_left=%"ISYM", m_right=%"ISYM", n_left=%"ISYM", n_right=%"ISYM"\n",
               m_left, m_right, n_left, n_right);*/
      }
      else {
        Start[dim] -= n_total - m_right;
        End[dim] += m_right;
        /*printf("HJ DEBUG: ShrinkToMinimumSize 1, m_left=%"ISYM", m_right=%"ISYM", n_left=%"ISYM", n_right=%"ISYM"\n",
               m_left, m_right, n_left, n_right);*/
      }
      if (End[dim] - Start[dim] == GridDimension[dim] - 1) {
        MoveFlag = FALSE;
        /*printf("HJ DEBUG: ShrinkToMinimumSize 2, m_left=%"ISYM", m_right=%"ISYM", n_left=%"ISYM", n_right=%"ISYM"\n",
               m_left, m_right, n_left, n_right);*/
      }
    }
  }
 
  /* Move, if necessary. */
 
  if (MoveFlag) {
 
    int size = 1;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      NewGridDim[dim] = End[dim] - Start[dim] + 1;
      size *= NewGridDim[dim];
    }
 
    /*
    if (debug)
      printf("Shrinking from %"ISYM" %"ISYM" %"ISYM" to %"ISYM" %"ISYM" %"ISYM"\n",
	     GridDimension[0], GridDimension[1], GridDimension[2],
	     NewGridDim[0], NewGridDim[1], NewGridDim[2]); */
 
    /* Create new buffer and copy selected region of the GridFlaggingField. */
 
    int *TempBuffer = new int[size];
 
    FORTRAN_NAME(copy3dint)(GridFlaggingField, TempBuffer,
			 GridDimension, GridDimension+1, GridDimension+2,
			 NewGridDim, NewGridDim+1, NewGridDim+2,
			 StartIndex, StartIndex+1, StartIndex+2,
			 Start, Start+1, Start+2);
 
    /* Delete old field and put new field in it's place. */
 
    delete [] GridFlaggingField;
    GridFlaggingField = TempBuffer;
 
    /* Copy valid parts of the Signatures. */
 
    for (dim = 0; dim < GridRank; dim++)
      if (StartIndex[dim] != Start[dim])
        for (i = 0; i < NewGridDim[dim]; i++)
          Signature[dim][i] = Signature[dim][i+Start[dim]-StartIndex[dim]];
 
    /* Set new index information. */
 
    for (dim = 0; dim < GridRank; dim++) {
      CellWidth = (GridRightEdge[dim] - GridLeftEdge[dim])/
	          FLOAT(GridDimension[dim]);
      GridDimension[dim] = NewGridDim[dim];
      GridLeftEdge[dim]  = GridLeftEdge[dim] +
	FLOAT(Start[dim] - StartIndex[dim])*CellWidth;
      GridRightEdge[dim] = GridLeftEdge[dim] +
	FLOAT(GridDimension[dim])*CellWidth;
      StartIndex[dim]    = Start[dim];
      EndIndex[dim]      = End[dim];
    }
 
  }
 
  return SUCCESS;
}
