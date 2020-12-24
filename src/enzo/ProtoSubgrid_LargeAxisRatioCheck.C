/***********************************************************************
/
/  PROTOSUBGRID CLASS (CHECK FOR LARGE AXIS RATIO AND SPLIT IF TOO LARGE)
/
/  written by: Greg Bryan
/  date:       Mar, 2014
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
 
 
int ProtoSubgrid::LargeAxisRatioCheck(int &SplitDim, int GridEnds[MAX_DIMENSION*2][2], float CriticalRatio)
{

  SplitDim = -1;

  if (GridRank == 1)
    return SUCCESS;

  /* Find longest and shortest dimension. */

  int DimLong = this->ReturnNthLongestDimension(0);
  int DimShort = this->ReturnNthLongestDimension(GridRank-1);

  /* If ratio of longest to shortest is above a critical value, then split grid
     in half along longest dimension. */

  if (float(GridDimension[DimLong])/float(GridDimension[DimShort]) > CriticalRatio) {

    int Center = (GridDimension[DimLong] - 1) / 2;
    if (Center >= MinimumSubgridEdge - 1 && EndIndex[DimLong] - StartIndex[DimLong] - Center >= MinimumSubgridEdge) {
      SplitDim = DimLong;
      //printf("GridDims are %d,\n",Center,SplitDim);
      GridEnds[SplitDim*2][0] = StartIndex[SplitDim];
      GridEnds[SplitDim*2][1] = StartIndex[SplitDim] + Center;
      GridEnds[SplitDim*2+1][0] = StartIndex[SplitDim] + Center + 1;
      GridEnds[SplitDim*2+1][1] = EndIndex[SplitDim];
    }
    /*else {
      printf("HJ DEBUG: LargeAxisRatioCheck, Center=%"ISYM", StartIndex[DimLong]=%"ISYM", EndIndex[DimLong]=%"ISYM"\n",
             Center, StartIndex[DimLong], EndIndex[DimLong]);
    }*/

    //    printf("Ori GridEnds: %d %d\n",StartIndex[SplitDim],EndIndex[SplitDim]);
    // printf("New GridEnds: %d %d %d %d\n",GridEnds[SplitDim*2][0],GridEnds[SplitDim*2][1],GridEnds[SplitDim*2+1][0],GridEnds[SplitDim*2+1][1]);
    // printf("Original grid ratio %g; new grid ratios %g\n",float(GridDimension[DimLong])/float(GridDimension[DimShort]),float(GridDimension[DimLong])/2.0/float(GridDimension[DimShort]));
    // printf("LeftEdge is %g %g %g\n",GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);
  }
 
  return SUCCESS;
}
