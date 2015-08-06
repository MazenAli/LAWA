/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
 
#include <lawa/methods/adaptive/datastructures/alignedindexset.h>
#include <lawa/methods/adaptive/datastructures/alignedcoefficients.h>
#include <lawa/methods/adaptive/datastructures/adaptiverhs.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
//#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/methods/adaptive/datastructures/matrixoperations.h>
#include <lawa/methods/adaptive/datastructures/tensorbasis2d.h>
#include <lawa/methods/adaptive/datastructures/tensorbasis3d.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

