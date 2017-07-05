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

#include <lawa/methods/adaptive/algorithms/lambdatilde.h>
#include <lawa/methods/adaptive/algorithms/linearsystemsolvers.h>
#include <lawa/methods/adaptive/algorithms/localrefinement.h>
#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>
#include <lawa/methods/adaptive/algorithms/security_zone.h>
#include <lawa/methods/adaptive/algorithms/thresh.h>
#include <lawa/methods/adaptive/algorithms/weightedapply1d.h>
#include <lawa/methods/adaptive/algorithms/indexset_generation.h>
#include <lawa/methods/adaptive/algorithms/sample.h>
#include <lawa/methods/adaptive/algorithms/indexops.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/optTTcore.h>
#include <lawa/methods/adaptive/algorithms/szoneres.h>
#include <lawa/methods/adaptive/algorithms/bulk.h>
