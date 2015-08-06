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
 
// Basis Operators
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivehelmholtzoperatoroptimized1d.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveconvectionoperator1d.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveidentityoperator1d.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivelaplaceoperator1d.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivepdeoperatoroptimized1d.h>
// PG Operators
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveconvectionoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveidentityoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivelaplaceoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivepdeoperator1d_pg.h>
// Weighted Operators
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveweightedpdeoperator1d.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveweightedpdeoperator1d_pg.h>
// Transposed Operators
#include <lawa/methods/adaptive/operators/pdeoperators1d/transposedadaptiveweightedpdeoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/transposedadaptiveidentityoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/transposedadaptiveconvectionoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/transposedadaptivelaplaceoperator1d_pg.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/transposedadaptivepdeoperator1d_pg.h>
