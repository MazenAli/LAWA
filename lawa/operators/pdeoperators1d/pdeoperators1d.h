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
 
#include <lawa/operators/pdeoperators1d/convectionoperator1d.h>
#include <lawa/operators/pdeoperators1d/helmholtzoperator1d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/pdeoperator1d.h>
#include <lawa/operators/pdeoperators1d/convectionoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/pdeoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/weightedhelmholtzoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedidentityoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedlaplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedconvectionoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedconvectionoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/weightedpdeoperator1d.h>
#include <lawa/operators/pdeoperators1d/weightedpdeoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/transposedweightedpdeoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/transposedconvectionoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/transposedidentityoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/transposedlaplaceoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/transposedpdeoperator1d_pg.h>
