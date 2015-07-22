/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdiis/diismanager.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void DCFTSolver::compute_gradient_UHF()
{
    // Print out the header
    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf(    "\t*                           DCFT Analytic Gradients Code                          *\n");
    outfile->Printf(    "\t*                     by Alexander Sokolov and Andy Simmonett                     *\n");
    outfile->Printf(    "\t***********************************************************************************\n\n");

    // Transform the one and two-electron integrals to the MO basis and write them into the DPD file
    gradient_init();

    if (!orbital_optimized_) {
        compute_gradient_dc();
    }
    else {
        compute_gradient_odc();
    }

}

}} // Namespace
