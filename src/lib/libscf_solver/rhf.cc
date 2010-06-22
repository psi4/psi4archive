/*
 *  rhf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/10/08.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include "rhf.h"
#include <psi4-dec.h>

#define TIME_SCF
#define _DEBUG

using namespace psi;
using namespace std;

namespace psi { namespace scf {

static void sort_rows_based_on_energies(SimpleMatrix* C, double *energies, int *order_mapping)
{
    unsigned int i, j;
    int itemp;
    double dtemp;
    int length = C->rows();

    // Populate order_mapping with original ordering.
    for (i=0; i< length; ++i)
        order_mapping[i] = i;

    // Sort using Quicksort algorithm
    for (i=0; i<length; ++i) {
        for (j = i+1; j<length; ++j) {
            if (energies[i] > energies[j]) {
                C->swap_rows(i, j);
    
                dtemp = energies[i];
                energies[i] = energies[j];
                energies[j] = dtemp;
    
                itemp = order_mapping[i];
                order_mapping[i] = order_mapping[j];
                order_mapping[j] = itemp;
            }
        }
    }
}

RHF::RHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : HF(options, psio, chkpt)
{
    common_init();
}

RHF::~RHF()
{
    if (pk_)
        delete[] pk_;
}

void RHF::common_init()
{
    Drms_ = 0.0;

    // Allocate matrix memory
    F_    = SharedMatrix(factory_.create_matrix("F"));
    D_    = SharedMatrix(factory_.create_matrix("D"));
    Dold_ = SharedMatrix(factory_.create_matrix("D old"));
    G_    = SharedMatrix(factory_.create_matrix("G"));
    J_    = SharedMatrix(factory_.create_matrix("J"));
    K_    = SharedMatrix(factory_.create_matrix("K"));

    if (scf_type_ == "L_DF") {
        Lref_ = SharedMatrix(factory_.create_matrix("Lref"));
        L_    = SharedMatrix(factory_.create_matrix("L"));
    }

    int nao = chkpt_->rd_nao();
    chkpt_->wt_nmo(nao);

    // PK super matrix for fast G
    pk_ = NULL;

    // Print DIIS status
    fprintf(outfile, "  SCF Algorithm Type is %s.\n", scf_type_.c_str());
    fprintf(outfile, "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");

    fflush(outfile);
    // Allocate memory for PK matrix
    if (scf_type_ == "PK")
        allocate_PK();
}

double RHF::compute_energy()
{
    //fprintf(outfile,"  Print = %d\n",print_);
    //print_ = options_.get_int("PRINT");
    bool converged = false, diis_iter = false;
    iteration_ = 0;

    // Do the initial work to get the iterations started.
    //form_multipole_integrals();  // handled by HF class
    timer_on("Core Hamiltonian");
    form_H(); //Core Hamiltonian
    timer_off("Core Hamiltonian");
    timer_on("Overlap Matrix");
    form_Shalf(); //Shalf Matrix
    timer_off("Overlap Matrix");
    //Form initial MO guess by user specified method
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them, unless the user disagrees.
    timer_on("Initial Guess");
    load_or_compute_initial_C();
    timer_off("Initial Guess");

    if (print_>3) {
        S_->print(outfile);
        Shalf_->print(outfile);
        if (canonical_X_)
            X_->print(outfile);
        H_->print(outfile);
    }
    if (print_>2) {
        fprintf(outfile,"  Initial Guesses:\n");
        C_->print(outfile);
        D_->print(outfile);
    }

    if (scf_type_ == "PK")
        form_PK();
    else if (scf_type_ == "CD" || scf_type_ == "1C_CD")
        form_CD();
    else if (scf_type_ == "DF")
        form_B();
    else if (scf_type_ == "L_DF") {
        form_A();
        I_ = block_matrix(basisset_->molecule()->natom(),doccpi_[0]);
        form_domain_bookkeeping();
    }

    fprintf(outfile, "                                  Total Energy            Delta E              Density RMS\n\n");
    fflush(outfile);

    // SCF iterations
    do {
        iteration_++;

        Dold_->copy(D_);  // Save previous density
        Eold_ = E_;       // Save previous energy

        //form_G_from_J_and_K(1.0);
        //D_->print(outfile);

        timer_on("Form G");
        if (scf_type_ == "PK")
            form_G_from_PK();
        else if (scf_type_ == "DIRECT")
            form_G_from_direct_integrals();
        else if (scf_type_ == "DF"||scf_type_ == "CD"||scf_type_ =="1C_CD")
           form_G_from_RI();
        else if (scf_type_ == "L_DF")
           form_G_from_RI_local_K();
        else if(scf_type_ == "OUT_OF_CORE")
           form_G();

        timer_off("Form G");

        if (print_>3) {
            J_->print(outfile);

            K_->print(outfile);

            G_->print(outfile);
        }

        form_F();
        if (print_>3) {
            F_->print(outfile);
        }
        if (diis_enabled_)
            save_fock();

        E_ = compute_E();

        timer_on("DIIS");
        if (diis_enabled_ == true && iteration_ >= num_diis_vectors_) {
            diis();
            diis_iter = true;
        } else {
            diis_iter = false;
        }
        timer_off("DIIS");

        if (print_>4 && diis_iter) {
            fprintf(outfile,"  After DIIS:\n");
            F_->print(outfile);
        }
        fprintf(outfile, "  @RHF iteration %3d energy: %20.14f    %20.14f %20.14f %s\n", iteration_, E_, E_ - Eold_, Drms_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);

        timer_on("Diagonalize H");
        form_C();
        timer_off("Diagonalize H");
        form_D();

        if (print_>2) {
            C_->print(outfile);
            D_->print(outfile);
        }

        converged = test_convergency();
    } while (!converged && iteration_ < maxiter_ ); 

    //Free the heavies pronto!
    if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD")
        free_B();
    else if (scf_type_ == "L_DF") {
        free_A();
        free_block(I_);
        free_domain_bookkeeping();
    }
    
    if (converged) {
        fprintf(outfile, "\n  Energy converged.\n");
        fprintf(outfile, "\n  @RHF Final Energy: %20.14f", E_);
        if (perturb_h_) {
            fprintf(outfile, " with %f perturbation", lambda_);
        }
        fprintf(outfile, "\n");
        save_information();
        if (options_.get_bool("DUAL_BASIS")) 
            save_dual_basis_projection();
        if (options_.get_str("SAPT") != "FALSE") //not a bool because it has types
            save_sapt_info();
    } else {
        fprintf(outfile, "\n  Failed to converged.\n");
        E_ = 0.0;
    }
    //if (save_grid_) {
        // DOWN FOR MAINTENANCE
        //fprintf(outfile,"\n  Saving Cartesian Grid\n");
        //save_RHF_grid(options_, basisset_, D_, C_);
    //}

    // Compute the final dipole.
    compute_multipole();

    //fprintf(outfile,"\nComputation Completed\n");
    fflush(outfile);
    return E_;
}

bool RHF::load_or_compute_initial_C()
{
    bool ret = false;
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them.
    string prefix(chkpt_->build_keyword(const_cast<char*>("MO coefficients")));
    //What does the user want?
    //Options will be:
    // ""-Either READ or CORE (try READ first)
    // "READ"-try to read MOs from checkpoint (restart style)
    // "BASIS2"-try to read MOs from checkpoint after basis cast up (in INPUT) NOT WORKING!
    // "DUAL_BASIS"-real the results of the DB computation from File 100, Temporary hack
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Denisties 
    string guess_type = options_.get_str("GUESS");
    if ((guess_type == "" || guess_type == "READ"||guess_type == "BASIS2")&&chkpt_->exist(const_cast<char*>(prefix.c_str()))) {
         //Try for existing MOs already. Deuces of loaded spades
        if (print_)
            fprintf(outfile, "  SCF Guess: Reading previous MOs.\n\n");

        // Read MOs from checkpoint and set C_ to them
        double **vectors = chkpt_->rd_scf();
        C_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        form_D();

        // Read SCF energy from checkpoint file.
        E_ = chkpt_->rd_escf();
        ret = true;
    } else if (guess_type == "DUAL_BASIS") {
         //Try for dual basis MOs, 
        if (print_)
            fprintf(outfile, "  SCF Guess: Dual-Basis. Reading from File 100.\n\n");

        double** C2 = block_matrix(basisset_->nbf(),nalpha_); 

        psio_->open(PSIF_SCF_DB_MOS,PSIO_OPEN_OLD);
        psio_address next_PSIF = PSIO_ZERO;
        psio_->read(PSIF_SCF_DB_MOS,"DB MO Integrals",(char *) &(C2[0][0]),sizeof(double)*basisset_->nbf()*nalpha_,next_PSIF,&next_PSIF);
        next_PSIF = PSIO_ZERO;
        psio_->read(PSIF_SCF_DB_MOS,"DB SCF Energy",(char *) &(E_),sizeof(double),next_PSIF,&next_PSIF);
        psio_->close(PSIF_SCF_DB_MOS,1);

        //fprintf(outfile, "nalpha_ = %d\n",nalpha_);
        //print_mat(C2,basisset_->nbf(),nalpha_,outfile);

        C_->zero();
        for (int m = 0; m<basisset_->nbf(); m++)
            for (int i = 0; i<nalpha_; i++)
                C_->set(0,m,i,C2[m][i]);
        
        free_block(C2); 
        
        //C_->print(outfile);
    
        memset((void*) doccpi_, '\0', factory_.nirreps()*sizeof(int));
        memset((void*) soccpi_, '\0', factory_.nirreps()*sizeof(int));
        memset((void*) nalphapi_, '\0', factory_.nirreps()*sizeof(int));
        memset((void*) nbetapi_, '\0', factory_.nirreps()*sizeof(int));

        doccpi_[0] = nalpha_;
        nalphapi_[0] = nalpha_;
        nbetapi_[0] = nalpha_;

        // Read MOs from checkpoint and set C_ to them
        form_D();

        //D_->print(outfile);

        ret = true;
        
    } else if (guess_type == "CORE" || guess_type == "") {
        //CORE is an old Psi standby, so we'll play this as spades
        if (print_)
            fprintf(outfile, "  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");
        F_->copy(H_); //Try the core Hamiltonian as the Fock Matrix
        form_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();
    } else if (guess_type == "SAD") {
        if (print_)
            fprintf(outfile, "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF.\n");
        //Superposition of Atomic Density (will be preferred when we can figure it out)
        compute_SAD_guess();
    } else if (guess_type == "GWH") {
        //Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        if (print_)
            fprintf(outfile, "  SCF Guess: Generalized Wolfsberg-Helmholtz.\n\n");
        F_->zero(); //Try F_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        S_->print(outfile);
        int *opi = S_->rowspi();
        int nirreps = S_->nirreps();
        for (h=0; h<nirreps; ++h) {
            for (i=0; i<opi[h]; ++i) {
                for (j=0; j<opi[h]; ++j) {
                    F_->set(h,i,j,0.5*S_->get(h,i,j)*(H_->get(h,i,i)+H_->get(h,j,j)));
                }
            }
        }
        form_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();

    } else if (guess_type == "READ" || guess_type == "BASIS2") {
        throw std::invalid_argument("Checkpoint MOs requested, but do not exist!");
    }
    if (print_)
        fprintf(outfile, "\n  Initial RHF energy: %20.14f\n\n", E_);
    fflush(outfile);
    return ret;
}

void RHF::save_dual_basis_projection()
{
    if (print_)
        fprintf(outfile,"\n  Computing dual basis set projection from %s to %s.\n  Results will be stored in File 100.\n",options_.get_str("BASIS").c_str(),options_.get_str("DUAL_BASIS_SCF").c_str()); 
    shared_ptr<BasisSet> dual_basis =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DUAL_BASIS_SCF"));
    SharedMatrix C_2 = dualBasisProjection(C_,doccpi_[0],basisset_,dual_basis);
    //C_->print(outfile);
    if(print_>3)
        C_2->print(outfile);
    //Write stuff to chkpt here
    double** C2 = C_2->to_block_matrix();
    //print_mat(C2,dual_basis->nbf(),doccpi_[0],outfile);
    //fflush(outfile);

    psio_->open(PSIF_SCF_DB_MOS,PSIO_OPEN_NEW);
    psio_address next_PSIF = PSIO_ZERO;
    psio_->write(PSIF_SCF_DB_MOS,"DB MO Integrals",(char *) &(C2[0][0]),sizeof(double)*dual_basis->nbf()*doccpi_[0],next_PSIF,&next_PSIF);
    next_PSIF = PSIO_ZERO;
    psio_->write(PSIF_SCF_DB_MOS,"DB SCF Energy",(char *) &(E_),sizeof(double),next_PSIF,&next_PSIF);
    psio_->close(PSIF_SCF_DB_MOS,1);
    free_block(C2); 
}

void RHF::compute_multipole()
{
    // Begin dipole
    double dex, dey, dez, dx, dy, dz;
    // Convert blocked density to a full block
    SimpleMatrix D(D_->to_simple_matrix());
    D.set_name("Full Block Density");

    dex = D.vector_dot(Dipole_[0]) * 2.0;
    dey = D.vector_dot(Dipole_[1]) * 2.0;
    dez = D.vector_dot(Dipole_[2]) * 2.0;

    dx = dex + nuclear_dipole_contribution_[0];
    dy = dey + nuclear_dipole_contribution_[1];
    dz = dez + nuclear_dipole_contribution_[2];

    double d;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    // End dipole

    // Begin quadrupole
    double qexx, qexy, qexz, qeyy, qeyz, qezz;
    double mexx, mexy, mexz, meyy, meyz, mezz;
    double texx, texy, texz, teyy, teyz, tezz;

    mexx = D.vector_dot(Quadrupole_[0]) * 2.0;
    mexy = D.vector_dot(Quadrupole_[1]) * 2.0;
    mexz = D.vector_dot(Quadrupole_[2]) * 2.0;
    meyy = D.vector_dot(Quadrupole_[3]) * 2.0;
    meyz = D.vector_dot(Quadrupole_[4]) * 2.0;
    mezz = D.vector_dot(Quadrupole_[5]) * 2.0;

    texx = mexx + nuclear_quadrupole_contribution_[0];
    texy = mexy + nuclear_quadrupole_contribution_[1];
    texz = mexz + nuclear_quadrupole_contribution_[2];
    teyy = meyy + nuclear_quadrupole_contribution_[3];
    teyz = meyz + nuclear_quadrupole_contribution_[4];
    tezz = mezz + nuclear_quadrupole_contribution_[5];

    qexx = texx - (teyy+tezz)/2.0;
    qeyy = teyy - (texx+tezz)/2.0;
    qezz = tezz - (texx+teyy)/2.0;
    qexy = 1.5 * texy;
    qexz = 1.5 * texz;
    qeyz = 1.5 * teyz;

    SimpleVector evals(3);
    SimpleMatrix evecs(3, 3), temp(3, 3);

    temp.set(0, 0, qexx);
    temp.set(1, 1, qeyy);
    temp.set(2, 2, qezz);
    temp.set(0, 1, qexy);
    temp.set(0, 2, qexz);
    temp.set(1, 2, qeyz);
    temp.set(1, 0, qexy);
    temp.set(2, 0, qexz);
    temp.set(2, 1, qeyz);
    temp.diagonalize(&evecs, &evals);
    // End Quadrupole

    fprintf(outfile, "\n  Electric dipole (a.u.):\n");
    fprintf(outfile, "\t%15s\t%15s\t%15s", "X", "Y", "Z");
    fprintf(outfile, "\n    Nuclear part:\n");
    fprintf(outfile, "\t%15.10f\t%15.10f\t%15.10f\n", nuclear_dipole_contribution_[0], nuclear_dipole_contribution_[1], nuclear_dipole_contribution_[2]);

    fprintf(outfile, "\n    Electronic part:\n");
    fprintf(outfile, "\t%15.10f\t%15.10f\t%15.10f\n", dex, dey, dez);

    fprintf(outfile, "\n    Dipole moments:\n");
    fprintf(outfile, "\t%15.10f\t%15.10f\t%15.10f\n", dx, dy, dz);

    fprintf(outfile, "\n    Total dipole: %15.10f a.u.  %15.10f Debye\n", d, d*_dipmom_au2debye);
    fprintf(outfile, "    Conversion: 1.0 a.u. = %15.10f Debye\n", _dipmom_au2debye);

    if (print_ > 1) {
        SimpleMatrix *C = C_->to_simple_matrix();
        // Transform dipole integrals to MO basis
        Dipole_[0]->transform(C);
        Dipole_[1]->transform(C);
        Dipole_[2]->transform(C);

        fprintf(outfile, "\n    Orbital contributions to dipole (a.u.)\n");
        fprintf(outfile, "\t%3s%15s  %15s  %15s\n", "MO", "X", "Y", "Z");
        double totx=0.0, toty=0.0, totz=0.0;
        for (int i=0; i<nalpha_; ++i) {
            double x_contrib = 0.0, y_contrib = 0.0, z_contrib = 0.0;
            x_contrib += Dipole_[0]->get(i, i) * 2.0;
            y_contrib += Dipole_[1]->get(i, i) * 2.0;
            z_contrib += Dipole_[2]->get(i, i) * 2.0;
            totx += x_contrib;
            toty += y_contrib;
            totz += z_contrib;
            fprintf(outfile, "\t%3d%15.10f  %15.10f  %15.10f\n", i+1,
                    x_contrib, y_contrib, z_contrib);
        }

        // fprintf(outfile, "\n  Electric quadrupole (a.u.):");
        // fprintf(outfile, "\n    Nuclear part:\n");
        // fprintf(outfile, "\txx=%15.10f\txy=%15.10f\txz=%15.10f\n", nuclear_quadrupole_contribution_[0], nuclear_quadrupole_contribution_[1], nuclear_quadrupole_contribution_[2]);
        // fprintf(outfile, "\tyy=%15.10f\tyz=%15.10f\tzz=%15.10f\n", nuclear_quadrupole_contribution_[3], nuclear_quadrupole_contribution_[4], nuclear_quadrupole_contribution_[5]);
        // fprintf(outfile, "\n    Electronic part (cross terms do not match oeprop):\n");
        // fprintf(outfile, "\txx=%15.10f\txy=%15.10f\txz=%15.10f\n", qexx, qexy, qexz);
        // fprintf(outfile, "\tyy=%15.10f\tyz=%15.10f\tzz=%15.10f\n", qeyy, qeyz, qezz);
        // fprintf(outfile, "\n    Principal values (a.u.) and axis:\n");
        // fprintf(outfile, "\tQ1=%15.10f\tV1=(%15.10f %15.10f %15.10f)\n", evals.get(0), evecs.get(0, 0), evecs.get(1, 0), evecs.get(2, 0));
        // fprintf(outfile, "\tQ2=%15.10f\tV2=(%15.10f %15.10f %15.10f)\n", evals.get(1), evecs.get(0, 1), evecs.get(1, 1), evecs.get(2, 1));
        // fprintf(outfile, "\tQ3=%15.10f\tV3=(%15.10f %15.10f %15.10f)\n", evals.get(2), evecs.get(0, 2), evecs.get(1, 2), evecs.get(2, 2));

        // Compute orbital extents
        SimpleMatrix orbital_extents("Orbital Extents", C->cols(), 4);
        for (int i=0; i<C->cols(); ++i) {
            double sumx=0.0, sumy=0.0, sumz=0.0;
            for (int k=0; k<C->rows(); ++k) {
                for (int l=0; l<C->rows(); ++l) {
                    double tmp = C->get(k, i) * C->get(l, i);
                    sumx += Quadrupole_[0]->get(k, l) * tmp;
                    sumy += Quadrupole_[3]->get(k, l) * tmp;
                    sumz += Quadrupole_[5]->get(k, l) * tmp;
                }
            }

            orbital_extents.set(i, 0, fabs(sumx));
            orbital_extents.set(i, 1, fabs(sumy));
            orbital_extents.set(i, 2, fabs(sumz));
            orbital_extents.set(i, 3, fabs(sumx + sumy + sumz));
        }

        // Sort orbital extent rows based on orbital energies
        double *orbital_energies = orbital_energies_->to_block_vector();
        int *order_mapping = new int[C->cols()];
        sort_rows_based_on_energies(&orbital_extents, orbital_energies, order_mapping);

        fprintf(outfile, "\n  Orbital extents (a.u.):\n");
        fprintf(outfile, "\t%3s%15s  %15s  %15s  %15s\n", "MO", "<x^2>", "<y^2>", "<z^2>", "<r^2>");

        for (int i=0; i<orbital_extents.rows(); ++i) {
            fprintf(outfile, "\t%3d%15.10f  %15.10f  %15.10f  %15.10f\n", i+1, 
                    orbital_extents.get(i, 0),orbital_extents.get(i, 1),orbital_extents.get(i, 2),orbital_extents.get(i, 3));
        }
        delete[] orbital_energies;
        delete[] order_mapping;
        delete C;
    }
}

/** CURRENTLY DOWN FOR MAINTENANCE
void RHF::save_RHF_grid(Options& opts, shared_ptr<BasisSet> basis, SharedMatrix D, SharedMatrix C)
{
    SharedProperties prop = SharedProperties(Properties::constructProperties(basis));
    int nmo_c;
    int* mo_inds_c;
    if (opts.get_int("N_CARTESIAN_MOS") != 0)
    {
        nmo_c = opts.get_int("N_CARTESIAN_MOS");
        mo_inds_c = init_int_array(nmo_c);
        for (int k = 0; k<nmo_c; k++)
        {
            mo_inds_c[k] = opts["CARTESIAN_MO_INDICES"][k].to_integer();
        }
        prop->setToComputeMOs(true,mo_inds_c,nmo_c);
    }
    double * extents = getCartesianGridExtents(opts, basis->molecule());
    double xmin = extents[0];
    double xmax = extents[1];
    double ymin = extents[2];
    double ymax = extents[3];
    double zmin = extents[4];
    double zmax = extents[5];
    free(extents);

    int * npoints = getCartesianGridResolution(opts);
    int npointsx = npoints[0];
    int npointsy = npoints[1];
    int npointsz = npoints[2];
    free(npoints);
    FILE * grid_file = fopen((opts.get_str("CARTESIAN_FILENAME")).c_str(),"w");
    shared_ptr<Molecule> mol = basisset_->molecule();
    fprintf(grid_file,"%d\n",mol->natom());
    fprintf(grid_file,"x,y,z,Z\n");
    for (int i=0; i<mol->natom(); i++)
        fprintf(grid_file,"%14.10f,%14.10f,%14.10f,%d\n",mol->x(i),mol->y(i), mol->z(i), mol->Z(i));
    fprintf(grid_file,"%10d,%10d,%10d\n",npointsx,npointsy,npointsz);
    fprintf(grid_file,"x,y,z,\\rho");
    if (nmo_c>0){
        for (int m=0; m<nmo_c; m++)
            fprintf(grid_file,",MO%d",mo_inds_c[m]);
    }
    fprintf(grid_file,"\n");
    double x,y,z,val;
    for (int i = 0; i<npointsx; i++) {
        if (npointsx == 1)
            x = (xmax+xmin)/2.0;
        else
            x = i*(xmax-xmin)/(npointsx-1.0)+xmin;
        for (int j = 0; j<npointsy; j++) {
            if (npointsy == 1)
                y = (ymax+ymin)/2.0;
            else
                y = j*(ymax-ymin)/(npointsy-1.0)+ymin;
            for (int k = 0; k<npointsz; k++) {
                if (npointsz == 1)
                    z = (zmax+zmin)/2.0;
                else
                    z = k*(zmax-zmin)/(npointsz-1.0)+zmin;
                Vector3 v(x,y,z);
                prop->computeProperties(v,D,C); //Needs to be updated as props changes
                val = prop->getDensity();
                fprintf(grid_file, "%14.10f,%14.10f,%14.10f,%14.10f",x,y,z,val);
                if (nmo_c>0) {
                    for (int m = 0; m<nmo_c; m++)
                        fprintf(grid_file,",%14.10f",prop->getMO(m));
                }
                fprintf(grid_file, "\n");
            }
        }
    }
    fclose(grid_file);

    if (opts.get_int("N_CARTESIAN_MOS") != 0)
        free(mo_inds_c);

}
int* RHF::getCartesianGridResolution(Options &opts)
{
    int *npoints = init_int_array(3);
    //Use individual resolutions
    if (opts["CARTESIAN_RESOLUTION_X"].has_changed()) {
     npoints[0] = opts.get_int("CARTESIAN_RESOLUTION_X");
     npoints[1] = opts.get_int("CARTESIAN_RESOLUTION_Y");
     npoints[2] = opts.get_int("CARTESIAN_RESOLUTION_Z");
    } else { //Use global or default resolutions
     npoints[0] = opts.get_int("CARTESIAN_RESOLUTION");
     npoints[1] = npoints[0];
     npoints[2] = npoints[0];
    }
    return npoints;
}
double* RHF::getCartesianGridExtents(Options &opts, shared_ptr<Molecule> mol)
{
    double *ext = init_array(6);
    if (opts["CARTESIAN_EXTENTS"].has_changed()) {
    //Use defined extents
    for (int i=0; i<6; i++)
        ext[i] = opts["CARTESIAN_EXTENTS"][i].to_double();
    } else {
        //Use molecule size plus overage
        double xmin = mol->x(0);
        double xmax = mol->x(0);
        double ymin = mol->y(0);
        double ymax = mol->y(0);
        double zmin = mol->z(0);
        double zmax = mol->z(0);
        for (int i=1;i<mol->natom(); i++) {
            if (xmin>mol->x(i))
                xmin = mol->x(i);
            if (ymin>mol->y(i))
                ymin = mol->y(i);
            if (zmin>mol->z(i))
                zmin = mol->z(i);
            if (xmax<mol->x(i))
                xmax = mol->x(i);
            if (ymax<mol->y(i))
                ymax = mol->y(i);
            if (zmax<mol->z(i))
                zmax = mol->z(i);
        }
        double overage = opts.get_double("CARTESIAN_OVERAGE");
        ext[0] = xmin - overage;
        ext[1] = xmax + overage;
        ext[2] = ymin - overage;
        ext[3] = ymax + overage;
        ext[4] = zmin - overage;
        ext[5] = zmax + overage;
    }
    return ext;
}**/

void RHF::save_information()
{
    // Print the final docc vector
    char **temp2 = chkpt_->rd_irr_labs();
    int nso = chkpt_->rd_nso();

    fprintf(outfile, "\n  Final occupation vector = (");
    for (int h=0; h<factory_.nirreps(); ++h) {
        fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
    }
    fprintf(outfile, ")\n");

    // Needed for a couple of places.
    SharedMatrix eigvector(factory_.create_matrix());
    SharedVector eigvalues(factory_.create_vector());

    F_->diagonalize(eigvector, eigvalues);

    int print_mos = false;
    print_mos = options_.get_bool("PRINT_MOS");
    if (print_mos) {
        fprintf(outfile, "\n  Molecular orbitals:\n");

        C_->eivprint(eigvalues);
    }

    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<eigvalues->nirreps(); ++h) {
        for (int i=0; i<eigvalues->dimpi()[h]; ++i)
            pairs.push_back(make_pair(eigvalues->get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());
    int ndocc = 0;
    for (int i=0; i<eigvalues->nirreps(); ++i)
        ndocc += doccpi_[i];

    fprintf(outfile, "\n  Orbital energies (a.u.):\n    Doubly occupied orbitals\n      ");
    for (int i=1; i<=ndocc; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first, temp2[pairs[i-1].second]);
        if (i % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "\n    Unoccupied orbitals\n      ");
    for (int i=ndocc+1; i<=nso; ++i) {
        fprintf(outfile, "%12.6f %3s  ", pairs[i-1].first, temp2[pairs[i-1].second]);
        if ((i-ndocc) % 4 == 0)
            fprintf(outfile, "\n      ");
    }
    fprintf(outfile, "\n");

    for (int i=0; i<eigvalues->nirreps(); ++i)
        free(temp2[i]);
    free(temp2);

    int *vec = new int[eigvalues->nirreps()];
    for (int i=0; i<eigvalues->nirreps(); ++i)
        vec[i] = 0;

    chkpt_->wt_nmo(nso);
    chkpt_->wt_ref(0);        // Only RHF right now
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_orbspi(eigvalues->dimpi());
    chkpt_->wt_openpi(vec);
    chkpt_->wt_phase_check(0);

    // Figure out frozen core orbitals
    int nfzc = chkpt_->rd_nfzc();
    int nfzv = chkpt_->rd_nfzv();
    int *frzcpi = compute_fcpi(nfzc, eigvalues);
    int *frzvpi = compute_fvpi(nfzv, eigvalues);
    chkpt_->wt_frzcpi(frzcpi);
    chkpt_->wt_frzvpi(frzvpi);
    delete[](frzcpi);
    delete[](frzvpi);

    // Save the Fock matrix
    // Need to recompute the Fock matrix as F_ is modified during the SCF interation
    form_F();
    double *ftmp = F_->to_lower_triangle();
    chkpt_->wt_fock(ftmp);
    delete[](ftmp);

    // This code currently only handles RHF
    chkpt_->wt_iopen(0);

    // Write eigenvectors and eigenvalue to checkpoint
    double *values = eigvalues->to_block_vector();
    chkpt_->wt_evals(values);
    delete[] values;
    double **vectors = C_->to_block_matrix();
    chkpt_->wt_scf(vectors);
    free_block(vectors);
}

void RHF::save_fock()
{
    static bool initialized_diis_manager = false;
    if (initialized_diis_manager == false) {
        diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, F_.get());
        diis_manager_->set_vector_size(1, DIISEntry::Matrix, F_.get());
        initialized_diis_manager = true;
    }

    // Determine error matrix for this Fock
    SharedMatrix FDS(factory_.create_matrix()), DS(factory_.create_matrix());
    SharedMatrix SDF(factory_.create_matrix()), DF(factory_.create_matrix());

    // FDS = F_ * D_ * S_;
    DS->gemm(false, false, 1.0, D_, S_, 0.0);
    FDS->gemm(false, false, 1.0, F_, DS, 0.0);
    // SDF = S_ * D_ * F_;
    DF->gemm(false, false, 1.0, D_, F_, 0.0);
    SDF->gemm(false, false, 1.0, S_, DF, 0.0);

    Matrix FDSmSDF;
    FDSmSDF.copy(FDS);
    FDSmSDF.subtract(SDF);

    // Orthonormalize the error matrix
    FDSmSDF.transform(Shalf_);
    //FDSmSDF.print(outfile);

    diis_manager_->add_entry(2, &FDSmSDF, F_.get());
}

void RHF::diis()
{
    diis_manager_->extrapolate(1, F_.get());
}

bool RHF::test_convergency()
{
    // energy difference
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix D_rms;
    D_rms.copy(D_);
    D_rms.subtract(Dold_);
    Drms_ = D_rms.rms();

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void RHF::allocate_PK()
{
    // The size of the pk matrix is determined in HF::form_indexing
    // Allocate memory for the PK matrix (using a vector)
    if (pk_size_ < (memory_ / sizeof(double))) {
        pk_ = new double[pk_size_];

        if (pk_ == NULL) {
            fprintf(outfile, "  Insufficient free system memory for in-core PK implementation.\n");
            fprintf(outfile, "  Switching to out-of-core algorithm.\n");
            scf_type_ = "OUT_OF_CORE";
        } else {
             // Zero out PK
            memset(pk_, 0, pk_size_*sizeof(double));
             fprintf(outfile, "  Allocated %lu elements (%lu pairs) for PK. (%5.2f MiB)\n\n", (unsigned long)pk_size_, (unsigned long)pk_pairs_, pk_size_ * 8.0 / 1048576.0);
        }
    }
    else {
        fprintf(outfile, "  Insufficient memory for in-core PK implementation.\n");
        fprintf(outfile, "  Would need %lu elements of double memory. (%5f MiB)\n", (unsigned long)pk_size_, pk_size_ * sizeof(double) / 1048576.0);
        fprintf(outfile, "  Switching to out-of-core algorithm.\n");
        scf_type_ = "OUT_OF_CORE";
    }
}
void RHF::form_F()
{
    F_->copy(H_);
    F_->add(G_);

#ifdef _DEBUG
    if (debug_) {
        F_->print(outfile);
    }
#endif
}
void RHF::form_C()
{
    if (!canonical_X_) {
        Matrix eigvec;
        Vector eigval;
        factory_.create_matrix(eigvec);
        factory_.create_vector(eigval);

        F_->transform(Shalf_);
        F_->diagonalize(eigvec, eigval);

        // Save the orbital energies
        orbital_energies_->copy(eigval);

        C_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

        find_occupation(eigval);

        // Save C to checkpoint file.
        //double **vectors = C_->to_block_matrix();
        //chkpt_->wt_scf(vectors);
        //free_block(vectors);

#ifdef _DEBUG
        if (debug_) {
            C_->eivprint(eigval);
        }
#endif
    } else {

        C_->zero();
        Vector eigval;
        factory_.create_vector(eigval);

        for (int h = 0; h<C_->nirreps(); h++) {

            int norbs = nsopi_[h];
            int nmos = nmopi_[h];

            //fprintf(outfile,"  Norbs = %d, Nmos = %d\n",norbs,nmos); fflush(outfile);
            
            double **X = block_matrix(norbs,nmos);
            for (int m = 0 ; m<norbs; m++)
                for (int i = 0; i<nmos; i++)
                    X[m][i] = X_->get(h,m,i);

            double **F = block_matrix(norbs,norbs);
            for (int m = 0 ; m<norbs; m++)
                for (int i = 0; i<norbs; i++)
                    F[m][i] = F_->get(h,m,i);
            
            double **C = block_matrix(norbs,nmos);
            double **Temp = block_matrix(nmos,norbs);
            double **Fp = block_matrix(nmos,nmos);
            double **Cp = block_matrix(nmos,nmos);
    
            //Form F' = X'FX for canonical orthogonalization 
            C_DGEMM('T','N',nmos,norbs,norbs,1.0,X[0],nmos,F[0],norbs,0.0,Temp[0],norbs);  
            C_DGEMM('N','N',nmos,nmos,norbs,1.0,Temp[0],norbs,X[0],nmos,0.0,Fp[0],nmos);  

            //Form C' = eig(F')
            double *eigvals = init_array(nmos);
            sq_rsp(nmos, nmos, Fp,  eigvals, 1, Cp, 1.0e-14);
            for (int i = 0; i<nmos; i++)
                eigval.set(h,i,eigvals[i]);    
            
            free(eigvals);    

            //Form C = XC'
            C_DGEMM('N','N',norbs,nmos,nmos,1.0,X[0],nmos,Cp[0],nmos,0.0,C[0],nmos);

            for (int m = 0 ; m<norbs; m++)
                for (int i = 0; i<nmos; i++)
                    C_->set(h,m,i,C[m][i]);

            free_block(X);
            free_block(F);
            free_block(C);
            free_block(Temp);
            free_block(Cp);
            free_block(Fp);
        }
        
        find_occupation(eigval);
    }
}

void RHF::form_D()
{
    int h, i, j;
    int *opi = D_->rowspi();
    int nirreps = D_->nirreps();
    int norbs = basisset_->nbf();

    double** C = C_->to_block_matrix();
    double** D = block_matrix(norbs,norbs);
    
    int offset = 0;
    for (h = 0; h<nirreps; ++h) {
        C_DGEMM('n','t',opi[h],opi[h],doccpi_[h],1.0,&C[offset][offset],norbs,&C[offset][offset],norbs,0.0,&D[offset][offset],norbs);
        
        for (i = 0; i<opi[h]; i++)
            for (int j = 0; j<opi[h]; j++)
                D_->set(h,i,j,D[offset+i][offset+j]);

        offset += opi[h]; 
    } 

    free_block(C);
    free_block(D);


#ifdef _DEBUG
    if (debug_) {
        D_->print(outfile);
    }
#endif
}

double RHF::compute_initial_E()
{
    double Etotal = nuclearrep_ + D_->vector_dot(H_);
    return Etotal;
}

double RHF::compute_E()
{
    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(F_);
    double Etotal = nuclearrep_ + D_->vector_dot(HplusF);
    return Etotal;
}

void RHF::form_PK()
{
    // struct iwlbuf ERIIN;
    int ilsti, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    size_t bra, ket, braket=0;
    int idx;
    int counter = 0;
    int pk_counter = 0;
    bool pk_flag = false;
    double value;

    // PK zeroed out during allocation
    fprintf(outfile, "  Forming PK matrix.\n");
    fflush(outfile);

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            if (ERIIN.labels()[fi] >= 0) {
                i = ERIIN.labels()[fi];
                pk_flag = false;
            }
            else {
                i = -ERIIN.labels()[fi];
                pk_flag = true;
            }
            i = ERIIN.labels()[fi] >= 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi];
            j = ERIIN.labels()[fi+1];
            k = ERIIN.labels()[fi+2];
            l = ERIIN.labels()[fi+3];
            value = ERIIN.values()[idx];
            fi += 4;

            // Get the symmetries
            is = so2symblk_[i];
            js = so2symblk_[j];
            ks = so2symblk_[k];
            ls = so2symblk_[l];

            // Get the offset of the SO index in its symblock
            ii = so2index_[i];
            jj = so2index_[j];
            kk = so2index_[k];
            ll = so2index_[l];

            // J
            if ((is == js) && (ks == ls)) {
                bra = INDEX2(ii, jj) + pk_symoffset_[is];
                ket = INDEX2(kk, ll) + pk_symoffset_[ks];
                // _pk_symoffset corrects for the symmetry offset in the _pk vector
                braket = INDEX2(bra, ket);
                pk_[braket] += value;
                // K/2 (2nd sort)
                if ((ii != jj) && (kk != ll)) {
                    if ((is == ls) && (js == ks)) {
                        bra = INDEX2(ii, ll) + pk_symoffset_[is];
                        ket = INDEX2(jj, kk) + pk_symoffset_[js];
                        braket = INDEX2(bra, ket);
                        if ((ii == ll) || (jj == kk))
                            pk_[braket] -= 0.5 * value;
                        else
                            pk_[braket] -= 0.25 * value;
                    }
                }
            }

            // K/2 (1st sort)
            if ((is == ks) && (js == ls)) {
                bra = INDEX2(ii, kk) + pk_symoffset_[is];
                ket = INDEX2(jj, ll) + pk_symoffset_[js];
                braket = INDEX2(bra, ket);
                if ((ii == kk) || (jj == ll))
                    pk_[braket] -= 0.5 * value;
                else
                    pk_[braket] -= 0.25 * value;
            }
            pk_counter++;
            counter++;

            if (pk_flag) {
                pk_counter = 0;
            }
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);

    // After stage two is complete, the elements of P must be halved for the case IJ=KL.
    for (size_t ij=0; ij < pk_pairs_; ++ij)
        pk_[INDEX2(ij,ij)] *= 0.5;

    fprintf(outfile, "  Processed %d two-electron integrals.\n\n", counter);
#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "_pk:\n");
        print_array(pk_, pk_pairs_, outfile);
    }
#endif
}

void RHF::form_G_from_PK()
{
    int nirreps = factory_.nirreps();
    int *opi = factory_.rowspi();
    size_t ij;
    double *D_vector = new double[pk_pairs_];
    double *G_vector = new double[pk_pairs_];

    G_->zero();
    memset(D_vector, 0, sizeof(double) * pk_pairs_);
    memset(G_vector, 0, sizeof(double) * pk_pairs_);

    ij=0;
    for (int h=0; h<nirreps; ++h) {
        for (int p=0; p<opi[h]; ++p) {
            for (int q=0; q<=p; ++q) {
                if (p != q)
                    D_vector[ij] = 2.0 * D_->get(h, p, q);
                else
                    D_vector[ij] = D_->get(h, p, q);
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "PK: ij = %lu\n", (unsigned long)ij);
        fflush(outfile);
        fprintf(outfile, "PK: D matrix:\n");
        D_->print(outfile);
        fprintf(outfile, "PK: D vector (appears to be OK):\n");
        for (ij=0; ij<pk_pairs_; ++ij)
            fprintf(outfile, "PK: D vector [%lu] = %20.14f\n", (unsigned long)ij, D_vector[ij]);
    }
#endif

    double G_pq,D_pq;
    double* D_rs;
    double* G_rs;
    int pq,rs;
    double* PK_block = pk_;
    int ts_pairs = pk_pairs_;
    for(pq = 0; pq < ts_pairs; ++pq){
        G_pq = 0.0;
        D_pq = D_vector[pq];
        D_rs = &D_vector[0];
        G_rs = &G_vector[0];
        for(rs = 0; rs <= pq; ++rs){
            G_pq += *PK_block * (*D_rs);
            *G_rs += *PK_block * D_pq;

            // fprintf(outfile, "pq=%d, rs=%d, G_pq=%f, PK_block=%f, D_rs=%f, G_rs=%f, D_pq=%f\n", pq, rs, G_pq, *PK_block, *D_rs, *G_rs, D_pq);

            ++D_rs;
            ++G_rs;
            ++PK_block;
        }
        G_vector[pq] += G_pq;
    }

    // Convert G to a matrix
    ij = 0;
    for(int h = 0; h < nirreps; ++h){
        for(int p = 0; p < opi[h]; ++p){
            for(int q = 0; q <= p; ++q){
                G_->set(h,p,q, 2.0 * G_vector[ij]);
                G_->set(h,q,p, 2.0 * G_vector[ij]);
                ij++;
            }
        }
    }

#ifdef _DEBUG
    if (debug_) {
        G_->print(outfile);
    }
#endif

    delete[](D_vector);
    delete[](G_vector);
}

void RHF::form_G_from_direct_integrals()
{
    timer_on("form_G_from_direct_integrals");
    double temp1, temp2, temp3, temp4, temp5, temp6;
    int itype;

    // Zero out the G matrix
    G_->zero();

    // Need to back-transform the density from SO to AO basis
    SimpleMatrix *D = D_->to_simple_matrix();

    // D->set_name("D (AO basis) pre-transform");
    // D->print();
    // D->back_transform(basisset_->uso_to_bf());
    // D->set_name("D (AO basis) post-transform");
    // D->print();

    // Need a temporary G in the AO basis
    SimpleMatrix G;
    factory_.create_simple_matrix(G, "G (AO basis)");
    G.zero();

    // Initialize an integral object
    // Begin factor ou

    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    if (eri_.get() == NULL) {
        eri_ = shared_ptr<TwoBodyInt>(integral.eri());
    }
    ShellCombinationsIterator iter = integral.shells_iterator();
    const double *buffer = eri_->buffer();
    // End factor out

    //fprintf(outfile, "      Computing integrals..."); fflush(outfile);
    int P, Q, R, S;
    int i, j, k, l;
    int index;
    double value;

    int count=0;
    for (iter.first(); !iter.is_done(); iter.next()) {
        P = iter.p();
        Q = iter.q();
        R = iter.r();
        S = iter.s();

        // Compute quartet
        timer_on("compute_shell");
        eri_->compute_shell(P, Q, R, S);
        timer_off("compute_shell");

        // From the quartet get all the integrals
        IntegralsIterator int_iter = iter.integrals_iterator();
        for (int_iter.first(); !int_iter.is_done(); int_iter.next()) {
            i = int_iter.i();
            j = int_iter.j();
            k = int_iter.k();
            l = int_iter.l();
            index = int_iter.index();
            value = buffer[index];

            // We only care about those greater that 1.0e-14
            if (fabs(value) > 1.0e-14) {
                itype = integral_type(i, j, k, l);
                switch(itype) {
                    case 1:
                    temp1 = D->get(i, i) * value;

                    G.add(i, i, temp1);
                    break;

                    case 2:
                    temp1 = D->get(k, k) * 2.0 * value;
                    temp2 = D->get(i, k) * value;
                    temp3 = D->get(i, i) * 2.0 * value;

                    G.add(i, i, temp1);
                    G.add(k, k, temp3);
                    G.add(i, k, -temp2);
                    G.add(k, i, -temp2);
                    break;

                    case 3:
                    temp1 = D->get(i, i) * value;
                    temp2 = D->get(i, l) * value * 2.0;

                    G.add(i, l, temp1);
                    G.add(l, i, temp1);
                    G.add(i, i, temp2);
                    break;

                    case 4:
                    temp1 = D->get(j, j) * value;
                    temp2 = D->get(i, j) * value * 2.0;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(j, j, temp2);
                    break;

                    case 5:
                    temp1 = D->get(i, j) * value * 3.0;
                    G.add(i, j, temp1);
                    G.add(j, i, temp1);

                    temp2 = D->get(i, i) * value;
                    temp3 = D->get(j, j) * value;
                    G.add(j, j, -temp2);
                    G.add(i, i, -temp3);
                    break;

                    case 6:
                    temp1 = D->get(k, l) * value * 4.0;
                    temp2 = D->get(i, l) * value;
                    temp3 = D->get(i, i) * value * 2.0;
                    temp4 = D->get(i, k) * value;

                    G.add(i, i, temp1);
                    G.add(i, k, -temp2);
                    G.add(k, i, -temp2);
                    G.add(k, l, temp3);
                    G.add(l, k, temp3);
                    G.add(i, l, -temp4);
                    G.add(l, i, -temp4);
                    break;

                    case 7:
                    temp1 = D->get(i, j) * value * 4.0;
                    temp2 = D->get(j, k) * value;
                    temp3 = D->get(i, k) * value;
                    temp4 = D->get(k, k) * value * 2.0;

                    G.add(k, k,  temp1);
                    G.add(i, k, -temp2);
                    G.add(k, i, -temp2);
                    G.add(j, k, -temp3);
                    G.add(k, j, -temp3);
                    G.add(i, j,  temp4);
                    G.add(j, i,  temp4);
                    break;

                    case 8:
                    temp1 = D->get(k, k) * value * 2.0;
                    temp2 = D->get(i, j) * value * 4.0;
                    temp3 = D->get(j, k) * value;
                    temp4 = D->get(i, k) * value;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(k, k, temp2);
                    G.add(i, k, -temp3);
                    G.add(k, i, -temp3);
                    G.add(j, k, -temp4);
                    G.add(k, j, -temp4);
                    break;

                    case 9:
                    temp1 = D->get(i, l) * value * 3.0;
                    temp2 = D->get(i, j) * value * 3.0;
                    temp3 = D->get(j, l) * value * 2.0;
                    temp4 = D->get(i, i) * value;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(i, l, temp2);
                    G.add(l, i, temp2);
                    G.add(i, i, -temp3);
                    G.add(j, l, -temp4);
                    G.add(l, j, -temp4);
                    break;

                    case 10:
                    temp1 = D->get(j, l) * value * 3.0;
                    temp2 = D->get(i, j) * value * 3.0;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, l) * value * 2.0;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(j, l, temp2);
                    G.add(l, j, temp2);
                    G.add(i, l, -temp3);
                    G.add(l, i, -temp3);
                    G.add(j, j, -temp4);
                    break;

                    case 11:
                    temp1 = D->get(k, j) * value * 3.0;
                    temp2 = D->get(i, j) * value * 3.0;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, k) * value * 2.0;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(k, j, temp2);
                    G.add(j, k, temp2);
                    G.add(i, k, -temp3);
                    G.add(k, i, -temp3);
                    G.add(j, j, -temp4);
                    break;

                    case 12:
                    case 13:
                    case 14:
                    temp1 = D->get(k, l) * value * 4.0;
                    temp2 = D->get(i, j) * value * 4.0;
                    temp3 = D->get(j, l) * value;
                    temp4 = D->get(i, k) * value;
                    temp5 = D->get(j, k) * value;
                    temp6 = D->get(i, l) * value;

                    G.add(i, j, temp1);
                    G.add(j, i, temp1);
                    G.add(k, l, temp2);
                    G.add(l, k, temp2);
                    G.add(i, k, -temp3);
                    G.add(k, i, -temp3);
                    G.add(j, l, -temp4);
                    G.add(l, j, -temp4);
                    G.add(i, l, -temp5);
                    G.add(l, i, -temp5);
                    G.add(j, k, -temp6);
                    G.add(k, j, -temp6);
                    break;
                };

                count++;
            }
        }
    }
    //fprintf(outfile, "done.  %d two-electron integrals.\n", count); fflush(outfile);
//    delete eri;

    // Set RefMatrix to RefSimpleMatrix handling symmetry blocking, if needed
    // Transform G back to symmetry blocking
    // G.transform(basisset_->uso_to_bf());
    // G.print();
    G_->set(&G);
    delete D;
    // G_->print();
    timer_off("form_G_from_direct_integrals");
}

void RHF::form_G()
{
    // struct iwlbuf ERIIN;
    int ilsti, nbuf;
    int i, j, k, l;
    int ii, jj, kk, ll;
    int is, js, ks, ls;
    int fi;
    double value=0.0;
    double temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0, temp5=0.0, temp6=0.0;
    int idx;
    int itype;
    int counter = 0;

    // Zero out the G matrix
    G_->zero();

    IWL ERIIN(psio_.get(), PSIF_SO_TEI, 0.0, 1, 1);

    do {
        ilsti = ERIIN.last_buffer();
        nbuf  = ERIIN.buffer_count();

        fi = 0;
        for (idx=0; idx<nbuf; ++idx) {
            i = ERIIN.labels()[fi] > 0 ? ERIIN.labels()[fi] : -ERIIN.labels()[fi];
            j = ERIIN.labels()[fi+1];
            k = ERIIN.labels()[fi+2];
            l = ERIIN.labels()[fi+3];
            value = ERIIN.values()[idx];
            fi += 4;

            itype = integral_type(i, j, k, l);
            switch(itype) {
                case 1:
                    ii = so2index_[i];
                    is = so2symblk_[i];
                    temp1 = D_->get(is, ii, ii) * value;

                    G_->add(is, ii, ii, temp1);

// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 1:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is, ii,ii));
//                     }
// #endif
                    break;

                    case 2:
                    ii = so2index_[i];
                    kk = so2index_[k];
                    is = so2symblk_[i];
                    ks = so2symblk_[k];

                    temp1 = D_->get(ks, kk, kk) * 2.0 * value;
                    temp2 = 0.0;
                    temp3 = D_->get(is, ii, ii) * 2.0 * value;

                    G_->add(is, ii, ii, temp1);
                    G_->add(ks, kk, kk, temp3);

                    if (is == ks) {
                        temp2 = D_->get(is, ii, kk) * value;
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 2:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, kk, G_->get(ks,kk,kk));
//                     }
// #endif
                    break;

                    case 3:
                    ii = so2index_[i];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    ls = so2symblk_[l];

                    temp1 = temp2 = 0.0;
                    if (is == ls) {
                        temp1 = D_->get(is, ii, ii) * value;
                        temp2 = D_->get(is, ii, ll) * value * 2.0;

                        G_->add(is, ii, ll, temp1);
                        G_->add(is, ll, ii, temp1);
                        G_->add(is, ii, ii, temp2);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 3:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                     }
// #endif
                    break;


                    case 4:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    temp1 = temp2 = 0.0;
                    if (is == js) {
                        temp1 = D_->get(js, jj, jj) * value;
                        temp2 = D_->get(is, ii, jj) * value * 2.0;

                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                        G_->add(js, jj, jj, temp2);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 4:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                     }
// #endif
                    break;

                    case 5:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    if (is == js) {
                        temp1 = D_->get(is, ii, jj) * value * 3.0;
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }

                    temp2 = D_->get(is, ii, ii) * value;
                    temp3 = D_->get(js, jj, jj) * value;
                    G_->add(js, jj, jj, -temp2);
                    G_->add(is, ii, ii, -temp3);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 5:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                     }
// #endif
                    break;

                    case 6:
                    ii = so2index_[i];
                    kk = so2index_[k];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    ks = so2symblk_[k];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (ks == ls)
                        temp1 = D_->get(ks, kk, ll) * value * 4.0;
                    if (is == ls)
                        temp2 = D_->get(is, ii, ll) * value;
                    temp3 = D_->get(is, ii, ii) * value * 2.0;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;

                    G_->add(is, ii, ii, temp1);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    if (ks == ls) {
                        G_->add(ks, kk, ll, temp3);
                        G_->add(ks, ll, kk, temp3);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp4);
                        G_->add(is, ll, ii, -temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 6:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, ll, G_->get(ks,kk,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                     }
// #endif
                    break;

                    case 7:
                    kk = so2index_[k];
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ks = so2symblk_[k];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (is == js)
                        temp1 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ks)
                        temp2 = D_->get(js, jj, kk) * value;
                    if (is == ks)
                        temp3 = D_->get(is, ii, kk) * value;
                    temp4 = D_->get(ks, kk, kk) * value * 2.0;

                    G_->add(ks, kk, kk, temp1);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp2);
                        G_->add(is, kk, ii, -temp2);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp3);
                        G_->add(js, kk, jj, -temp3);
                    }
                    if (is == js) {
                        G_->add(is, ii, jj, temp4);
                        G_->add(is, jj, ii, temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 7:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, kk, G_->get(ks,kk,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, kk, G_->get(js,jj,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                     }
// #endif
                    break;

                    case 8:
                    kk = so2index_[k];
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ks = so2symblk_[k];
                    is = so2symblk_[i];
                    js = so2symblk_[j];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    temp1 = D_->get(ks, kk, kk) * value * 2.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ks)
                        temp3 = D_->get(js, jj, kk) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    G_->add(ks, kk, kk, temp2);
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp4);
                        G_->add(js, kk, jj, -temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 8:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, kk, G_->get(ks,kk,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, kk, G_->get(js,jj,kk));
//                     }
// #endif
                    break;

                    case 9:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    js = so2symblk_[j];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (is == ls)
                        temp1 = D_->get(is, ii, ll) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    if (js == ls)
                        temp3 = D_->get(js, jj, ll) * value * 2.0;
                    temp4 = D_->get(is, ii, ii) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, temp2);
                        G_->add(is, ll, ii, temp2);
                    }
                    G_->add(is, ii, ii, -temp3);
                    if (js == ls) {
                        G_->add(js, jj, ll, -temp4);
                        G_->add(js, ll, jj, -temp4);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 9:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ii, G_->get(is,ii,ii));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, ll, G_->get(js,jj,ll));
//                     }
// #endif
                    break;

                    case 10:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    js = so2symblk_[j];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (js == ls)
                        temp1 = D_->get(js, jj, ll) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    temp3 = D_->get(js, jj, jj) * value;
                    if (is == ls)
                        temp4 = D_->get(is, ii, ll) * value * 2.0;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (js == ls) {
                        G_->add(js, jj, ll, temp2);
                        G_->add(js, ll, jj, temp2);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp3);
                        G_->add(is, ll, ii, -temp3);
                    }
                    G_->add(js, jj, jj, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 10:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, ll, G_->get(js,jj,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                     }
// #endif
                    break;

                    case 11:
                    ii = so2index_[i];
                    kk = so2index_[k];
                    jj = so2index_[j];
                    is = so2symblk_[i];
                    ks = so2symblk_[k];
                    js = so2symblk_[j];

                    temp1 = temp2 = temp3 = temp4 = 0.0;
                    if (ks == js)
                        temp1 = D_->get(ks, kk, jj) * value * 3.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 3.0;
                    temp3 = D_->get(js, jj, jj) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value * 2.0;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (ks == js) {
                        G_->add(ks, kk, jj, temp2);
                        G_->add(ks, jj, kk, temp2);
                    }
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    G_->add(js, jj, jj, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 11:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, jj, G_->get(ks,kk,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, jj, G_->get(js,jj,jj));
//                     }
// #endif
                    break;

                    case 12:
                    case 13:
                    case 14:
                    ii = so2index_[i];
                    jj = so2index_[j];
                    kk = so2index_[k];
                    ll = so2index_[l];
                    is = so2symblk_[i];
                    js = so2symblk_[j];
                    ks = so2symblk_[k];
                    ls = so2symblk_[l];

                    temp1 = temp2 = temp3 = temp4 = temp5 = temp6 = 0.0;
                    if (ks == ls)
                        temp1 = D_->get(ks, kk, ll) * value * 4.0;
                    if (is == js)
                        temp2 = D_->get(is, ii, jj) * value * 4.0;
                    if (js == ls)
                        temp3 = D_->get(js, jj, ll) * value;
                    if (is == ks)
                        temp4 = D_->get(is, ii, kk) * value;
                    if (js == ks)
                        temp5 = D_->get(js, jj, kk) * value;
                    if (is == ls)
                        temp6 = D_->get(is, ii, ll) * value;

                    if (is == js) {
                        G_->add(is, ii, jj, temp1);
                        G_->add(is, jj, ii, temp1);
                    }
                    if (ks == ls) {
                        G_->add(ks, kk, ll, temp2);
                        G_->add(ks, ll, kk, temp2);
                    }
                    if (is == ks) {
                        G_->add(is, ii, kk, -temp3);
                        G_->add(is, kk, ii, -temp3);
                    }
                    if (js == ls) {
                        G_->add(js, jj, ll, -temp4);
                        G_->add(js, ll, jj, -temp4);
                    }
                    if (is == ls) {
                        G_->add(is, ii, ll, -temp5);
                        G_->add(is, ll, ii, -temp5);
                    }
                    if (js == ks) {
                        G_->add(js, jj, kk, -temp6);
                        G_->add(js, kk, jj, -temp6);
                    }
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 12,13,14:\n");
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, jj, G_->get(is,ii,jj));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", ks, kk, ll, G_->get(ks,kk,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, kk, G_->get(is,ii,kk));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, ll, G_->get(js,jj,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", is, ii, ll, G_->get(is,ii,ll));
//                         fprintf(outfile, "G[%d][%d][%d] %.11f\n", js, jj, kk, G_->get(js,jj,kk));
//                     }
// #endif
                    break;
            };
            counter++;
        }

        if (!ilsti)
            ERIIN.fetch();
    } while (!ilsti);

    // Going out of scope will close the buffer
    // iwl_buf_close(&ERIIN, 1);
    fprintf(outfile, "  Processed %6d two-electron integrals.\n", counter);
}

void RHF::form_G_from_J_and_K(double scale_K_by)
{
    form_J_and_K();
    //J_->print();
    //K_->print();
    J_->scale(2.0);
    K_->scale(-1.0);
    G_->copy(K_);
    G_->scale(scale_K_by);
    G_->add(J_);
    //G_->print();
}

void RHF::form_J_and_K()
{
    if (scf_type_ == "PK") {
        //form_J_from_PK();
        //form_K_from_PK();
    }
    else if (scf_type_ == "DIRECT") {
        form_J_and_K_from_direct_integrals();
    }
    else if (scf_type_ == "DF" || scf_type_ == "CD" || scf_type_ == "1C_CD") {
        form_J_from_RI();
        form_K_from_RI();
    }
    else {
        //form_J();
        //form_K();
    }

}

void RHF::form_J_and_K_from_direct_integrals()
{
    double temp1, temp2, temp3, temp4, temp5, temp6;
    int itype;

    // Zero out the J and K matrices
    J_->zero();
    K_->zero();

    //D_->zero();
    //D_->set(0,0,0,2.0000008);
    //D_->set(0,5,5,2.0000008);


    // Need to back-transform the density from SO to AO basis
    SimpleMatrix *D = D_->to_simple_matrix();

    // D->set_name("D (AO basis) pre-transform");
    // D->print();
    // D->back_transform(basisset_->uso_to_bf());
    // D->set_name("D (AO basis) post-transform");
    // D->print();

    // Need a temporary J in the AO basis
    SimpleMatrix J;
    factory_.create_simple_matrix(J, "J (AO basis)");
    J.zero();
    // Need a temporary K in the AO basis
    SimpleMatrix K;
    factory_.create_simple_matrix(K, "K (AO basis)");
    K.zero();

    // Initialize an integral object
    // Begin factor out
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    if (eri_.get() == NULL) {
        eri_ = shared_ptr<TwoBodyInt>(integral.eri());
    }
    ShellCombinationsIterator iter = integral.shells_iterator();
    const double *buffer = eri_->buffer();
    // End factor out

    //fprintf(outfile, "\n      Computing integrals..."); fflush(outfile);
    int P, Q, R, S;
    int i, j, k, l;
    int index;
    double value;

    int count=0;
    for (iter.first(); !iter.is_done(); iter.next()) {
        P = iter.p();
        Q = iter.q();
        R = iter.r();
        S = iter.s();

        // Compute quartet
        eri_->compute_shell(P, Q, R, S);

        // fprintf(outfile, "Doing shell ( %d %d | %d %d )\n", P, Q, R, S); fflush(outfile);

        // From the quartet get all the integrals
        IntegralsIterator int_iter = iter.integrals_iterator();
        for (int_iter.first(); !int_iter.is_done(); int_iter.next()) {
            i = int_iter.i();
            j = int_iter.j();
            k = int_iter.k();
            l = int_iter.l();
            index = int_iter.index();
            value = buffer[index];

            //fprintf(outfile, "\tDoing integral ( %d %d | %d %d )\n", i, j, k, l); fflush(outfile);
            //fprintf(outfile, "\n (%d, %d| %d, %d) = %20.10f", i, j, k, l, value); fflush(outfile);
            // We only care about those greater that 1.0e-14
            if (fabs(value) > 1.0e-14) {
// #ifdef _DEBUG
//                 if (debug_)

// #endif
                itype = integral_type(i, j, k, l);
                switch(itype) {
                    case 1:
                    temp1 = D->get(i, i) * value;

                    J.add(i, i, 1.0*temp1);
                    K.add(i, i, -1.0*temp1);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 1:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, i, G.get(i,i), D->get(i, i));
//                     }
// #endif
                    break;

                    case 2:
                    temp1 = D->get(k, k) * 1.0 * value;
                    temp2 = D->get(i, k) * value;
                    temp3 = D->get(i, i) * 1.0 * value;

                    J.add(i, i, temp1);
                    J.add(k, k, temp3);
                    K.add(i, k, -temp2);
                    K.add(k, i, -temp2);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 2:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, i, G.get(i,i), D->get(k, k));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, k, G.get(i,k), D->get(i, k));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", k, k, G.get(k,k), D->get(i, i));
//                     }
// #endif
                    break;

                    case 3:
                    temp1 = D->get(i, i) * value;
                    temp2 = D->get(i, l) * value * 2.0;

                    J.add(i, l, 1.0*temp1);
                    J.add(l, i, 1.0*temp1);
                    K.add(i, l, -temp1);
                    K.add(l, i, -temp1);
                    J.add(i, i, 1.0*temp2);
                    K.add(i, i, -temp2);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 3:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, l, G.get(i,l), D->get(i, i));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, i, G.get(i,i), D->get(i, l));
//                     }
// #endif
                    break;

                    case 4:
                    temp1 = D->get(j, j) * value;
                    temp2 = D->get(i, j) * value * 2.0;

                    J.add(i, j, 1.0*temp1);
                    J.add(j, i, 1.0*temp1);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    J.add(j, j, 1.0*temp2);
                    K.add(j, j, -temp2);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 4:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", i, j, G.get(i,j), D->get(j, j));
//                         fprintf(outfile, "G[%d][%d] %.11f %.11f\n", j, j, G.get(j,j), D->get(i, j));
//                     }
// #endif
                    break;

                    case 5:
                    temp1 = D->get(i, j) * value;
                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);

                    temp2 = D->get(i, i) * value;
                    temp3 = D->get(j, j) * value;
                    K.add(j, j, -temp2);
                    K.add(i, i, -temp3);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 5:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, j, G.get(j,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, i, G.get(i,i));
//                     }
// #endif
                    break;

                    case 6:
                    temp1 = D->get(k, l) * value * 2.0;
                    temp2 = D->get(i, l) * value;
                    temp3 = D->get(i, i) * value * 1.0;
                    temp4 = D->get(i, k) * value;

                    J.add(i, i, temp1);
                    K.add(i, k, -temp2);
                    K.add(k, i, -temp2);
                    J.add(k, l, temp3);
                    J.add(l, k, temp3);
                    K.add(i, l, -temp4);
                    K.add(l, i, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 6:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, i, G.get(i,i));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, l, G.get(k,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                     }
// #endif
                    break;

                    case 7:
                    temp1 = D->get(i, j) * value * 2.0;
                    temp2 = D->get(j, k) * value;
                    temp3 = D->get(i, k) * value;
                    temp4 = D->get(k, k) * value * 1.0;

                    J.add(k, k,  temp1);
                    K.add(i, k, -temp2);
                    K.add(k, i, -temp2);
                    K.add(j, k, -temp3);
                    K.add(k, j, -temp3);
                    J.add(i, j,  temp4);
                    J.add(j, i,  temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 7:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, k, G.get(k,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, k, G.get(j,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                     }
// #endif
                    break;

                    case 8:
                    temp1 = D->get(k, k) * value * 1.0;
                    temp2 = D->get(i, j) * value * 2.0;
                    temp3 = D->get(j, k) * value;
                    temp4 = D->get(i, k) * value;

                    J.add(i, j, temp1);
                    J.add(j, i, temp1);
                    J.add(k, k, temp2);
                    K.add(i, k, -temp3);
                    K.add(k, i, -temp3);
                    K.add(j, k, -temp4);
                    K.add(k, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 8:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, k, G.get(k,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, k, G.get(j,k));
//                     }
// #endif
                    break;

                    case 9:
                    temp1 = D->get(i, l) * value;
                    temp2 = D->get(i, j) * value;
                    temp3 = D->get(j, l) * value * 2.0;
                    temp4 = D->get(i, i) * value;

                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    J.add(i, l, 2.0*temp2);
                    J.add(l, i, 2.0*temp2);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    K.add(i, l, -temp2);
                    K.add(l, i, -temp2);
                    K.add(i, i, -temp3);
                    K.add(j, l, -temp4);
                    K.add(l, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 9:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, i, G.get(i,i));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, l, G.get(j,l));
//                     }
// #endif
                    break;

                    case 10:
                    temp1 = D->get(j, l) * value;
                    temp2 = D->get(i, j) * value;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, l) * value * 2.0;

                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    J.add(j, l, 2.0*temp2);
                    J.add(l, j, 2.0*temp2);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    K.add(j, l, -temp2);
                    K.add(l, j, -temp2);
                    K.add(i, l, -temp3);
                    K.add(l, i, -temp3);
                    K.add(j, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 10:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, l, G.get(j,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, j, G.get(j,j));
//                     }
// #endif
                    break;

                    case 11:
                    temp1 = D->get(k, j) * value;
                    temp2 = D->get(i, j) * value;
                    temp3 = D->get(j, j) * value;
                    temp4 = D->get(i, k) * value * 2.0;

                    J.add(i, j, 2.0*temp1);
                    J.add(j, i, 2.0*temp1);
                    J.add(k, j, 2.0*temp2);
                    J.add(j, k, 2.0*temp2);
                    K.add(i, j, -temp1);
                    K.add(j, i, -temp1);
                    K.add(k, j, -temp2);
                    K.add(j, k, -temp2);
                    K.add(i, k, -temp3);
                    K.add(k, i, -temp3);
                    K.add(j, j, -temp4);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 11:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, j, G.get(k,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, j, G.get(j,j));
//                     }
// #endif
                    break;

                    case 12:
                    case 13:
                    case 14:
                    temp1 = D->get(k, l) * value * 2.0;
                    temp2 = D->get(i, j) * value * 2.0;
                    temp3 = D->get(j, l) * value;
                    temp4 = D->get(i, k) * value;
                    temp5 = D->get(j, k) * value;
                    temp6 = D->get(i, l) * value;

                    J.add(i, j, temp1);
                    J.add(j, i, temp1);
                    J.add(k, l, temp2);
                    J.add(l, k, temp2);
                    K.add(i, k, -temp3);
                    K.add(k, i, -temp3);
                    K.add(j, l, -temp4);
                    K.add(l, j, -temp4);
                    K.add(i, l, -temp5);
                    K.add(l, i, -temp5);
                    K.add(j, k, -temp6);
                    K.add(k, j, -temp6);
// #ifdef _DEBUG
//                     if (debug_) {
//                         fprintf(outfile, "INTEGRAL CASE 12,13,14:\n");
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, j, G.get(i,j));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", k, l, G.get(k,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, k, G.get(i,k));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, l, G.get(j,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", i, l, G.get(i,l));
//                         fprintf(outfile, "G[%d][%d] %.11f\n", j, k, G.get(j,k));
//                     }
// #endif
                    break;
                };

                count++;
            }
        }
    }
    //fprintf(outfile, "done. %d two-electron integrals.\n", count); fflush(outfile);
//    delete eri;

    // Set RefMatrix to RefSimpleMatrix handling symmetry blocking, if needed
    // Transform G back to symmetry blocking
    // G.transform(basisset_->uso_to_bf());
    // G.print();
    J_->set(&J);
    K_->set(&K);
    J_->scale(1.0);
    K_->scale(-1.0);
    delete D;
    // G_->print();
}
void RHF::save_sapt_info()
{
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1. Period.\n"); fflush(outfile);
        abort();
    }
    if (soccpi_[0] != 0)
    {
        fprintf(outfile,"Aren't we in RHF Here? Pair those electrons up cracker!\n"); fflush(outfile);
        abort();
    }

    int fileno;
    char* body_type = (char*)malloc(400*sizeof(char));
    char* key_buffer = (char*)malloc(4000*sizeof(char));

    string type = options_.get_str("SAPT");
    if (type == "2-DIMER") {
        fileno = PSIF_SAPT_DIMER;
        sprintf(body_type,"Dimer");
    } else if (type == "2-MONOMER_A") {
        fileno = PSIF_SAPT_MONOMERA;
        sprintf(body_type,"Monomer");
    } else if (type == "2-MONOMER_B") {
        fileno = PSIF_SAPT_MONOMERB;
        sprintf(body_type,"Monomer");
    } else if (type == "3-TRIMER") {
        fileno = PSIF_3B_SAPT_TRIMER;
        sprintf(body_type,"Trimer");
    } else if (type == "3-DIMER_AB") {
        fileno = PSIF_3B_SAPT_DIMER_AB;
        sprintf(body_type,"Dimer");
    } else if (type == "3-DIMER_BC") {
        fileno = PSIF_3B_SAPT_DIMER_BC;
        sprintf(body_type,"Dimer");
    } else if (type == "3-DIMER_AC") {
        fileno = PSIF_3B_SAPT_DIMER_AC;
        sprintf(body_type,"Dimer");
    } else if (type == "3-MONOMER_A") {
        fileno = PSIF_3B_SAPT_MONOMER_A;
        sprintf(body_type,"Monomer");
    } else if (type == "3-MONOMER_B") {
        fileno = PSIF_3B_SAPT_MONOMER_B;
        sprintf(body_type,"Monomer");
    } else if (type == "3-MONOMER_C") {
        fileno = PSIF_3B_SAPT_MONOMER_C;
        sprintf(body_type,"Monomer");
    } else {
        throw std::domain_error("SAPT Output option invalid");
    }

    psio_open(fileno,0);

    int sapt_nso = basisset_->nbf();
    int sapt_nmo = basisset_->nbf();
    int sapt_nocc = doccpi_[0];
    int sapt_nvir = sapt_nso-sapt_nocc;
    int sapt_ne = 2*sapt_nocc;
    double sapt_E_HF = E_;
    double sapt_E_nuc = nuclearrep_;
    double *sapt_evals = chkpt_->rd_evals();
    double **sapt_C = C_->to_block_matrix();
    //print_mat(sapt_C,sapt_nso,sapt_nso,outfile);
    SharedMatrix potential(factory_.create_matrix("Potential Integrals"));
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    shared_ptr<OneBodyInt> V(integral.potential());
    V->compute(potential);
    double *sapt_V_ints = potential->to_lower_triangle();
    double *sapt_S_ints = S_->to_lower_triangle();

    int errcod;

    sprintf(key_buffer,"%s NSO",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_nso, sizeof(int));
    sprintf(key_buffer,"%s NMO",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_nmo, sizeof(int));
    sprintf(key_buffer,"%s NOCC",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_nocc, sizeof(int));
    sprintf(key_buffer,"%s NVIR",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_nvir, sizeof(int));
    sprintf(key_buffer,"%s Number of Electrons",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_ne,sizeof(int));
    sprintf(key_buffer,"%s HF Energy",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_E_HF,sizeof(double));
    sprintf(key_buffer,"%s Nuclear Repulsion Energy",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &sapt_E_nuc, sizeof(double));
    sprintf(key_buffer,"%s HF Eigenvalues",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &(sapt_evals[0]),sizeof(double)*sapt_nmo);
    sprintf(key_buffer,"%s Nuclear Attraction Integrals",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &(sapt_V_ints[0]), sizeof(double)*sapt_nso*(sapt_nso+1)/2);
    sprintf(key_buffer,"%s Overlap Integrals",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &(sapt_S_ints[0]), sizeof(double)*sapt_nso*(sapt_nso+1)/2);
    sprintf(key_buffer,"%s HF Coefficients",body_type);
    errcod = psio_write_entry(fileno,key_buffer,(char *) &(sapt_C[0][0]),sizeof(double)*sapt_nmo*sapt_nso);

    psio_close(fileno,1);

    free(sapt_evals);
    free(sapt_V_ints);
    free(sapt_S_ints);
    free_block(sapt_C);

    free(body_type);
    free(key_buffer);
}
void RHF::compute_SAD_guess() 
{

    shared_ptr<Molecule> mol = basisset_->molecule();
    std::vector<shared_ptr<BasisSet> > atomic_bases;

    if (print_ > 6) {
        fprintf(outfile,"\n  Constructing atomic basis sets\n  Molecule:\n");
        mol->print();
    }

    //Build the atomic basis sets for libmints use in UHF
    for (int A = 0; A<mol->nallatom(); A++) {
        atomic_bases.push_back(basisset_->atomic_basis_set(A));
        if (print_>6) {
            fprintf(outfile,"  Atomic Basis Set %d\n", A);
            atomic_bases[A]->molecule()->print();
            fprintf(outfile,"\n");
            atomic_bases[A]->print(outfile);
            fprintf(outfile,"\n");
        }
    }

    //Spin occupations per atom, to be determined by Hund's Rules
    //or user input
    int* nalpha = init_int_array(mol->nallatom());
    int* nbeta = init_int_array(mol->nallatom());
    int* nelec = init_int_array(mol->nallatom());
    int* nhigh = init_int_array(mol->nallatom());

    //Ground state high spin occupency array, atoms 0 to 36 (see Giffith's Quantum Mechanics, pp. 217)
    const int reference_S[] = {0,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,6,5,4,3,2,1,0,1,2,3,2,1,0};

    //At the moment, we'll assume no ions, and Hund filling
    //Improvemnts to SAD can be made simply by changing the setup of nalpha and nbeta
    if (print_>3)
        fprintf(outfile,"\n  Determining atomic occupations:\n");
    for (int A = 0; A<mol->nallatom(); A++) {
        int Z = mol->fZ(A);
        if (Z>36) {
            throw std::domain_error("Atoms up to 36 (Kr) are currently supported with SAD");
        }
        nhigh[A] = reference_S[Z];
        nelec[A] = Z;
        nbeta[A] = (nelec[A]-nhigh[A])/2;
        nalpha[A] = nelec[A]-nbeta[A];
        if (print_>6)
            fprintf(outfile,"  Atom %d, Z = %d, nelec = %d, nohigh = %d, nalpha = %d, nbeta = %d\n",A,Z,nelec[A],nhigh[A],nalpha[A],nbeta[A]);
    }

    timer_on("Atomic UHF");

    //Atomic D matrices within the atom specific AO basis
    double*** atomic_D = (double***)malloc(mol->nallatom()*sizeof(double**));
    for (int A = 0; A<mol->nallatom(); A++)
        atomic_D[A] = block_matrix(atomic_bases[A]->nbf(),atomic_bases[A]->nbf());

    if (print_>2)
        fprintf(outfile,"\n  Performing Atomic UHF Computations:\n");
    for (int A = 0; A<mol->nallatom(); A++) {
        if (print_>4)
            fprintf(outfile,"\n  UHF Computation for Atom %d:\n",A);
        getUHFAtomicDensity(atomic_bases[A],nelec[A],nhigh[A],atomic_D[A]);
    }
    timer_off("Atomic UHF");

    //Add atomic_D into D (scale by 1/2, we like pairs in RHF)
    D_->zero();
    for (int A = 0, offset = 0; A < mol->nallatom(); A++) {
        int norbs = atomic_bases[A]->nbf();
        for (int m = 0; m<norbs; m++)
            for (int n = 0; n<norbs; n++)
                D_->set(0,m+offset,n+offset,0.5*atomic_D[A][m][n]);
        offset+=norbs;
    }
    if (print_>6)
        D_->print(outfile);

    //D_->print(outfile);    

    //A C matrix is needed. Do one of:
    //   --An integral direct step (expensive)
    //   --A Cholesky orbital approximation (cheap, but not as accurate) 
    if (options_.get_str("SAD_C") == "ID") {
    //>>>>>>>>ID SAD GUESS 
    //Compute a rough Fock matrix via integral direct
    //Note, my convention is backwards, l and s are index variables
    // m and n are zip variables
    timer_on("SAD Fock");
    double SAD_Schwarz = options_.get_double("SAD_SCHWARZ_CUTOFF");
    G_->zero();
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    TwoBodyInt *TEI = integral.eri(0, SAD_Schwarz);
    const double* buffer = TEI->buffer();
    for (int A = 0, shell_offset = 0; A<mol->nallatom(); A++) {
    for (int MU = shell_offset; MU < shell_offset+atomic_bases[A]->nshell(); MU++) {
    int numMU = basisset_->shell(MU)->nfunction();
    for (int NU = shell_offset; NU < shell_offset+atomic_bases[A]->nshell() ; NU++) {
    int numNU = basisset_->shell(NU)->nfunction();
    for (int LA = 0; LA < basisset_->nshell(); LA++) {
    int numLA = basisset_->shell(LA)->nfunction();
    for (int SI = 0; SI < basisset_->nshell(); SI++) {
    int numSI = basisset_->shell(SI)->nfunction();

    //J
    if (!TEI->shell_is_zero(MU,NU,LA,SI)) {
    TEI->compute_shell(MU,NU,LA,SI);
    for (int m = 0, index = 0; m<numMU; m++) {
    int omu = basisset_->shell(MU)->function_index() + m;
    for (int n = 0; n<numNU; n++) {
    int onu = basisset_->shell(NU)->function_index() + n;
    for (int l = 0; l<numLA; l++) {
    int ola = basisset_->shell(LA)->function_index() + l;
    for (int s = 0; s<numSI; s++,index++) {
    int osi = basisset_->shell(SI)->function_index() + s;
        G_->add(0,ola,osi,2.0*D_->get(0,omu,onu)*buffer[index]);
    }
    }
    }
    }
    }

    //K
    if (!TEI->shell_is_zero(LA,MU,NU,SI)) {
    TEI->compute_shell(LA,MU,NU,SI);
    for (int l = 0, index = 0; l<numLA; l++) {
    int ola = basisset_->shell(LA)->function_index() + l;
    for (int m = 0; m<numMU; m++) {
    int omu = basisset_->shell(MU)->function_index() + m;
    for (int n = 0; n<numNU; n++) {
    int onu = basisset_->shell(NU)->function_index() + n;
    for (int s = 0; s<numSI; s++,index++) {
    int osi = basisset_->shell(SI)->function_index() + s;
        G_->add(0,ola,osi,-D_->get(0,omu,onu)*buffer[index]);
    }
    }
    }
    }
    }
    //fprintf(outfile,"  Cholesky Tensor: ntri_ = %d, ri_nbf_ = %d:\n",ntri_,ri_nbf_);

    }
    }
    }
    }
    shell_offset+= atomic_bases[A]->nshell();
    }

    F_->copy(H_);
    F_->add(G_);

    //Compute initial E for reference
    E_ = compute_E();

    //Form the C matrix from the rough Fock matrix
    form_C();
    //Form the D matrix from the resultant C matrix
    form_D();
    timer_off("SAD Fock");

    if (print_>7) {
        G_->print(outfile);
        F_->print(outfile);
        C_->print(outfile);
        D_->print(outfile);
    }//<<<< END ID SAD GUESS
    } else if (options_.get_str("SAD_C") == "CHOLESKY") {
        timer_on("SAD Cholesky");
        if (print_) {
            fprintf(outfile,"  Approximating occupied orbitals via Partial Cholesky Decomposition.\n");
            fprintf(outfile,"  NOTE: The first SCF iteration will not be variational.\n");
        }
        int norbs = nso_;
        int ndocc = nalpha_;

        double** D = D_->to_block_matrix();
        double* Temp = init_array(norbs);    
        int* P = init_int_array(norbs);
        for (int i = 0; i<norbs; i++)
            P[i] = i;
    
        //fprintf(outfile,"  D:\n");
        //print_mat(D,norbs,norbs,outfile);        
        //Pivot
        double max;
        int Temp_p;
        int pivot;
        for (int i = 0; i<norbs-1; i++) {
            max = 0.0;
            //Where's the pivot diagonal?
            for (int j = i; j<norbs; j++)
                if (max <= fabs(D[j][j])) {
                    max = fabs(D[j][j]);
                    pivot = j;
                }
        
            //Rows
            C_DCOPY(norbs,&D[pivot][0],1,Temp,1);            
            C_DCOPY(norbs,&D[i][0],1,&D[pivot][0],1);            
            C_DCOPY(norbs,Temp,1,&D[i][0],1);            

            //Columns
            C_DCOPY(norbs,&D[0][pivot],norbs,Temp,1);            
            C_DCOPY(norbs,&D[0][i],norbs,&D[0][pivot],norbs);            
            C_DCOPY(norbs,Temp,1,&D[0][i],norbs);            

            Temp_p = P[i];
            P[i] = P[pivot];
            P[pivot] = Temp_p;
        }

        //fprintf(outfile,"  D (pivoted):\n");
        //print_mat(D,norbs,norbs,outfile);        
        //for (int i = 0; i<norbs; i++)
        //    fprintf(outfile,"  Pivot %d is %d.\n",i,P[i]); 

        //Cholesky Decomposition
        int status = C_DPOTRF('U',norbs,D[0],norbs);
        if (status < 0) {
            fprintf(outfile,"  Cholesky Decomposition Failed");
            fflush(outfile);
            exit(PSI_RETURN_FAILURE);
        }
        for (int i = 0; i<norbs-1; i++)
            for (int j = i+1; j<norbs; j++)
                D[i][j] = 0.0; 
 
        //fprintf(outfile,"  C Guess (Cholesky Unpivoted):\n");
        //print_mat(D,norbs,ndocc,outfile);        
        
        //Unpivot
        double** C = block_matrix(norbs,ndocc);
        for (int m = 0; m < norbs; m++) {
            C_DCOPY(ndocc,&D[m][0],1,&C[P[m]][0],1);
        }

        //Set C
        C_->zero();
        doccpi_[0] = nalpha_;
        for (int m = 0; m<norbs; m++)
            for (int i = 0; i<ndocc; i++)
                C_->set(0,m,i,C[m][i]);


        //fprintf(outfile,"  C Guess (Cholesky):\n");
        //print_mat(C,norbs,ndocc,outfile);        

        free(P);
        free(Temp);
        free_block(D);
        free_block(C);

        E_ = 0.0; //For now

        timer_off("SAD Cholesky");

    } else {
        throw std::invalid_argument("ID or CHOLESKY are the only two C matrix generators at the moment");
    }
    /**
    //Compute a rough Fock matrix via integral direct
    //Note, my convention is backwards, l and s are index variables
    // m and n are zip variables
    timer_on("SAD Fock");
    double SAD_Schwarz = options_.get_double("SAD_SCHWARZ_CUTOFF");
    G_->zero();
    J_->zero();
    K_->zero();
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
    TwoBodyInt *TEI = integral.eri(0, SAD_Schwarz);
    const double* buffer = TEI->buffer();
    
    long kshells = 0;
    long kdone = 0;
    long jshells = 0;
    long jdone = 0;

    //J
    for (int A = 0, shell_offset = 0; A<mol->nallatom(); A++) {
    for (int MU = shell_offset; MU < shell_offset+atomic_bases[A]->nshell(); MU++) {
    int numMU = basisset_->shell(MU)->nfunction();
    for (int NU = shell_offset; NU <= MU ; NU++) {
    int numNU = basisset_->shell(NU)->nfunction();
    for (int LA = 0; LA < basisset_->nshell(); LA++) {
    int numLA = basisset_->shell(LA)->nfunction();
    for (int SI = 0; SI <= LA; SI++) {
    int numSI = basisset_->shell(SI)->nfunction();

    jshells++;

    if (!TEI->shell_is_zero(MU,NU,LA,SI)) {
    jdone++;
    TEI->compute_shell(MU,NU,LA,SI);
    for (int m = 0, index = 0; m<numMU; m++) {
    int omu = basisset_->shell(MU)->function_index() + m;
    for (int n = 0; n<numNU; n++) {
    int onu = basisset_->shell(NU)->function_index() + n;
    for (int l = 0; l<numLA; l++) {
    int ola = basisset_->shell(LA)->function_index() + l;
    for (int s = 0; s<numSI; s++,index++) {
    int osi = basisset_->shell(SI)->function_index() + s;
        //J_->add(0,ola,osi,2.0*D_->get(0,omu,onu)*buffer[index]); //Redundant shells
        if (onu == omu) {
            if (osi == ola) {
                J_->add(0,ola,osi,2.0*D_->get(0,omu,onu)*buffer[index]);
            } else if (osi < ola) {
                J_->add(0,ola,osi,2.0*D_->get(0,omu,onu)*buffer[index]);
                J_->add(0,osi,ola,2.0*D_->get(0,omu,onu)*buffer[index]);
            }
        }
        else if (onu < omu) {
            if (osi == ola) {
                J_->add(0,ola,osi,4.0*D_->get(0,omu,onu)*buffer[index]);
            } else if (osi < ola) {
                J_->add(0,ola,osi,4.0*D_->get(0,omu,onu)*buffer[index]);
                J_->add(0,osi,ola,4.0*D_->get(0,omu,onu)*buffer[index]);
            }
        }  
    }
    }
    }
    }
    }


    }
    }
    }
    }
    shell_offset+= atomic_bases[A]->nshell();
    }
    //End J
    
    //K
    for (int A = 0, shell_offset = 0; A<mol->nallatom(); A++) {
    for (int LA = 0; LA < basisset_->nshell(); LA++) {
    int numLA = basisset_->shell(LA)->nfunction();
    for (int MU = shell_offset; MU < shell_offset+atomic_bases[A]->nshell() && MU <= LA; MU++) {
    int numMU = basisset_->shell(MU)->nfunction();
    for (int SI = 0; SI <= LA; SI++) {
    int numSI = basisset_->shell(SI)->nfunction();
    for (int NU = shell_offset; NU < shell_offset+atomic_bases[A]->nshell() && NU <= SI; NU++) {
    int numNU = basisset_->shell(NU)->nfunction();

    kshells++;

    if (!TEI->shell_is_zero(LA,MU,SI,NU)) {
    kdone++;
    TEI->compute_shell(LA,MU,SI,NU);
    for (int l = 0, index = 0; l<numLA; l++) {
    int ola = basisset_->shell(LA)->function_index() + l;
    for (int m = 0; m<numMU; m++) {
    int omu = basisset_->shell(MU)->function_index() + m;
    for (int s = 0; s<numSI; s++) {
    int osi = basisset_->shell(SI)->function_index() + s;
    for (int n = 0; n<numNU; n++, index++) {
    int onu = basisset_->shell(NU)->function_index() + n;
        if (ola == osi && omu == onu) {
            K_->add(0,ola,osi,-D_->get(0,omu,onu)*buffer[index]);
        } else if (ola == osi && omu > onu) {
            K_->add(0,ola,osi,-2.0*D_->get(0,omu,onu)*buffer[index]);       
        } else if (ola > osi && omu == onu) {
            K_->add(0,ola,osi,-D_->get(0,omu,onu)*buffer[index]);       
            K_->add(0,osi,ola,-D_->get(0,omu,onu)*buffer[index]);       
        } else if (ola > osi && omu > onu) {
            K_->add(0,ola,osi,-2.0*D_->get(0,omu,onu)*buffer[index]);       
            K_->add(0,osi,ola,-2.0*D_->get(0,omu,onu)*buffer[index]);        
        }
    }
    }
    }
    }
    }

    }
    }
    }
    }
    shell_offset+= atomic_bases[A]->nshell();
    }
    //End K
    timer_off("SAD Fock");


    if (print_>1) {
        fprintf(outfile,"  SAD Fock Computation:\n");
        fprintf(outfile,"  %ld of %ld shells computed for Coulomb.\n",jdone,jshells);
        fprintf(outfile,"  %ld of %ld shells computed for Exchange.\n",kdone,kshells);
    }

    G_->copy(J_);
    G_->add(K_);

    J_->print(outfile);
    K_->print(outfile);

    F_->copy(H_);
    F_->add(G_);

    //Compute initial E for reference
    E_ = compute_E();

    //Form the C matrix from the rough Fock matrix
    form_C();
    //Form the D matrix from the resultant C matrix
    form_D();

    if (print_>7) {
        G_->print(outfile);
        F_->print(outfile);
        C_->print(outfile);
        D_->print(outfile);
    }

    **/
    /**
    timer_on("SAD C");


    //Form the C matrix in reverse via eigendecoposition of D.
    C_->zero();
    int norbs = nso_;
    int ndocc = nalpha_;
    doccpi_[0] = nalpha_;
    //fprintf(outfile,"  Ndocc %d\n",ndocc);
    double **D = D_->to_block_matrix();
    double **V = block_matrix(norbs,norbs);    

    //print_mat(D,norbs,norbs,outfile);

    double *eigvals = init_array(norbs);
    sq_rsp(norbs, norbs, D,  eigvals, 1, V, 1.0e-14);

    for (int i = ndocc; i < norbs; i++)
        fprintf(outfile,"  Eigval %d, %14.10f\n",i,eigvals[i]);

    print_mat(V,norbs,norbs,outfile);
    //Yep
    for (int i = norbs-ndocc-1; i<norbs; i++)
        C_DSCAL(norbs, sqrt(eigvals[i]), V[i], 1);

    print_mat(V,norbs,norbs,outfile);

    for (int m = 0; m<norbs; m++)
        for (int i = 0; i<ndocc; i++)
            C_->set(0,m,i,V[m][norbs-i-1]);
    C_->print(outfile);

    //Frees
    free(D);    
    free(V);    

    timer_off("SAD C");
    **/
    for (int A = 0; A<mol->nallatom(); A++)
        free_block(atomic_D[A]);
    free(atomic_D);

    free(nelec);
    free(nhigh);
    free(nalpha);
    free(nbeta);
    //abort();
}
}}
