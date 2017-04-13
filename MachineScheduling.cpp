#ifndef MACHINE_SCHEDULING_CPP
#define MACHINE_SCHEDULING_CPP

#define ZBASED_1D_POS(r,c,num_cols) (r*num_cols+c)
#define GLPK_1D_POS(r,c,num_cols) ((r-1)*num_cols+c)

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <glpk.h>

int i, j, size;
std::vector<int> times;
std::vector<double> schedule;

// reads a well formatted file into grid
void loadFromFile(char* filename) {
    std::ifstream infile;
    infile.open(filename);

    if (infile.fail())
        throw std::invalid_argument("File cannot be read.");

    infile >> j; //number jobs     -rows
    infile >> i; //number machines -cols
    size = i*j;
    char ch;

    std::vector<int> temp(size);

    int c = 0;
    for (int r = 0; r < j; ++r) {
        for (c = 0; c < i - 1; ++c) {
            infile >> temp.at(ZBASED_1D_POS(r, c, i));
            infile >> ch; //eat commas
        }
        infile >> temp.at(ZBASED_1D_POS(r, c, i));
    }

    //transpose the matrix to make things easier later
    for (int c = 0; c < i; ++c)
        for (int r = 0; r < j; ++r)
            times.push_back(temp.at(ZBASED_1D_POS(r, c, i)));

    infile.close();
};



void findOpt() {
    glp_prob *lp;
      /*
      MATRIX:

      minimize t subject to

      x11 + x12 +      ...              +0xij + s1 + ...      - t = 0 /
      .                                           .                 / i equations representing
      .                                           .                 / sum of xij <= t
      0  +     +      ...              + xij +      ... + si - t = 0 /

      x11 +     + ... + x21 +0x22 + ... +0xij                     = 1 /
      .                                           .                 / j equations representing
      .                                           .                 / sum of xij =1
      0x11 + ... + x1j +     ...         + xij                     = 1

      */
    
    int ia[40000], ja[40000];
    double ar[40000];
    lp = glp_create_prob();

    glp_add_rows(lp, i + j);

    for (int ii = 1; ii <= i; ++ii)
        glp_set_row_bnds(lp, ii, GLP_FX, 0, 0);

    for (int jj = 1; jj <= j; ++jj)
        glp_set_row_bnds(lp, i + jj, GLP_FX, 1, 1);

    int numCols = size + i + 1;
    glp_add_cols(lp, numCols);

    // slack vars >=0
    for (int n = i*j + 1; n <= numCols - 1; ++n)
        glp_set_col_bnds(lp, n, GLP_LO, 0.0, 0.0);

    // t>=0
    glp_set_col_bnds(lp, numCols, GLP_LO, 0.0, 0.0);

    // xij >= 0
    for (int n = 0; n < size; ++n)
        glp_set_col_bnds(lp, n + 1, GLP_LO, 0.0, 0.0);

    // obj func values all 0 except t
    for (int n = 1; n <= numCols - 1; ++n)
        glp_set_obj_coef(lp, n, 0);
    glp_set_obj_coef(lp, numCols, 1);

    // set all coefficients to 0 first; A is sparse
    for (int r = 1; r <= (i + j); ++r)
        for (int c = 1; c <= numCols; ++c)
        {
        int p = GLPK_1D_POS(r, c, numCols);
        ia[p] = r, ja[p] = c, ar[p] = 0;
        }

    // xij <= t
    int c = 1;
    for (int r = 1; r <= i; ++r)
        for (int n = 0; n < j; ++n, ++c)
            ar[GLPK_1D_POS(r, c, numCols)] = times.at(ZBASED_1D_POS((r - 1), n, j));

    // xij = 1
    for (int r = i + 1; r <= i + j; ++r) //bottom j rows
        for (int n = 0; n < i; ++n)
            ar[GLPK_1D_POS(r, (r - i + j*n), numCols)] = 1;

    // slack variables
    for (int r = 1; r <= i; ++r)
    {
        c = i*j + r;
        ar[GLPK_1D_POS(r, c, numCols)] = 1;
    }

    // subtracted t
    for (int n = 1; n <= i; ++n)
        ar[numCols*n] = -1;

    glp_load_matrix(lp, numCols*(i + j), ia, ja, ar);

    // sort the job times
    std::vector<int> sorted_times = times;
    std::sort(sorted_times.begin(), sorted_times.end());

    int idx = size - 1;
    double T = 0, t = 0, Tprev;
    double min = std::numeric_limits<double>::infinity();
    bool done = false;
    glp_smcp* smcp = new glp_smcp();
    glp_init_smcp(smcp);
    smcp->msg_lev = GLP_MSG_OFF;

    glp_prob *oldlp;
    oldlp = glp_create_prob();
    glp_copy_prob(oldlp, lp, GLP_OFF);

    do {
        // cuts unneccessary simplexes if there are multiple edges with the same weight
        while (t == T && idx != -1) {
            t = sorted_times.at(idx);
            --idx;
        }

        T = t;

        // xij = 0 if dij > T
        for (int n = 0; n < size; ++n)
            if (times.at(n) > T)
                glp_set_col_bnds(lp, n + 1, GLP_FX, 0.0, 0.0);

        if (glp_simplex(lp, smcp) == 0 && glp_get_status(lp) == GLP_OPT) {
            double val = glp_get_obj_val(lp) + T;

            if (val <= min) {
                // if the most recent solution is good, save it to oldlp
                // in case the next solution is infeasible or worse
                glp_copy_prob(oldlp, lp, GLP_OFF);
                min = val;
            }
        }
        else
            done = true;

    } while (!done && idx != -1);

    for (int n = 1; n <= size; ++n)
        schedule.push_back(glp_get_col_prim(oldlp, n));

    glp_delete_prob(lp);
    glp_delete_prob(oldlp);
};


void matchG()
{
    // these vectors are used to keep track of our graph
    // all are 0 based
    // a 'true' in position n means that machine or job n is included
    std::vector<bool> machines(i, false);
    std::vector<bool> jobs(j, false);
    // edges contain a tuple of their weight, their machine #, and their job #
    // edges' machines & jobs are enumerated starting with 1
    std::vector<std::vector<double> > edges;

    for (int n = 0; n < size; ++n) {
        double v = schedule.at(n);
        if (v != 1.0 && v != 0.0) {
            int machine = n / j;
            machines.at(machine) = true;

            int job = n - machine*j;
            jobs.at(job) = true;

            // since we are only approximating the optimal solution,
            // minimizing times is a naive attempt to push the solution in the right direction

            edges.push_back(std::vector < double > {1.0 / (double)times.at(n), (double)machine + 1.0, (double)job + 1});
        }
    }

    int numEdges = edges.size();

    if (numEdges != 0) {
        int numMachines = 0, numJobs = 0;
        for (int n = 0; n < machines.size(); ++n)
            if (machines.at(n) == true)
                ++numMachines;

        for (int n = 0; n < jobs.size(); ++n)
            if (jobs.at(n) == true)
                ++numJobs;

        int numVerticies = numMachines + numJobs;
        // there is a column for each edge and for each slack variable
        // there are numVerticies equations and that many slack variables
        int numCols = numVerticies + numEdges;

        glp_prob *lp;

        int ia[10000], ja[10000];
        double ar[10000];
        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MAX);

        glp_add_rows(lp, numVerticies);

        // sum of all edges into all nodes is equal to 1
        for (int n = 1; n <= numVerticies; ++n)
            glp_set_row_bnds(lp, n, GLP_FX, 1, 1);

        glp_add_cols(lp, numCols);

        // all edges have constraints >=0
        for (int n = 1; n <= numEdges; ++n)
            glp_set_col_bnds(lp, n, GLP_LO, 0, 0);

        // slack vars have constraints >=0
        for (int n = numEdges; n <= numCols; ++n)
            glp_set_col_bnds(lp, n, GLP_LO, 0, 0);

        // obj func values for edges
        for (int n = 0; n < numEdges; ++n)
            glp_set_obj_coef(lp, n + 1, edges.at(n).at(0));

        // obj func values for slack variables
        for (int n = numEdges + 1; n <= numCols; ++n)
            glp_set_obj_coef(lp, n, 0);

        // set all vals to 0 first; A is sparse
        for (int r = 1; r <= numVerticies; ++r)
            for (int c = 1; c <= numCols; ++c)
            {
            int p = GLPK_1D_POS(r, c, numCols);
            ia[p] = r, ja[p] = c, ar[p] = 0;
            }

        // xij <= t
        int m, machineNo = 1, idx = 0;
        for (int r = 1; r <= numMachines; ++r) {
            // find all edges that share a machine vertex
            machineNo = edges.at(idx)[1];
            for (; idx < edges.size(); ++idx) {
                m = edges.at(idx)[1];
                if (m != machineNo)
                    break;
                else
                    ar[GLPK_1D_POS(r, (idx + 1), numCols)] = 1;
            }
        }

        int jobNo = 0;
        // xij = 1
        for (int r = numMachines + 1; r <= numVerticies; ++r) { //bottom j rows
            ++jobNo;
            // find all edges with the same job number
            for (int idx = 0; idx < numEdges; ++idx)
                if (edges.at(idx)[2] == jobNo)
                    ar[GLPK_1D_POS(r, (idx + 1), numCols)] = 1;
        }

        // slack variables
        int c = numEdges + 1;
        for (int r = 1; r <= numVerticies; ++r) {
            ar[GLPK_1D_POS(r, c, numCols)] = 1;
            c++;
        }

        glp_load_matrix(lp, numVerticies*numCols, ia, ja, ar);

        glp_smcp* smcp = new glp_smcp();
        glp_init_smcp(smcp);
        smcp->msg_lev = GLP_MSG_OFF;

        if (glp_simplex(lp, smcp) != 0 || glp_get_status(lp) != GLP_OPT)
            std::cerr << "ERROR SOMETHING IS BAD";

        idx = 0;
        for (int n = 0; n < size; ++n) {
            double v = schedule.at(n);
            if (v != 1.0 && v != 0.0) {
                schedule.at(n) = glp_get_col_prim(lp, idx + 1);
                ++idx;
            }
        }
        glp_delete_prob(lp);
    }
};

void writeSolution(char* filename) {
    std::vector<int> output(j, 0);
    for (int n = 0; n < size; ++n)
        if (schedule.at(n) == 1)
            output.at(n%j) = n / j;

    std::ofstream outfile;
    outfile.open(filename);
    outfile << (output.at(0) + 1);
    for (int n = 1; n < j; ++n)
        outfile << "," << (output.at(n) + 1);

    outfile.close();
};

int main(int argc, char*argv[]) {
    if (argc != 3)
        throw std::invalid_argument("Incorrect amount of parameters passed to main.");

    loadFromFile(argv[1]);
    findOpt();
    matchG();
    writeSolution(argv[2]);
}

#endif
