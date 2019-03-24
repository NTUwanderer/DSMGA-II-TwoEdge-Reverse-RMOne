/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#define MAX_GEN 200

int step = 30;

using namespace std;

struct Record {
    int n;
    double nfe;
    double gen;
    double nfe_std;

    int rmSuccess = 0;
    int rmFail = 0;
    int bmSuccess = 0;
    int bmFail = 0;
    double rmRate = 0;
    double bmRate = 0;
    int successCnt = 0;
};

void initCnt(Record& rec) {
    rec.rmSuccess = 0;
    rec.rmFail = 0;
    rec.bmSuccess = 0;
    rec.bmFail = 0;
    rec.rmRate = 0;
    rec.bmRate = 0;
    rec.successCnt = 0;
}

void addGaCnt(Record& rec, const DSMGA2& ga) {
    rec.rmSuccess += ga.rmSuccess;
    rec.rmFail += ga.rmFail;
    rec.bmSuccess += ga.bmSuccess;
    rec.bmFail += ga.bmFail;
    rec.rmRate += 1.0 * ga.rmSuccess / max(ga.rmFail+ga.rmSuccess, 1);
    rec.bmRate += 1.0 * ga.bmSuccess / max(ga.bmFail+ga.bmSuccess, 1);
    rec.successCnt++;
}

bool parseLine(string line, int popu, bool& foundOptima, double& nfe) {
    int pos = line.find("]");
    if (pos == -1)
        return false;
    if (to_string(popu) != line.substr(1, pos - 1))
        return false;

    int firstColon = line.find(":");
    if (firstColon == -1)
        return false;
    int secondColon = line.find(":", firstColon + 2);
    if (secondColon == -1)
        return false;

    foundOptima = line[secondColon - 2];
    nfe = stod(line.substr(secondColon + 2));

    return true;
}

int main (int argc, char *argv[]) {

    if (argc != 4 && argc!=5 && argc !=6 && argc != 7) {
        printf ("sweep ell numConvergence function(0~3, 7, 8)\n");
        printf ("sweep ell numConvergence 4 [step #] [nk problem #]\n");
        printf ("sweep ell numConvergence 5 [spin problem #]\n");
        printf ("sweep ell numConvergence 6 [sat problem #]\n");
        printf ("sweep ell numConvergence 9 [mkp problem #]\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");
        printf ("     L_FTRAP  : 7\n");
        printf ("     L_2FTRAP : 8\n");
        printf ("     MKP   :  9\n");
        return -1;
    }

    int ell = atoi (argv[1]);
    int numConvergence = atoi (argv[2]); // problem size
    int fffff = atoi(argv[3]);

    int problemNum = 0;
    int neighborNum = 0;
    int stepNum = 0;

    printf (" \n");

    if (fffff == 4) {
        neighborNum = 4;
        stepNum = atoi (argv[4]);
        problemNum = atoi (argv[5]);
    }

    if (fffff == 5 || fffff == 6 || fffff == 9) {
        problemNum = atoi (argv[4]);
    }


    int nInitial = 10;


    // for debug
    // myRand.seed(2);


    Statistics st;

    Statistics stGen, stLS, stNFE;


    string line;
    ifstream record;
    bool failLoad = false;
    if (fffff == 5) {
	char filename[200];
        sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSPIN(filename, &mySpinGlassParams);

        if (argc > 5) {
            record.open(argv[5], ifstream::in);
            while (record.good()) {
                getline(record, line);
                if (line == "Bisection phase 1")
                    break;
            }
        } else {
            failLoad = true;
        }
    }

    if (fffff == 4) {
        char filename[200];
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, neighborNum, stepNum, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (fffff == 6) {
        char filename[200];
        sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf",ell,ell,problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);

        if (argc > 5) {
            record.open(argv[5], ifstream::in);
            while (true) {
                getline(record, line);
                if (line == "Bisection phase 1")
                    break;
            }
        }
    }
    
    if (fffff == 9) {
        char filename[200];
        sprintf(filename, "./All-MKP-Instances/sac94/weish/weish%02d.dat", problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadMKP(filename, &myMKP);
        ell = myMKP.var;
    }
    
    bool foundOptima;
    Record rec[3];
    rec[0].n = nInitial;
    rec[1].n = nInitial+step;
    rec[2].n = nInitial+step+step;

    int popu;
    Record q1, q3;

    if (SHOW_BISECTION && (!record.is_open()) && !failLoad) printf("Bisection phase 1\n");
    fflush(NULL);

    for (int i=0; i<3; ++i) {
        popu = rec[i].n;

        initCnt(rec[i]);

        bool loaded = false;
        if (record.is_open() && !failLoad) {
            getline(record, line);

            loaded = parseLine(line, popu, foundOptima, rec[i].nfe);
        } 
        if (!loaded) {
            failLoad = true;
            if (SHOW_BISECTION) printf("[%d]: ", popu);

            foundOptima = true;

            stGen.reset();
            stNFE.reset();
            stLS.reset();

            for (int j=0; j<numConvergence; j++) {

                DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
                ga.doIt(false);

                stGen.record(ga.getGeneration());
                stNFE.record(Chromosome::hitnfe);
                stLS.record(Chromosome::lsnfe);


                if (!ga.foundOptima()) {

                    foundOptima = false;

                    if (SHOW_BISECTION) {
                        printf("-");
                        fflush(NULL);
                    }
                    break;
                }

                addGaCnt(rec[i], ga);

                if (SHOW_BISECTION) {
                    printf("+");
                    fflush(NULL);
                }
            }


            rec[i].gen = stGen.getMean();

            if (!foundOptima)
                rec[i].nfe = INF;
            else {
                rec[i].nfe = stNFE.getMean();
                rec[i].nfe_std = stNFE.getStdev();
            }

            if (SHOW_BISECTION) printf(" : %f \n", rec[i].nfe);
        }

        if (!foundOptima)
            rec[i].nfe = INF;
    }

    while (rec[0].nfe < rec[1].nfe  && ((rec[2].n-rec[0].n)*20 > rec[1].n)) {

        rec[2] = rec[1];
        rec[1].n = (rec[0].n + rec[2].n) / 2;
        step /= 2;
        popu = rec[1].n;

        initCnt(rec[1]);

        if (SHOW_BISECTION) printf("[%d]: ", popu);

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt(false);

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::hitnfe);
            stLS.record(Chromosome::lsnfe);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }
            addGaCnt(rec[1], ga);

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }


        rec[1].gen = stGen.getMean();

        if (!foundOptima)
            rec[1].nfe = INF;
        else {
            rec[1].nfe = stNFE.getMean();
            rec[1].nfe_std = stNFE.getStdev();
        }

        if (SHOW_BISECTION) printf(" : %f \n", rec[1].nfe);
    }


    while ( (rec[1].nfe >= rec[0].nfe) || (rec[1].nfe >= rec[2].nfe)) {

        popu = rec[2].n + step;

        rec[0] = rec[1];
        rec[1] = rec[2];
        rec[2].n = popu;

        initCnt(rec[2]);

        bool loaded = false;
        if (record.is_open() && !failLoad) {
            getline(record, line);

            loaded = parseLine(line, popu, foundOptima, rec[2].nfe);
        } 

        if (!loaded) {
            failLoad = true;

            if (SHOW_BISECTION) printf("[%d]: ", popu);

            foundOptima = true;

            stGen.reset();
            stNFE.reset();
            stLS.reset();

            for (int j=0; j<numConvergence; j++) {

                DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
                ga.doIt(false);

                stGen.record(ga.getGeneration());
                stNFE.record(Chromosome::hitnfe);
                stLS.record(Chromosome::lsnfe);


                if (!ga.foundOptima()) {

                    foundOptima = false;

                    if (SHOW_BISECTION) {
                        printf("-");
                        fflush(NULL);
                    }
                    break;
                }
                addGaCnt(rec[2], ga);

                if (SHOW_BISECTION) {
                    printf("+");
                    fflush(NULL);
                }
            }

            rec[2].gen = stGen.getMean();

            if (!foundOptima)
                rec[2].nfe = INF;
            else {
                rec[2].nfe = stNFE.getMean();
                rec[2].nfe_std = stNFE.getStdev();
            }

            if (SHOW_BISECTION) printf(" : %f \n", rec[2].nfe);
        }


        if (!foundOptima)
            rec[2].nfe = INF;
    }


    if (record.is_open())
        getline(record, line);

    if (SHOW_BISECTION && (!record.is_open() || !record.good())) printf("Bisection phase 2\n");

    while ( ((rec[2].n-rec[0].n)*20 > rec[1].n) && (rec[2].n>rec[1].n+1) && (rec[1].n>rec[0].n+1)) {

        q1.n = (rec[0].n + rec[1].n) / 2;
        initCnt(q1);

        bool loaded = false;
        if (record.is_open() && !failLoad) {
            getline(record, line);

            loaded = parseLine(line, q1.n, foundOptima, q1.nfe);
        } 

        if (!loaded) {
            failLoad = true;

            if (SHOW_BISECTION) printf("[%d]: ", q1.n);

            foundOptima = true;

            for (int j=0; j<numConvergence; j++) {

                DSMGA2 ga(ell, q1.n, MAX_GEN, -1, fffff);
                ga.doIt(false);

                if (!ga.foundOptima()) {
                    foundOptima = false;
                    if (SHOW_BISECTION) {
                        printf("-");
                        fflush(NULL);
                    }
                    break;
                }
                addGaCnt(q1, ga);

                if (SHOW_BISECTION) {
                    printf("+");
                    fflush(NULL);
                }
                if (j==0) {
                    stGen.reset();
                    stLS.reset();
                    stNFE.reset();
                }
                stGen.record(ga.getGeneration());
                stNFE.record(Chromosome::hitnfe);
                stLS.record(Chromosome::lsnfe);
            }

            q1.gen = stGen.getMean();
            if (foundOptima) {
                q1.nfe = stNFE.getMean();
                q1.nfe_std = stNFE.getStdev();
            }
            else
                q1.nfe = INF;


            if (SHOW_BISECTION) printf(" : %f \n", q1.nfe);

        }

        if (!foundOptima)
            q1.nfe = INF;

        q3.n = (rec[1].n + rec[2].n) / 2;
        initCnt(q3);

        loaded = false;
        if (record.is_open() && !failLoad) {
            getline(record, line);

            loaded = parseLine(line, q3.n, foundOptima, q3.nfe);
        } 

        if (!loaded) {
            failLoad = true;

            if (SHOW_BISECTION) printf("[%d]: ", q3.n);

            foundOptima = true;

            for (int j=0; j<numConvergence; j++) {

                DSMGA2 ga(ell, q3.n, MAX_GEN, -1, fffff);
                ga.doIt(false);

                if (!ga.foundOptima()) {
                    foundOptima = false;
                    if (SHOW_BISECTION) {
                        printf("-");
                        fflush(NULL);
                    }
                    break;
                }
                addGaCnt(q3, ga);

                if (SHOW_BISECTION) {
                    printf("+");
                    fflush(NULL);
                }
                if (j==0) {
                    stGen.reset();
                    stLS.reset();
                    stNFE.reset();
                }
                stGen.record(ga.getGeneration());
                stNFE.record(Chromosome::hitnfe);
                stLS.record(Chromosome::lsnfe);
            }

            q3.gen = stGen.getMean();
            if (foundOptima) {
                q3.nfe = stNFE.getMean();
                q3.nfe_std = stNFE.getStdev();
            } else
                q3.nfe = INF;

            if (SHOW_BISECTION) printf(" : %f \n", q3.nfe);
        }
        if (!foundOptima)
            q3.nfe = INF;

        if (rec[1].nfe < q1.nfe && rec[1].nfe < q3.nfe) {
            rec[0] = q1;
            rec[2] = q3;
        } else if (q1.nfe < rec[1].nfe && q1.nfe < q3.nfe) {
            rec[2] = rec[1];
            rec[1] = q1;
        } else { // q3nfe smallest
            rec[0] = rec[1];
            rec[1] = q3;
        }
    };



    if (fffff == 4)
        freeNKWAProblem(&nkwa);

    int rmSuccess = rec[1].rmSuccess;
    int rmFail = rec[1].rmFail;
    int bmSuccess = rec[1].bmSuccess;
    int bmFail = rec[1].bmFail;
    double rmRate = rec[1].rmRate / rec[1].successCnt;
    double bmRate = rec[1].bmRate / rec[1].successCnt;

    printf("population: %d\n", rec[1].n);
    printf("generation: %f\n", rec[1].gen);
    printf ("RM_Success: %i %i %.5f %%\n", rmSuccess, rmFail, (100.0 * rmRate));
    printf ("BM_Success: %i %i %.5f %%\n", bmSuccess, bmFail, (100.0 * bmRate));
    printf("NFE: %f\n", rec[1].nfe);
    printf("N_std: %f\n", rec[1].nfe_std);


    return EXIT_SUCCESS;

}

