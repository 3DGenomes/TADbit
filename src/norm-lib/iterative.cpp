#include "iterative.h"
using namespace std;

struct Files
{
    string HiC_in;
    string HiC_iterative_fo;
};

struct Computation
{
    int nbin_low_bound; // lower bound of bin to consider (default: 1)

    bool drop_bin; // drop bins according to frac_low_bound at first step (default: 1)
    double freq_low_bound; // fraction of bins to throw away depending on their total number of contacts

    int nloop_init;
    int nloop_max;

    bool t_chrom; // iterative analysis constrained to a specific chromosome (default: 0)
    string chrom; // iterative analysis constrained to a specific chromosome (default: 0)    
    
    bool B; // compute and write out the biases
    
    string chrom_to_rw;
};

int compare_dec_a1b1(const double ** a, const double ** b)
{
    if ((*a)[1] - (*b)[1] < 0) return 1;
    else return -1;
}

void read_arg(int argc, const char *argv[], Files *files, Computation* comp);
void init_W(map<string, map<int, map<string, map<int, double> > > > *W, map<string, map<int, map<string, map<int, bool> > > > *pairW, Files *file, Computation *comp);
void drop_the_bins(Computation *comp, map<string, map<int, map<string, map<int, double> > > > *W,
        map<string, map<int, map<string, map<int, bool> > > > *pairW);
void init_B(map<string, map<int, double> > *B, map<string, map<int, map<string, map<int, double> > > > &W);
void update_S(map<string, map<int, double> > *S, double *meanS, map<string, map<int, map<string, map<int, double> > > > &W);
void update_DB(map<string, map<int, double> > *DB, map<string, map<int, double> > &S, double meanS);
void update_W(map<string, map<int, map<string, map<int, double> > > > *W, map<string, map<int, double> > &DB);
void update_B(map<string, map<int, double> > *B, map<string, map<int, double> > &DB);
// print_out_W: update meanTotalW for the nomalisation of B: mtW
void print_out_W(map<string, map<int, map<string, map<int, double> > > > &W, double *mtW,
        map<string, map<int, map<string, map<int, bool> > > > &pairW, Files *file, int nloop);
void print_out_B(map<string, map<int, double> > &B, double mtW, Files *file, int nloop);

int main(int argc, const char *argv[])
{
    Files *files = new Files;
    Computation *comp = new Computation;

    read_arg(argc, argv, files, comp);
    
    map<string, map<int, map<string, map<int, double> > > > W; // Imakaev annotation ([chrom1][bin1][chrom2][bin2] => freq)
    map<string, map<int, map<string, map<int, bool> > > > pairW; // test only half pairs
    init_W(&W, &pairW, files, comp);
    
    // drop the bins with too few contacts        
    if (comp->drop_bin) drop_the_bins(comp, &W, &pairW); 
    
    // Biases
    map<string, map<int, double> > B;
    init_B(&B, W);
    
    // mean total weight for the normalisation of B
    double mtW = 0;
    for (int nloop = comp->nloop_init; nloop <= comp->nloop_max; nloop++) {
        print_out_W(W, &mtW,pairW, files, nloop);
        if(comp->B) print_out_B(B,mtW,files, nloop);
        
        // Coverages
        map<string, map<int, double> > S;
        double meanS;
        update_S(&S, &meanS, W);
        
        // Additional Biases
        map<string, map<int, double> > DB; 
        update_DB(&DB, S, meanS);
        
        if(comp->B) update_B(&B, DB);
        update_W(&W, DB);
    }
    
    return 0;
}

void read_arg(int argc, const char *argv[], Files *files, Computation* comp)
{
     // default values
    comp->freq_low_bound = 0.02; // 2%: Ismakaev
    comp->drop_bin = 1;
    comp->t_chrom = 0;
    comp->nloop_init = 0;
    comp->nloop_max = 20;
    comp->B = 0;
    
    for (int i = 0; i < argc; i++) {
        // <editor-fold defaultstate="collapsed" desc="help command">
        if (!strcmp(argv[i], "-h")) {
            cout << endl << "-i [input file]" << endl;
            cout << "structure of the input file:" << endl;
            cout << "chrom1\tbin1\tchrom2\tnin2\tHiC_bin" << endl;
            cout << "18\t2\t18\t4\t0.23" << endl;
            cout << "!! check the redondance pb (forgotten one) !!"<< endl << endl;

            cout << "-iterative [output folder]" << endl;
            cout << "output: initial input file + corresponding file after n iterations" << endl;            
            cout << "[optional] -nmax: number of loops" <<endl;
            cout << "[optional] -B: returns the biases" <<endl << endl;
            exit(0);
        }// </editor-fold>
        
        if (!strcmp(argv[i], "-i")) files->HiC_in = argv[i + 1];
        if (!strcmp(argv[i], "-n0")) comp->nloop_init = atoi(argv[i + 1]);        
        if (!strcmp(argv[i], "-nmax")) comp->nloop_max = atoi(argv[i + 1]);       
        if (!strcmp(argv[i], "-iterative")){
            files->HiC_iterative_fo = argv[i + 1];
            // mkdir(files->HiC_iterative_fo);
        }
        if (!strcmp(argv[i], "-freq_low")) comp->freq_low_bound = atof(argv[i + 1]);
        if (!strcmp(argv[i], "-no_drop")) comp->drop_bin = 0;
        if (!strcmp(argv[i], "-chrom")) {
            comp->chrom = argv[i + 1];
            comp->t_chrom = 1;
        }
        if (!strcmp(argv[i], "-B")) comp->B = 1;
    }
    return;
}

void init_W(map<string, map<int, map<string, map<int, double> > > > *W, map<string, map<int, map<string, map<int, bool> > > > *pairW,
        Files *file, Computation *comp)
{
    ifstream ifstr(file->HiC_in.c_str());
    string line;
    getline(ifstr, line); // headers
    int nline = 1;
    while (getline(ifstr, line)) {
        if (!(nline % 10000)) {
            cout << "Downloading line " << nline << "             \r";
            cout.flush();
        }        
        vector<string> terms;
        boost::split(terms, line, boost::is_any_of("\t"));
        if (!comp->t_chrom || ((terms[0] == comp->chrom)&&(terms[2] == comp->chrom))) {
            (*W)[terms[0]][atoi(terms[1].c_str())][terms[2]][atoi(terms[3].c_str())] = atof(terms[4].c_str());
            (*pairW)[terms[0]][atoi(terms[1].c_str())][terms[2]][atoi(terms[3].c_str())] = 1;
            (*W)[terms[2]][atoi(terms[3].c_str())][terms[0]][atoi(terms[1].c_str())] = atof(terms[4].c_str());
        }
        nline++;
    }
    ifstr.close();
    return;
}

void drop_the_bins(Computation *comp, map<string, map<int, map<string, map<int, double> > > > *W,
        map<string, map<int, map<string, map<int, bool> > > > *pairW)
{
    map<string, map<int, double> > contact_nb; // [chrom][bin] =< contact nb
    vector<double> list_contact_nb;
    vector<string> list_chrom; // chrom correspondant
    vector<int> list_bin; // bin correspondant

    for (map<string, map<int, map<string, map<int, double> > > >::iterator itW = (*W).begin(); itW != (*W).end(); itW++) {
        string chrom1 = (*itW).first;
        cout << "Droping the bins chrom " << chrom1 << "          \r";
        cout.flush();
        for (map<int, map<string, map<int, double> > >::iterator itC1 = (*W)[chrom1].begin(); itC1 != (*W)[chrom1].end(); itC1++) {
            int bin1 = (*itC1).first;
            contact_nb[chrom1][bin1] = 0;

            for (map<string, map<int, double> >::iterator itCb1 = (*W)[chrom1][bin1].begin(); itCb1 != (*W)[chrom1][bin1].end(); itCb1++) {
                string chrom2 = (*itCb1).first;
                for (map<int, double>::iterator itC2 = (*W)[chrom1][bin1][chrom2].begin(); itC2 != (*W)[chrom1][bin1][chrom2].end(); itC2++) {
                    contact_nb[chrom1][bin1] += (*itC2).second;
                }
            }
            list_contact_nb.push_back(contact_nb[chrom1][bin1]);
            list_chrom.push_back(chrom1);
            list_bin.push_back(bin1);
        }
    }

    cout << "Sorting init...                            \r";
    cout.flush();

    long long nlist = list_contact_nb.size();
    double **to_be_sorted = new double*[nlist];
    for (long long i = 0; i < nlist; i++) {
        to_be_sorted[i] = new double[2];
        to_be_sorted[i][0] = (double) i;
        to_be_sorted[i][1] = list_contact_nb[i];
    }

    cout << "Sorting core of a list of size " << nlist << "                            \r";
    cout.flush();
    qsort(to_be_sorted, nlist, sizeof (double), (int(*)(const void*, const void*)) & compare_dec_a1b1);

    int i0 = (int) ((1 - comp->freq_low_bound) * nlist);
    for (long long i = i0; i < nlist; i++) {
        cout << "Removing bins " << i << " of " << nlist << "                           \r";
        cout.flush();

        long long itrue = (long long) to_be_sorted[i][0];
        string chrom = list_chrom[itrue];
        int bin = list_bin[itrue];
        (*W)[chrom].erase(bin);
        if ((*pairW)[chrom].count(bin)) (*pairW)[chrom].erase(bin);
    }

    for (long long i = i0; i < nlist; i++) {
        cout << "Removing (bis) bins " << i << " of " << nlist << "                           \r";
        cout.flush();

        long itrue = (long long) to_be_sorted[i][0];
        string chrom = list_chrom[itrue];
        int bin = list_bin[itrue];

        for (map<string, map<int, map<string, map<int, double> > > >::iterator itW = (*W).begin(); itW != (*W).end(); itW++) {
            string chrom1 = (*itW).first;
            for (map<int, map<string, map<int, double> > >::iterator itC1 = (*W)[chrom1].begin(); itC1 != (*W)[chrom1].end(); itC1++) {
                int bin1 = (*itC1).first;
                if ((*W)[chrom1][bin1].count(chrom) && (*W)[chrom1][bin1][chrom].count(bin)) {
                    (*W)[chrom1][bin1][chrom].erase(bin);
                }

                if ((*pairW)[chrom1].count(bin1) && (*pairW)[chrom1][bin1].count(chrom) && (*pairW)[chrom1][bin1][chrom].count(bin))
                    (*pairW)[chrom1][bin1][chrom].erase(bin);
            }
        }
    }

    return;
}

void init_B(map<string, map<int, double> > *B, map<string, map<int, map<string, map<int, double> > > > &W)
{
    for (map<string, map<int, map<string, map<int, double> > > >::iterator itW = W.begin(); itW != W.end(); itW++) {
        string chrom = (*itW).first;
        cout << "initializing B chrom " << chrom << "          \r";
        cout.flush();

        for (map<int, map<string, map<int, double> > >::iterator itC = W[chrom].begin(); itC != W[chrom].end(); itC++) {
            int bin = (*itC).first;
            (*B)[chrom][bin] = 1;
        }
    }
    return;
}

void update_S(map<string, map<int, double> > *S, double *meanS, map<string, map<int, map<string, map<int, double> > > > &W)
{
    (*meanS) = 0;
    int nstat = 0;
    for (map<string, map<int, map<string, map<int, double> > > >::iterator itW = W.begin(); itW != W.end(); itW++) {
        string chrom1 = (*itW).first;

        cout << "Updating S chrom " << chrom1 << "          \r";
        cout.flush();

        for (map<int, map<string, map<int, double> > >::iterator itC1 = W[chrom1].begin(); itC1 != W[chrom1].end(); itC1++) {
            int bin1 = (*itC1).first;
            (*S)[chrom1][bin1] = 0;

            for (map<string, map<int, double> >::iterator itCb1 = W[chrom1][bin1].begin(); itCb1 != W[chrom1][bin1].end(); itCb1++) {
                string chrom2 = (*itCb1).first;
                for (map<int, double>::iterator itC2 = W[chrom1][bin1][chrom2].begin(); itC2 != W[chrom1][bin1][chrom2].end(); itC2++) {
                    (*S)[chrom1][bin1] += (*itC2).second;
                }
            }

            (*meanS) += (*S)[chrom1][bin1];
            nstat++;
        }
    }

    (*meanS) /= nstat;

    return;
}

void update_DB(map<string, map<int, double> > *DB, map<string, map<int, double> > &S, double meanS)
{
    for (map<string, map<int, double> >::iterator itS = S.begin(); itS != S.end(); itS++) {
        string chrom = (*itS).first;
        cout << "Updating DB chrom " << chrom << "          \r";
        cout.flush();

        for (map<int, double>::iterator itC = S[chrom].begin(); itC != S[chrom].end(); itC++) {
            int bin = (*itC).first;
            (*DB)[chrom][bin] = S[chrom][bin] / meanS;
            // cout << S[chrom][bin] <<"\t" << meanS << endl;
        }
    }
    return;
}

void update_W(map<string, map<int, map<string, map<int, double> > > > *W, map<string, map<int, double> > &DB)
{
    for (map<string, map<int, map<string, map<int, double> > > >::iterator itW = (*W).begin(); itW != (*W).end(); itW++) {
        string chrom1 = (*itW).first;
        cout << "Updating W chrom " << chrom1 << "          \r";
        cout.flush();

        for (map<int, map<string, map<int, double> > >::iterator itC1 = (*W)[chrom1].begin(); itC1 != (*W)[chrom1].end(); itC1++) {
            int bin1 = (*itC1).first;

            for (map<string, map<int, double> >::iterator itCb1 = (*W)[chrom1][bin1].begin(); itCb1 != (*W)[chrom1][bin1].end(); itCb1++) {
                string chrom2 = (*itCb1).first;
                for (map<int, double>::iterator itC2 = (*W)[chrom1][bin1][chrom2].begin(); itC2 != (*W)[chrom1][bin1][chrom2].end(); itC2++) {
                    int bin2 = (*itC2).first;

                    (*W)[chrom1][bin1][chrom2][bin2] = (*W)[chrom1][bin1][chrom2][bin2] / DB[chrom1][bin1] / DB[chrom2][bin2];
                }
            }
        }
    }

    return;
}

void update_B(map<string, map<int, double> > *B, map<string, map<int, double> > &DB)
{
    for (map<string, map<int, double> >::iterator itB = (*B).begin(); itB != (*B).end(); itB++) {
        string chrom = (*itB).first;
        cout << "Updating DB chrom " << chrom << "          \r";
        cout.flush();

        for (map<int, double>::iterator itC = (*B)[chrom].begin(); itC != (*B)[chrom].end(); itC++) {
            int bin = (*itC).first;
            (*B)[chrom][bin] *= DB[chrom][bin];
        }
    }
    return;
}

void print_out_W(map<string, map<int, map<string, map<int, double> > > > &W, double *mtW,
        map<string, map<int, map<string, map<int, bool> > > > &pairW, Files *file, int nloop)
{
    stringstream ss_out;
    ss_out << file->HiC_iterative_fo << "/W_iteration" << nloop << ".txt";

    stringstream ss_norm_out;
    ss_norm_out << file->HiC_iterative_fo << "/Wnorm_iteration" << nloop << ".txt";
    
    ofstream ofstr(ss_out.str().c_str());
    ofstr << "Chrom1\tBin1" << "\t" << "Chrom2\tBin2" << "\t" << "HiC_bin" << endl;
    ofstream ofstr_norm(ss_norm_out.str().c_str());
    ofstr_norm << "Chrom1\tBin1" << "\t" << "Chrom2\tBin2" << "\t" << "HiC_bin" << endl;
    
    // mean total weight
    *mtW = 0.; 
    int nmtW = 0;
    for (map < string, map<int, map < string, map<int, bool> > > >::iterator itW = pairW.begin(); itW != pairW.end(); itW++) {
        string chrom1 = (*itW).first;

        cout << "Printing out intermediate W chrom " << chrom1 << "          \r";
        cout.flush();
                
        for (map<int, map < string, map<int, bool> > >::iterator itC1 = pairW[chrom1].begin(); itC1 != pairW[chrom1].end(); itC1++) {
            int bin1 = (*itC1).first;

            double norm = 0.;
            for (map < string, map<int, bool> >::iterator itCb1 = pairW[chrom1][bin1].begin(); itCb1 != pairW[chrom1][bin1].end(); itCb1++) {
                string chrom2 = (*itCb1).first;                                
                for (map<int, bool>::iterator itC2 = pairW[chrom1][bin1][chrom2].begin(); itC2 != pairW[chrom1][bin1][chrom2].end(); itC2++) {
                    int bin2 = (*itC2).first;
                    norm += W[chrom1][bin1][chrom2][bin2];
                }                                
            }
            *mtW += norm;
            nmtW++;
            
            for (map < string, map<int, bool> >::iterator itCb1 = pairW[chrom1][bin1].begin(); itCb1 != pairW[chrom1][bin1].end(); itCb1++) {
                string chrom2 = (*itCb1).first;
                
                for (map<int, bool>::iterator itC2 = pairW[chrom1][bin1][chrom2].begin(); itC2 != pairW[chrom1][bin1][chrom2].end(); itC2++) {
                    int bin2 = (*itC2).first;
                    ofstr << chrom1 << "\t" << bin1 << "\t" << chrom2 << "\t" << bin2 << "\t" << W[chrom1][bin1][chrom2][bin2] << endl;
                    ofstr_norm << chrom1 << "\t" << bin1 << "\t" << chrom2 << "\t" << bin2 << "\t" << W[chrom1][bin1][chrom2][bin2] / norm << endl;
                }
            }
        }
    }
     *mtW /= nmtW;
    
    ofstr.close();
    ofstr_norm.close();
    
    return;
}

void print_out_B(map<string, map<int, double> > &B, double mtW, Files *file, int nloop)
{
    // rem: we compute the normalization again in order to avoid ob of memory; can be optimized
    stringstream ss_out;
    ss_out << file->HiC_iterative_fo << "/B_iteration" << nloop << ".txt";
    
    ofstream ofstr(ss_out.str().c_str());
    ofstr << "Chrom\tBin\tB"<< endl;
    
    for (map<string, map<int, double> >::iterator itB = B.begin(); itB != B.end(); itB++) {
        string chrom = (*itB).first;
        cout << "Printing out intermediate B values " << chrom << "          \r";
        cout.flush();

        for (map<int, double>::iterator itC = B[chrom].begin(); itC != B[chrom].end(); itC++) {
            int bin = (*itC).first;
            ofstr << chrom << "\t" << bin << "\t" <<  B[chrom][bin]*sqrt(mtW) << endl;            
        }
    }    
    ofstr.close();
    return;
}
