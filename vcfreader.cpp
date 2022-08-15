// #define TEST

#include "vcfreader.h"
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <iostream>
#include <cmath>
//#include "calcGL.h"
#include <vector>

std::vector<std::vector<float>> calcGL(char* vcfFN, char* chrom){
    std::vector<std::vector<float>> gls;
    int32_t *PL = NULL;
    int nPL = 0;
    int plRet = 0;
    long old_pos = -1;

    htsFile *vcf= NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = bcf_init();
    vcf = bcf_open(vcfFN, "r");
    if(vcf == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    hdr = bcf_hdr_read(vcf);
    if(hdr == NULL) {
        throw std::runtime_error("Unable to read header.");
    }


    std::cout << "chromosome\tposition\tnum_alleles" << std::endl;
    while(bcf_read(vcf, hdr, rec) == 0) {

        if(rec == NULL) {
            throw std::runtime_error("Unable to read record.");
        }

        if (strcmp(bcf_hdr_id2name(hdr, rec->rid), chrom) != 0) continue;

        nPL = 0;
        plRet = bcf_get_format_int32(hdr, rec, "PL", &PL, &nPL );
        if (plRet < 0){
            std::cout << plRet << std::endl;
            throw std::runtime_error("Error in accessing vcfs with ");
        }


        std::cout << bcf_hdr_id2name(hdr, rec->rid) << "\t" <<
                  rec->pos << "\t" <<
                  rec->n_allele << "\t" << nPL << "\t" <<
                  PL[0] << "\t" <<
                  std::endl;


        std::vector<float> pl_vec;
        int hom = 0, het = 0;
        for(int i = 0; i < nPL; i++){
            if (i == 0 or i == 2 or i == 5 or i == 9) hom += PL[i];
            else het += PL[i];
        }
        float homReduced = 0, hetReduced = 0;
        homReduced = (float)hom/4;
        hetReduced = (float)het/6;
        pl_vec.push_back(homReduced);
        pl_vec.push_back(hetReduced);

        if (old_pos != -1 and rec->pos != old_pos+1 ){

            for (int k = 0; k < rec->pos - old_pos; k++){
                gls.push_back({255, 0});
            }
        }


        gls.push_back(pl_vec);

        old_pos = rec->pos;

    }
    bcf_hdr_destroy(hdr);
    bcf_destroy(rec);
    bcf_close(vcf);
    return gls;
}





/*
For name of vcf file and chromosome name(???ID) returns struct rawdata of gls scores

PARAMETERS
@ inf: char*
file name string
@ chr: char*
chrom name

RETURNS
@ret: rawdata
structure




*/
rawdata readvcf(char* inf,char* chr){//infile

    std::vector<std::vector<float>> tmpgls = calcGL(inf, chr);

    rawdata ret;

    
    ret.pos = new int[1];
    ret.len = tmpgls.size();
    ret.gls = new double[ret.len];
//    ret.pos[0] = pos;

    for (int i = 0; i<ret.len; i++) {
        ret.gls[i] = log(0.0);// = -infinity
        if (tmpgls[i][0] != tmpgls[i][1]) {
//            fprintf(stderr,"@tmps: %f,%f\n",tmpgls[i][0], tmpgls[i][1]);
            double mmax = std::max(tmpgls[i][0], tmpgls[i][1]);
            double val =
                    std::min(tmpgls[i][0], tmpgls[i][1]) - mmax; //??? Why fill it min - max ANS:PREVENT UNDERFLOW

            ret.gls[i] = val;
            if (tmpgls[i][0] < tmpgls[i][1])
                ret.gls[i] = -ret.gls[i];//????Why reverse

//        }

//            fprintf(stderr, "@ret.gls[%d] = %f\n", i, ret.gls[i]);

        }
    }
    ret.lastp = 0;
    return ret;

}

// int main(){
//     rawdata check;
//     check = readvcf("NA12878.chr22.vcf.gz","22");
// }
