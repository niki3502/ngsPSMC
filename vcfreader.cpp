#include "vcfreader.h"
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <iostream>
rawdata readvcf(char* inf,char* chr){//infile

    int *pl = NULL;
    int npl_arr = 0;
    int npl = 0;

    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    int nad_arr = 0;
    int nad     = 0;
    int *ad     = NULL;

    rawdata ret;
    htsFile *fp = bcf_open(inf, "r");
    bcf1_t *rec = bcf_init();
    int *tmpgls;
    if(fp == NULL) {
        throw std::runtime_error("Unable to open file.");
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);

    if(hdr == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
    int rid = bcf_hdr_name2id(hdr,chr);
    bcf_idpair_t *ctg =hdr->id[BCF_DT_CTG];
    printf("file %s is successfully read\n",inf);
    int length;
    printf("Information from header:\n");
    for (int i = 0; i < hdr->n[BCF_DT_CTG]; ++i){
        printf("%s\t%d\n", ctg[i].key, ctg[i].val->info[0]);
        if (strcmp(ctg[i].key,chr)==0){
            length = ctg[i].val->info[0];
        }



    }

    

    tmpgls = (int*)calloc(length,sizeof(int));

    int k = 0;
    int pos = -1;
    while(bcf_read(fp, hdr, rec) == 0){
        if (rec->rid!=rid)continue;
        if (pos == -1)pos = k;
        bcf_unpack(rec, BCF_UN_ALL);
        npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
        tmpgls[k] = pl[0]; 
        k+=1;      
    }


    ret.pos = new int[1];
    ret.len = k+1; 
    ret.gls = new double[k+1];
    ret.pos[0] = pos;
    for (int i = 0; i<k+1;i++)
    ret.gls[i] = (mygltype)tmpgls[i];


    return ret;

}
