#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "fastas.h"
#include "header.h"
#include "psmcreader.h"
#include <vector>

rawdata readvcf(char* filename, char* chr);
std::vector<std::vector<float>> calcGL(char* vcfFN, char* chrom);