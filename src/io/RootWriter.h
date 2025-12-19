#pragma once

#include <string>
#include "TFile.h"
#include "TTree.h"

class RootWriter {
public:
    RootWriter(const std::string &fname);
    ~RootWriter();

    void initialize();
    void fill_event(int evtid, double x, double Q2, double t,
                    double weight, double pdf_k, double gpd);
    void write();

private:
    std::string fname_;
    TFile *tfile_;
    TTree *tree_;

    int evtid_;
    double x_;
    double Q2_;
    double t_;
    double weight_;
    double pdf_k_;
    double gpd_;
};

