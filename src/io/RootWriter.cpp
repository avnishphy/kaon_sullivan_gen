#include "RootWriter.h"
#include <iostream>

RootWriter::RootWriter(const std::string &fname)
    : fname_(fname), tfile_(nullptr), tree_(nullptr)
{}

RootWriter::~RootWriter()
{
    if(tfile_) tfile_->Close();
}

void RootWriter::initialize()
{
    tfile_ = TFile::Open(fname_.c_str(), "RECREATE");
    if (!tfile_ || tfile_->IsZombie()) {
        std::cerr << "ERROR: Could not open ROOT file " << fname_ << "\n";
        return;
    }

    tree_ = new TTree("events", "kaon sullivan events");

    tree_->Branch("evtid",    &evtid_,   "evtid/I");
    tree_->Branch("x",        &x_,       "x/D");
    tree_->Branch("Q2",       &Q2_,      "Q2/D");
    tree_->Branch("t",        &t_,       "t/D");
    tree_->Branch("weight",   &weight_,  "weight/D");
    tree_->Branch("pdf_k",    &pdf_k_,   "pdf_k/D");
    tree_->Branch("gpd",      &gpd_,     "gpd/D");
}

void RootWriter::fill_event(int evtid, double x, double Q2, double t,
                            double weight, double pdf_k, double gpd)
{
    evtid_ = evtid;
    x_ = x;
    Q2_ = Q2;
    t_ = t;
    weight_ = weight;
    pdf_k_ = pdf_k;
    gpd_ = gpd;

    tree_->Fill();
}

void RootWriter::write()
{
    tfile_->Write();
}
