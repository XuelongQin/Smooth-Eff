void smooth_th3(TH3 *h){
// Sample kernel
const Int_t ksize_x=5;
const Int_t ksize_y=5;
const Int_t ksize_z=5;
Double_t kernel[ksize_x][ksize_y][ksize_z] = {
{
	{0,0,0,0,0},
	{0,0,0,0,0},
	{0,0,1,0,0},
	{0,0,0,0,0},
	{0,0,0,0,0}
},
{
	{0,0,0,0,0},
	{0,0,1,0,0},
	{0,2,2,2,0},
	{0,0,1,0,0},
	{0,0,0,0,0}
},
{
	{0,0,1,0,0},
	{0,2,2,2,0},
	{1,2,5,2,1},
	{0,2,2,2,0},
	{0,0,1,0,0}
},
{
   {0,0,0,0,0},
   {0,0,1,0,0},
   {0,2,2,2,0},
   {0,0,1,0,0},
   {0,0,0,0,0}
},
{
   {0,0,0,0,0},
   {0,0,0,0,0},
   {0,0,1,0,0},
   {0,0,0,0,0},
   {0,0,0,0,0}
}
};

// Determine the size of the bin buffer(s) needed
Int_t lowest_bin=h->GetBin(0,0,0);
Int_t highest_bin=h->GetBin(h->GetNbinsX()+1,h->GetNbinsY()+1,h->GetNbinsZ()+1);
Int_t bufSize = highest_bin-lowest_bin+1;
Double_t *buf = new Double_t[bufSize];
Double_t *ebuf = new Double_t[bufSize];

// Copy all the data to the temporry buffers
for (Int_t i=0; i<=h->GetNbinsX(); i++){
for (Int_t j=0; j<=h->GetNbinsY(); j++){
for (Int_t k=0; k<=h->GetNbinsZ(); k++){
Int_t bin = h->GetBin(i,j,k);
buf[bin] =h->GetBinContent(bin);
ebuf[bin]=h->GetBinError(bin);
};
};
};

// Kernel tail sizes (kernel sizes must be odd for this to work!)
Int_t x_push = (ksize_x-1)/2;
Int_t y_push = (ksize_y-1)/2;
Int_t z_push = (ksize_z-1)/2;

// main work loop
for (Int_t i=1; i<=h->GetNbinsX(); i++){
for (Int_t j=1; j<=h->GetNbinsY(); j++){
for (Int_t l=1; l<=h->GetNbinsZ(); l++){
Double_t content = 0.0;
Double_t error = 0.0;
Double_t norm = 0.0;

  for (Int_t n=0; n<ksize_x; n++){
    for (Int_t m=0; m<ksize_y; m++){
		 for (Int_t o=0; o<ksize_z; o++){
      Int_t xb = i+(n-x_push);
      Int_t yb = j+(m-y_push);
		Int_t zb = l+(o-z_push);
      if ( (xb >= 1) && (xb <= h->GetNbinsX()) &&
           (yb >= 1) && (yb <= h->GetNbinsY()) &&
			  (zb >= 1) && (zb <= h->GetNbinsZ())){
        Int_t bin = h->GetBin(xb,yb,zb);
        Double_t k = kernel[n][m][o];
        if ( (k != 0.0 ) && (buf[bin] != 0.0) ) { // General version probably does not want the second condition
          content += k*buf[bin];
          error   += k*k*buf[bin]*buf[bin];
          norm    += k;
        };
      };
		  };
  };
  };

  content /= norm;
  error /= (norm*norm);
  if ( content != 0.0 ) { // General version probably does not want this condition
    h->SetBinContent(i,j,l,content);
    h->SetBinError(i,j,l,sqrt(error));
  };
};
};
};
delete[] buf;
delete[] ebuf;
};


void finalsmooth(){
	for (int i=0;i<8;i++){
		if (i==4||i==6)
			continue;
	TFile *fin = new TFile(Form(("/Users/xuelong/workdir/B0KstMuMu/eff/files_forSimFit/KDEeff_b%i_od_2017.root"),i));
	TH3D *we = (TH3D*)fin->Get(Form(("effWHist_b%ip1"),i));
	smooth_th3(we);
	TH3D *ce = (TH3D*)fin->Get(Form(("effCHist_b%ip1"),i));
	TH1D *mc = (TH1D*)fin->Get(Form(("MCint_b%ip1t1"),i));
	TH1D *mw = (TH1D*)fin->Get(Form(("MCint_b%ip1t0"),i));

	TFile *fout = new TFile(Form(("smeff/KDEeff_b%i_od_2017.root"),i),"recreate");
	we->Write();
	ce->Write();
	mc->Write();
	mw->Write();
	}
}
